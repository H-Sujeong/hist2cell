"""WSI -> Hist2Cell inference Data, v2 (framework-based tissue masking).

Differences from v1:
  - tissue mask uses WSI_tile_sampling_framework.ForegroundMasker (multi-channel
    YUV/HSV + paraffin/background filters), much better than v1's HSV+Otsu on
    slides with stickers, fat, or uneven staining
  - thumbnail comes straight from openslide.get_thumbnail (no matplotlib
    savefig path that silently rescales the saved PNG)
  - tile sampling = simple stride loop with min_tissue_frac, no contour-based
    post-filter (the framework's TileSampler drops sparse off-center tissue
    regions; we want them in)
  - graph = kNN on tile centers (k=6 + self-loop) — preserves the ~6-neighbor
    structure the GAT was trained on, regardless of tile pattern

Outputs in --output:
  <slide>.pt              PyG Data(x, edge_index, pos, spot_id) for the model
  <slide>_coords.h5       full-res tile coords + metadata (reference format)
  spots.csv               spot_id, X, Y, tile_x_topleft, tile_y_topleft
  tissue_mask.png         tissue mask at thumbnail resolution
  spot_view.jpg           thumbnail with red tile rectangles at TRUE scale

Usage:
  python prep/prepare_wsi_for_inference_v2.py \\
      --input  /path/to/slide.svs \\
      --output ./inference/MY_SLIDE \\
      --tile-size 400
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import openslide
import pandas as pd
import torch
from PIL import Image, ImageDraw
from torch_geometric.data import Data
from torchvision import transforms

# Resolve the framework relative to repo root so this works from anywhere.
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "WSI_tile_sampling_framework"))
from ForegroundMasking import ForegroundMasker  # noqa: E402

Image.MAX_IMAGE_PIXELS = None

IMAGENET_TRANSFORM = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406],
                         std=[0.229, 0.224, 0.225]),
])


# ---------- thumbnail + mask ----------

def make_thumbnail(slide, max_side):
    """openslide thumbnail straight to PIL — no matplotlib savefig roundtrip."""
    W, H = slide.dimensions
    scale = max(1.0, max(W, H) / float(max_side))
    tw, th = int(round(W / scale)), int(round(H / scale))
    thumb = slide.get_thumbnail((tw, th)).convert("RGB")
    return thumb, scale


def tissue_mask_via_framework(thumb_pil):
    """Multi-channel tissue mask via framework's ForegroundMasker.
    Returns uint8 0/1 array shaped like thumb (H, W)."""
    masker = ForegroundMasker()
    _, mask = masker.get_foreground(np.asarray(thumb_pil), is_norm=False)
    return (mask > 0).astype(np.uint8)


# ---------- tile coord generation ----------

def sample_tile_coords(wsi_dims, mask, scale, tile_size, min_tissue_frac):
    """Iterate full-res grid with stride=tile_size; keep tiles whose mapped
    tissue mask region has >= min_tissue_frac coverage. No contour post-filter.

    wsi_dims : (W, H) at level 0
    mask     : 0/1 array at thumbnail resolution
    scale    : full-res / thumbnail px ratio (>=1)
    """
    W, H = wsi_dims
    mh, mw = mask.shape
    coords = []
    for y in range(0, H - tile_size + 1, tile_size):
        ty0 = int(round(y / scale))
        ty1 = int(round((y + tile_size) / scale))
        ty0 = max(0, min(mh, ty0))
        ty1 = max(0, min(mh, ty1))
        if ty1 <= ty0:
            continue
        for x in range(0, W - tile_size + 1, tile_size):
            tx0 = int(round(x / scale))
            tx1 = int(round((x + tile_size) / scale))
            tx0 = max(0, min(mw, tx0))
            tx1 = max(0, min(mw, tx1))
            if tx1 <= tx0:
                continue
            sub = mask[ty0:ty1, tx0:tx1]
            if sub.size == 0:
                continue
            if float(sub.mean()) >= min_tissue_frac:
                coords.append((x, y))
    return np.asarray(coords, dtype=np.int64) if coords else np.zeros((0, 2), np.int64)


# ---------- kNN graph ----------

def build_knn_edges(centers_xy, k=6):
    """Mutual kNN graph (union) + self-loops, returned as edge_index (2, E)."""
    n = len(centers_xy)
    P = centers_xy.astype(np.float64)
    edges = []
    chunk = 1024
    for i0 in range(0, n, chunk):
        i1 = min(n, i0 + chunk)
        d2 = ((P[i0:i1, None, :] - P[None, :, :]) ** 2).sum(-1)  # (chunk, n)
        kk = min(k + 1, n)
        nn = np.argpartition(d2, kk - 1, axis=1)[:, :kk]
        for r, row in enumerate(nn):
            i = i0 + r
            for j in row:
                if int(j) != i:
                    edges.append((i, int(j)))
    if not edges:
        edges = np.zeros((0, 2), np.int64)
    else:
        edges = np.asarray(edges, np.int64)
    self_loops = np.stack([np.arange(n), np.arange(n)], axis=1)
    if len(edges):
        edges = np.unique(np.concatenate([edges, edges[:, [1, 0]], self_loops], axis=0), axis=0)
    else:
        edges = self_loops
    return torch.from_numpy(np.ascontiguousarray(edges.T)).long()


# ---------- main ----------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="path to WSI (.svs etc)")
    ap.add_argument("--output", required=True, help="output directory")
    ap.add_argument("--slide-name", default=None,
                    help="defaults to input stem with spaces/commas replaced")
    ap.add_argument("--tile-size", type=int, default=400,
                    help="grid spacing in level-0 px; pick to match Visium ~100um at this WSI's mpp")
    ap.add_argument("--patch-size", type=int, default=224,
                    help="model input crop centered on each tile (must be 224 for Hist2Cell)")
    ap.add_argument("--min-tissue-frac", type=float, default=0.10,
                    help="min tissue fraction in a tile region; 0.10 is reference default")
    ap.add_argument("--thumb-max-side", type=int, default=4000,
                    help="cap on thumbnail max side; larger = finer mask, slower mask compute")
    ap.add_argument("--knn", type=int, default=6,
                    help="k for the spatial graph (matches Visium hex 6-neighbor + self-loop)")
    ap.add_argument("--save-patches", action="store_true",
                    help="also dump individual jpg patches to <output>/patches/")
    return ap.parse_args()


def main():
    args = parse_args()
    in_path = Path(args.input)
    if not in_path.is_file():
        raise FileNotFoundError(in_path)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.slide_name:
        slide_name = args.slide_name
    else:
        slide_name = in_path.stem.replace(" ", "_").replace(",", "_")

    print(f"[1/6] open WSI: {in_path}")
    slide = openslide.OpenSlide(str(in_path))
    W, H = slide.dimensions
    mpp_x = slide.properties.get("openslide.mpp-x")
    mpp_y = slide.properties.get("openslide.mpp-y")
    obj = slide.properties.get("openslide.objective-power")
    print(f"     dims={W}x{H}  mpp=({mpp_x},{mpp_y})  obj={obj}")

    print(f"[2/6] thumbnail + tissue mask (framework ForegroundMasker)")
    thumb, scale = make_thumbnail(slide, args.thumb_max_side)
    print(f"     thumb size = {thumb.size} (scale = {scale:.2f})")
    mask = tissue_mask_via_framework(thumb)
    print(f"     tissue fraction on thumb = {mask.mean():.3f}")
    Image.fromarray(mask * 255).save(out_dir / "tissue_mask.png")

    print(f"[3/6] sample tiles (tile_size={args.tile_size}, min_frac={args.min_tissue_frac})")
    coords = sample_tile_coords((W, H), mask, scale, args.tile_size, args.min_tissue_frac)
    n = len(coords)
    if n == 0:
        raise SystemExit("No tile coords. Lower --min-tissue-frac or check the mask.")
    print(f"     kept {n} tiles")

    h5_path = out_dir / f"{slide_name}_coords.h5"
    with h5py.File(h5_path, "w") as f:
        f.create_dataset("coords", data=coords, compression="gzip")
        m = f.create_group("metadata")
        m.attrs["source_slide"] = str(in_path)
        m.attrs["slide_name"] = slide_name
        m.attrs["tile_size"] = int(args.tile_size)
        m.attrs["patch_size"] = int(args.patch_size)
        m.attrs["overlap"] = 0
        m.attrs["stride"] = int(args.tile_size)
        m.attrs["min_tissue_fraction"] = float(args.min_tissue_frac)
        m.attrs["total_tiles"] = int(n)
        m.attrs["coord_format"] = "x_y_level0_topleft"
        m.attrs["wsi_width"] = int(W)
        m.attrs["wsi_height"] = int(H)
        if mpp_x: m.attrs["mpp_x"] = str(mpp_x)
        if mpp_y: m.attrs["mpp_y"] = str(mpp_y)
        if obj:   m.attrs["objective_power"] = str(obj)
    print(f"     wrote {h5_path}")

    # Overlay at TRUE thumbnail scale (no matplotlib).
    sx = thumb.size[0] / W
    sy = thumb.size[1] / H
    overlay = thumb.copy()
    draw = ImageDraw.Draw(overlay)
    for x, y in coords:
        x0, y0 = x * sx, y * sy
        x1, y1 = (x + args.tile_size) * sx, (y + args.tile_size) * sy
        draw.rectangle([x0, y0, x1, y1], outline="red", width=1)
    overlay.save(out_dir / "spot_view.jpg", quality=85)

    print(f"[4/6] extract {args.patch_size}x{args.patch_size} patches at level 0")
    half = args.patch_size // 2
    centers = coords + np.array([args.tile_size // 2, args.tile_size // 2], dtype=np.int64)
    x_tensor = torch.zeros((n, 3, args.patch_size, args.patch_size), dtype=torch.float32)
    pos = np.zeros((n, 2), dtype=np.float32)
    spot_ids = []
    if args.save_patches:
        (out_dir / "patches").mkdir(exist_ok=True)
    for i, (cx, cy) in enumerate(centers):
        x0, y0 = int(cx - half), int(cy - half)
        patch = slide.read_region((x0, y0), 0, (args.patch_size, args.patch_size)).convert("RGB")
        if patch.size != (args.patch_size, args.patch_size):
            patch = patch.resize((args.patch_size, args.patch_size), Image.BILINEAR)
        sid = f"{slide_name}_x{int(cx)}y{int(cy)}"
        spot_ids.append(sid)
        pos[i] = (float(cx), float(cy))
        x_tensor[i] = IMAGENET_TRANSFORM(patch)
        if args.save_patches:
            patch.save(out_dir / "patches" / f"{sid}.jpg", quality=90)
        if (i + 1) % 1000 == 0 or (i + 1) == n:
            print(f"     {i+1}/{n}")

    print(f"[5/6] build kNN graph (k={args.knn})")
    edge_index = build_knn_edges(centers.astype(np.float64), k=args.knn)
    print(f"     edges = {edge_index.shape[1]}  (~{edge_index.shape[1]/n:.1f} per node incl. self-loop & symmetrization)")

    print(f"[6/6] save PyG Data")
    data = Data(x=x_tensor, edge_index=edge_index,
                pos=torch.from_numpy(pos), spot_id=spot_ids)
    pt_path = out_dir / f"{slide_name}.pt"
    torch.save(data, pt_path)

    df = pd.DataFrame({
        "spot_id": spot_ids,
        "X": pos[:, 0].astype(int),
        "Y": pos[:, 1].astype(int),
        "tile_x_topleft": coords[:, 0],
        "tile_y_topleft": coords[:, 1],
    })
    df.to_csv(out_dir / "spots.csv", index=False)

    slide.close()
    print(f"\nDone:")
    print(f"  {pt_path}            ({n} nodes, {edge_index.shape[1]} edges)")
    print(f"  {out_dir/'spots.csv'}")
    print(f"  {h5_path}")
    print(f"  {out_dir/'tissue_mask.png'}, {out_dir/'spot_view.jpg'}")


if __name__ == "__main__":
    main()
