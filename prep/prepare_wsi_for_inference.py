"""Prepare a histology WSI for Hist2Cell inference.

Given a WSI path, this script:
  1. detects tissue regions on a downsampled thumbnail
  2. lays a hexagonal (Visium-like) grid of spot centers across the slide
  3. keeps spots whose 224x224 patch has enough tissue
  4. extracts patches at full resolution and normalizes them with ImageNet stats
  5. builds the spot-to-spot neighborhood graph used by Hist2Cell
  6. saves a torch_geometric Data(.pt) ready to feed into the trained model,
     plus spots.csv and QC overlays.

The output .pt has the same fields the trained model expects (x, edge_index,
pos, spot_id) but no labels y, since we are predicting them.

Usage:
    python prep/prepare_wsi_for_inference.py \
        --input /path/to/slide.svs \
        --output ./prep_out/MY_SLIDE
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from PIL import Image, ImageDraw
from scipy.ndimage import binary_closing, binary_fill_holes, binary_opening
from torch_geometric.data import Data
from torch_geometric.utils import dense_to_sparse
from torchvision import transforms

Image.MAX_IMAGE_PIXELS = None

try:
    import openslide
    HAS_OPENSLIDE = True
except Exception:
    HAS_OPENSLIDE = False


IMAGENET_TRANSFORM = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406],
                         std=[0.229, 0.224, 0.225]),
])


# ---------- WSI I/O ----------

def open_wsi(path):
    """Return (reader, (width, height), kind) where kind in {'openslide','pil'}."""
    ext = Path(path).suffix.lower()
    if HAS_OPENSLIDE and ext in {".svs", ".ndpi", ".mrxs", ".tif", ".tiff", ".vms", ".vmu", ".scn", ".bif"}:
        try:
            slide = openslide.OpenSlide(path)
            return slide, slide.dimensions, "openslide"
        except Exception as e:
            print(f"[warn] openslide failed ({e}); falling back to PIL")
    img = Image.open(path).convert("RGB")
    return img, img.size, "pil"


def get_thumbnail(reader, dims, kind, max_side):
    w, h = dims
    scale = max(1.0, max(w, h) / float(max_side))
    new_w, new_h = int(round(w / scale)), int(round(h / scale))
    if kind == "openslide":
        thumb = reader.get_thumbnail((new_w, new_h)).convert("RGB")
    else:
        thumb = reader.resize((new_w, new_h), Image.BILINEAR)
    return thumb, scale


def crop_patch(reader, kind, X, Y, patch_size):
    half = patch_size // 2
    x0, y0 = int(X - half), int(Y - half)
    if kind == "openslide":
        patch = reader.read_region((x0, y0), 0, (patch_size, patch_size)).convert("RGB")
    else:
        patch = reader.crop((x0, y0, x0 + patch_size, y0 + patch_size))
    if patch.size != (patch_size, patch_size):
        patch = patch.resize((patch_size, patch_size), Image.BILINEAR)
    return patch


# ---------- tissue detection ----------

def otsu_threshold(values_u8):
    hist, _ = np.histogram(values_u8.ravel(), bins=256, range=(0, 256))
    total = values_u8.size
    if total == 0:
        return 0
    sum_total = (np.arange(256) * hist).sum()
    w_b = 0.0
    sum_b = 0.0
    var_max = 0.0
    threshold = 0
    for i in range(256):
        w_b += hist[i]
        if w_b == 0:
            continue
        w_f = total - w_b
        if w_f == 0:
            break
        sum_b += i * hist[i]
        m_b = sum_b / w_b
        m_f = (sum_total - sum_b) / w_f
        var_between = w_b * w_f * (m_b - m_f) ** 2
        if var_between > var_max:
            var_max = var_between
            threshold = i
    return threshold


def compute_tissue_mask(thumb_pil):
    """HSV-saturation + Otsu mask. Tissue=True, background=False."""
    arr = np.asarray(thumb_pil, dtype=np.uint8)
    mx = arr.max(axis=-1).astype(np.float32)
    mn = arr.min(axis=-1).astype(np.float32)
    sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1.0), 0.0)
    sat_u8 = (sat * 255).astype(np.uint8)
    t = otsu_threshold(sat_u8)
    mask = (sat_u8 > max(t, 8)) & (mx < 230)
    mask = binary_closing(mask, iterations=2)
    mask = binary_fill_holes(mask).astype(bool)
    mask = binary_opening(mask, iterations=2)
    return mask


# ---------- hex grid ----------

def hex_grid_spots(width, height, spot_distance, patch_size):
    """Generate Visium-like hex grid of spot centers entirely inside the slide.

    Returns list of (array_col, array_row, X_pixel, Y_pixel).
    Even rows -> array_col in {0,2,4,...}; odd rows -> {1,3,5,...}.
    Using the same |Δcol|<3 and |Δrow|<2 neighborhood rule as the training data
    therefore yields ~6 neighbors per interior spot.
    """
    dx = float(spot_distance)
    dy = dx * np.sqrt(3) / 2.0
    half = patch_size / 2.0
    spots = []
    row = 0
    while True:
        y = half + row * dy
        if y + half > height:
            break
        offset = (row % 2) * (dx / 2.0)
        col_start_array = row % 2
        col_idx = 0
        while True:
            x = half + offset + col_idx * dx
            if x + half > width:
                break
            array_col = col_start_array + col_idx * 2
            array_row = row
            spots.append((array_col, array_row, x, y))
            col_idx += 1
        row += 1
    return spots


def filter_by_tissue(spots, mask, scale, patch_size, min_frac):
    """Keep spots whose patch footprint on the thumbnail mask has >= min_frac tissue."""
    mh, mw = mask.shape
    half = patch_size / 2.0 / scale
    kept = []
    for ac, ar, X, Y in spots:
        cx, cy = X / scale, Y / scale
        x0 = max(0, int(cx - half))
        x1 = min(mw, int(cx + half))
        y0 = max(0, int(cy - half))
        y1 = min(mh, int(cy + half))
        if x1 <= x0 or y1 <= y0:
            continue
        sub = mask[y0:y1, x0:x1]
        if sub.size == 0:
            continue
        if sub.mean() >= min_frac:
            kept.append((ac, ar, X, Y))
    return kept


def build_hex_adjacency(spots):
    coords = np.array([(ac, ar) for ac, ar, _, _ in spots], dtype=np.int32)
    dc = np.abs(coords[:, 0:1] - coords[:, 0:1].T)
    dr = np.abs(coords[:, 1:2] - coords[:, 1:2].T)
    adj = ((dc < 3) & (dr < 2)).astype(np.float32)
    return adj


# ---------- main ----------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="path to WSI (jpg/png/tif/svs/...)")
    ap.add_argument("--output", required=True, help="output directory")
    ap.add_argument("--slide-name", default=None,
                    help="slide identifier; defaults to input file stem")
    ap.add_argument("--patch-size", type=int, default=224,
                    help="crop side in pixels; must be 224 to match Hist2Cell weights")
    ap.add_argument("--spot-distance", type=float, default=200.0,
                    help="center-to-center spacing of hex spots in WSI pixels. "
                         "Pick this so it matches Visium spacing (~150um) at your WSI's mpp.")
    ap.add_argument("--min-tissue-frac", type=float, default=0.30,
                    help="drop spots whose patch is less than this fraction tissue")
    ap.add_argument("--thumb-max-side", type=int, default=4000,
                    help="max side length of the thumbnail used for tissue detection")
    ap.add_argument("--save-patches", action="store_true",
                    help="also dump individual patch jpgs to <output>/patches/")
    return ap.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.input)
    if not in_path.is_file():
        raise FileNotFoundError(in_path)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    slide = args.slide_name or in_path.stem
    patch_size = args.patch_size

    print(f"[1/6] opening WSI: {in_path}")
    reader, (W, H), kind = open_wsi(str(in_path))
    print(f"      reader={kind}, dims={W}x{H}")

    print(f"[2/6] computing tissue mask")
    thumb, scale = get_thumbnail(reader, (W, H), kind, args.thumb_max_side)
    mask = compute_tissue_mask(thumb)
    print(f"      thumb={thumb.size}, scale={scale:.2f}, tissue_frac={mask.mean():.3f}")

    print(f"[3/6] generating hex grid (spot_distance={args.spot_distance})")
    candidates = hex_grid_spots(W, H, args.spot_distance, patch_size)
    print(f"      candidate spots: {len(candidates)}")

    print(f"[4/6] filtering by tissue (min_frac={args.min_tissue_frac})")
    spots = filter_by_tissue(candidates, mask, scale, patch_size, args.min_tissue_frac)
    if len(spots) == 0:
        raise RuntimeError("No spots passed tissue filter. "
                           "Check --min-tissue-frac and --spot-distance.")
    print(f"      kept: {len(spots)}")

    print(f"[5/6] extracting patches and normalizing")
    n = len(spots)
    x_tensor = torch.zeros((n, 3, patch_size, patch_size), dtype=torch.float32)
    pos = np.zeros((n, 2), dtype=np.float32)
    spot_ids = []
    patch_dir = out_dir / "patches"
    if args.save_patches:
        patch_dir.mkdir(parents=True, exist_ok=True)

    for i, (ac, ar, X, Y) in enumerate(spots):
        sid = f"{slide}_r{ar}c{ac}"
        spot_ids.append(sid)
        pos[i] = (X, Y)
        patch = crop_patch(reader, kind, X, Y, patch_size)
        x_tensor[i] = IMAGENET_TRANSFORM(patch)
        if args.save_patches:
            patch.save(patch_dir / f"{sid}.jpg", quality=90)
        if (i + 1) % 200 == 0 or (i + 1) == n:
            print(f"      {i+1}/{n}")

    print(f"[6/6] building graph and writing outputs")
    adj = build_hex_adjacency(spots)
    edge_index, _ = dense_to_sparse(torch.from_numpy(adj))
    pos_t = torch.from_numpy(pos)

    data = Data(x=x_tensor, edge_index=edge_index, pos=pos_t, spot_id=spot_ids)

    pt_path = out_dir / f"{slide}.pt"
    torch.save(data, pt_path)

    df = pd.DataFrame({
        "spot_id": spot_ids,
        "X": pos[:, 0].astype(int),
        "Y": pos[:, 1].astype(int),
        "array_col": [s[0] for s in spots],
        "array_row": [s[1] for s in spots],
    })
    df.to_csv(out_dir / "spots.csv", index=False)

    Image.fromarray((mask.astype(np.uint8) * 255)).save(out_dir / "tissue_mask.png")

    overlay = thumb.copy()
    draw = ImageDraw.Draw(overlay)
    r_thumb = max(2, int(patch_size / scale / 2))
    for _, _, X, Y in spots:
        cx, cy = X / scale, Y / scale
        draw.ellipse([cx - r_thumb, cy - r_thumb, cx + r_thumb, cy + r_thumb],
                     outline="red", width=1)
    overlay.save(out_dir / "spot_view.jpg", quality=85)

    print(f"\nDone. Wrote:")
    print(f"  {pt_path}                    ({n} nodes, {edge_index.shape[1]} edges)")
    print(f"  {out_dir/'spots.csv'}")
    print(f"  {out_dir/'tissue_mask.png'}, {out_dir/'spot_view.jpg'}")
    if args.save_patches:
        print(f"  {patch_dir}/*.jpg")


if __name__ == "__main__":
    main()
