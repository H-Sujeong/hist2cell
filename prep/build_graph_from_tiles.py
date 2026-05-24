"""
Tiling HDF5 (WSI_tile_sampling_framework) -> PyG Data(.pt) graph 빌드.

각 <slide>.h5 (tile 좌표 + metadata) + 원본 SVS 로부터:
  - tile 좌상단에서 patch 추출 (h5 의 tile_size px, level 0)
  - tile_size != 224 이면 224 로 resize (모델 입력 고정 크기)
  - ImageNet 정규화 -> x [N, 3, 224, 224] float32
  - tile 중심으로 kNN(k) 그래프 -> edge_index (symmetric union + self-loop)
  - Data(x, edge_index, pos, spot_id) 를 <slide>.pt 로 저장 + <slide>_spots.csv

112 타일(56µm)은 추출 후 224 로 ×2 업샘플된다 — 40× 학습 모델(HEX)용.
224 타일(112µm)은 resize 없이 그대로 — Hist2Cell(20× 학습)용.

Usage:
  python prep/build_graph_from_tiles.py \
      --tiles-dir lung_pilot/tilitng_output/TCGA-LUAD \
      --wsi-dir   /mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi \
      --output    lung_pilot/graph_output/224
"""
import argparse
import glob
import os
import time

import h5py
import numpy as np
import openslide
import pandas as pd
import torch
from PIL import Image
from sklearn.neighbors import NearestNeighbors
from torch_geometric.data import Data
from torchvision import transforms

MODEL_INPUT = 224  # ResNet18 / 대부분 backbone 의 고정 입력 크기
IMAGENET = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406],
                         std=[0.229, 0.224, 0.225]),
])


def build_knn_edges(centers, k=6):
    """tile 중심에 대한 kNN 그래프. symmetric union + self-loop, edge_index (2,E)."""
    n = len(centers)
    kk = min(k + 1, n)  # self 포함
    _, idx = NearestNeighbors(n_neighbors=kk).fit(centers).kneighbors(centers)
    src = np.repeat(np.arange(n), kk)
    dst = idx.reshape(-1)
    e = np.stack([src, dst], axis=0)
    e = np.concatenate([e, e[::-1]], axis=1)                       # 양방향
    e = np.concatenate([e, np.stack([np.arange(n)] * 2)], axis=1)  # self-loop
    e = np.unique(e, axis=1)
    return torch.from_numpy(np.ascontiguousarray(e)).long()


def find_svs(wsi_dir, slide):
    for pat in (f"{slide}.*.svs", f"{slide}.svs"):
        hits = sorted(glob.glob(os.path.join(wsi_dir, pat)))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"{slide}: SVS 를 {wsi_dir} 에서 못 찾음")


def process(h5_path, wsi_dir, out_dir, knn):
    slide = os.path.splitext(os.path.basename(h5_path))[0]
    with h5py.File(h5_path, "r") as hf:
        coords = hf["coords"][:].astype(np.int64)        # (N,2) 좌상단, level 0
        tile_size = int(hf["metadata"].attrs["tile_size"])
    svs = find_svs(wsi_dir, slide)
    sl = openslide.OpenSlide(svs)
    n = len(coords)
    centers = coords + tile_size // 2
    resize = tile_size != MODEL_INPUT
    print(f"  {slide}: N={n}  tile_size={tile_size}px  "
          f"{'resize->224' if resize else 'no-resize'}  svs={os.path.basename(svs)}")

    x = torch.zeros((n, 3, MODEL_INPUT, MODEL_INPUT), dtype=torch.float32)
    spot_ids = []
    t0 = time.time()
    for i, (tx, ty) in enumerate(coords):
        patch = sl.read_region((int(tx), int(ty)), 0,
                               (tile_size, tile_size)).convert("RGB")
        if patch.size != (MODEL_INPUT, MODEL_INPUT):
            patch = patch.resize((MODEL_INPUT, MODEL_INPUT), Image.BILINEAR)
        x[i] = IMAGENET(patch)
        spot_ids.append(f"{slide}_x{int(centers[i, 0])}y{int(centers[i, 1])}")
        if (i + 1) % 4000 == 0 or i + 1 == n:
            print(f"    {slide}: {i+1}/{n}  ({time.time()-t0:.0f}s)", flush=True)
    sl.close()

    edge_index = build_knn_edges(centers.astype(np.float64), k=knn)
    data = Data(x=x,
                edge_index=edge_index,
                pos=torch.from_numpy(centers.astype(np.float32)),
                spot_id=spot_ids)
    pt_path = os.path.join(out_dir, f"{slide}.pt")
    torch.save(data, pt_path)
    pd.DataFrame({
        "spot_id": spot_ids,
        "X": centers[:, 0], "Y": centers[:, 1],
        "tile_x_topleft": coords[:, 0], "tile_y_topleft": coords[:, 1],
        "tile_size": tile_size,
    }).to_csv(os.path.join(out_dir, f"{slide}_spots.csv"), index=False)

    gb = os.path.getsize(pt_path) / 1e9
    print(f"  -> {pt_path}  ({n} nodes, {edge_index.shape[1]} edges, {gb:.1f} GB)")
    return slide, n, tile_size, edge_index.shape[1], gb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles-dir", required=True, help="framework <slide>.h5 들이 있는 폴더")
    ap.add_argument("--wsi-dir", required=True, help="원본 .svs 폴더")
    ap.add_argument("--output", required=True, help="출력 폴더")
    ap.add_argument("--knn", type=int, default=6, help="kNN 그래프 k (Visium hex 6-이웃)")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)
    h5s = sorted(glob.glob(os.path.join(args.tiles_dir, "*.h5")))
    if not h5s:
        raise SystemExit(f"h5 없음: {args.tiles_dir}")
    print(f"{len(h5s)} h5 -> {args.output}  (knn={args.knn})")

    rows = []
    for h5 in h5s:
        print(f"\n== {os.path.basename(h5)}")
        rows.append(process(h5, args.wsi_dir, args.output, args.knn))

    print(f"\n=== SUMMARY ({args.output}) ===")
    for slide, n, ts, e, gb in rows:
        print(f"  {slide:26s}  nodes={n:6d}  tile={ts}px  edges={e:8d}  {gb:5.1f} GB")


if __name__ == "__main__":
    main()
