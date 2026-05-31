"""DINO cluster(=Hist2Cell dominant cell type) 별 centroid-최근접 패치 grid + UMAP overlay.

해상도(224 또는 146) 하나에 대해:
  1) 3 슬라이드 합쳐 dominant cell type = argmax(prediction) 산출.
  2) 4 representation (prediction_log1p / features_fused / features_resnet / features_dinov2)
     cross-slide UMAP → dominant cell type 으로 색칠한 4-panel PNG.
     (DINO morphology cluster 가 Hist2Cell cell-type 과 모이는지/분산되는지 확인)
  3) cluster = dominant cell type. spot 수 >= MIN_GRID 인 cluster 마다:
       centroid = 그 cluster 의 features_dinov2 평균.
       centroid 최근접 min(100, m) 패치를 거리 오름차순(=대표성 순)으로,
       무간격 정방형 grid (좌상단=최근접) PNG 로 저장. 패치는 SVS 재crop.
  4) 패치 추적 CSV: cluster, grid_index(row-major), slide, spot_id, x, y, dist_to_centroid.
     grid_index 순 = 좌상단→우하단.

Usage:
    .venv/bin/python lung_pilot/dino_cluster_patches.py \
        --infer-dir lung_pilot/inference_output \
        --dino-dir  lung_pilot/dino_output \
        --graph-dir lung_pilot/graph_output/224 \
        --wsi-dir   /mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi \
        --out-dir   lung_pilot/dino_cluster_output/224 \
        --res-label 224
"""
import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openslide
import umap
from PIL import Image

SLIDES = [
    "TCGA-05-4245-01A-01-BS1",
    "TCGA-05-4245-01A-01-TS1",
    "TCGA-05-4390-01A-01-BS1",
]
UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)
MIN_GRID = 36          # 정방형 grid 를 그릴 최소 spot 수 (>=6x6)
MAX_PATCH = 100        # centroid 최근접 최대 패치 수 (10x10)
DISP_PX = 128          # grid 내 패치 한 변 픽셀


def short(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def safe(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s)


def find_svs(wsi_dir, slide):
    for pat in (f"{slide}.*.svs", f"{slide}.svs"):
        hits = sorted(glob.glob(os.path.join(wsi_dir, pat)))
        if hits:
            return hits[0]
    raise FileNotFoundError(slide)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infer-dir", required=True)
    ap.add_argument("--dino-dir", required=True)
    ap.add_argument("--graph-dir", required=True, help="<slide>_spots.csv 들이 있는 폴더")
    ap.add_argument("--wsi-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--res-label", required=True)
    a = ap.parse_args()
    INFER, DINO, GRAPH = Path(a.infer_dir), Path(a.dino_dir), Path(a.graph_dir)
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    EMB = OUT / "embeddings"; EMB.mkdir(exist_ok=True)

    # ---- load (3 슬라이드 합치기) ----
    pred_l, ff_l, fr_l, dn_l, slide_l, sid_l, x_l, y_l = ([] for _ in range(8))
    ct_cols = None
    spots_by_slide = {}
    for s in SLIDES:
        df = pd.read_csv(INFER / s / "predictions.csv")
        ct_cols = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
        pred = np.load(INFER / s / "predictions.npy")
        ff = np.load(INFER / s / "features_fused.npy")
        fr = np.load(INFER / s / "features_resnet.npy")
        dn = np.load(DINO / s / "features_dinov2.npy")
        sp = pd.read_csv(GRAPH / f"{s}_spots.csv")
        spots_by_slide[s] = sp
        pred_l.append(pred); ff_l.append(ff); fr_l.append(fr); dn_l.append(dn)
        slide_l += [s] * len(pred)
        sid_l += list(df["spot_id"])
        x_l += list(df["X"]); y_l += list(df["Y"])
        print(f"  loaded {short(s)}: {len(pred)} spots")
    pred = np.vstack(pred_l); ff = np.vstack(ff_l); fr = np.vstack(fr_l); dn = np.vstack(dn_l)
    slide_arr = np.array(slide_l); sid_arr = np.array(sid_l)
    x_arr = np.array(x_l); y_arr = np.array(y_l)
    dom_idx = pred.argmax(1)
    dom_ct = np.array([ct_cols[i] for i in dom_idx])
    N = len(pred)
    print(f"[{a.res_label}] total {N} spots, distinct dominant types: {len(set(dom_ct))}")

    REPS = [
        ("prediction_log1p", np.log1p(pred)),
        ("features_fused", ff),
        ("features_resnet", fr),
        ("features_dinov2", dn),
    ]

    # ---- cluster (=dominant ct) 선정: spot 수 desc, >= MIN_GRID ----
    from collections import Counter
    cnt = Counter(dom_ct.tolist())
    grid_cts = [ct for ct, c in cnt.most_common() if c >= MIN_GRID]
    rare_cts = [(ct, c) for ct, c in cnt.most_common() if c < MIN_GRID]
    print(f"  grid clusters ({len(grid_cts)}): " +
          ", ".join(f"{ct}({cnt[ct]})" for ct in grid_cts))
    print(f"  rare (<{MIN_GRID}, grid 생략): " + ", ".join(f"{ct}({c})" for ct, c in rare_cts))

    # 색: grid cluster 는 tab20, 나머지는 회색
    cmap = plt.get_cmap("tab20", max(20, len(grid_cts)))
    color_map = {ct: cmap(i % cmap.N) for i, ct in enumerate(grid_cts)}
    OTHER = (0.7, 0.7, 0.7, 1.0)

    def point_colors():
        return [color_map.get(ct, OTHER) for ct in dom_ct]

    # ---- (2) 4-panel UMAP, color = dominant ct ----
    fig, axes = plt.subplots(1, len(REPS), figsize=(5.5 * len(REPS), 6.5))
    pcs = point_colors()
    for ax, (name, X) in zip(axes, REPS):
        cache = EMB / f"combined_{name}_umap2d.npy"
        if cache.exists():
            z = np.load(cache)
        else:
            print(f"  [umap] {name} ...", flush=True)
            z = umap.UMAP(**UMAP_KW).fit_transform(X)
            np.save(cache, z)
        ax.scatter(z[:, 0], z[:, 1], c=pcs, s=3, alpha=0.5, linewidths=0)
        ax.set_title(f"{name}  ({X.shape[1]}-d)")
        ax.set_xticks([]); ax.set_yticks([])
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=color_map[ct],
                          markersize=8, label=f"{ct} ({cnt[ct]})") for ct in grid_cts]
    handles.append(plt.Line2D([], [], marker="o", linestyle="", color=OTHER,
                              markersize=8, label=f"other rare (<{MIN_GRID})"))
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5),
               fontsize=8, title="Hist2Cell dominant cell type")
    fig.suptitle(
        f"[{a.res_label}] cross-slide UMAP (4 rep), color = Hist2Cell dominant cell type — "
        f"{N} spots. DINO cluster 가 cell-type 과 모이면 morphology↔cell-type 일치.",
        fontsize=12)
    fig.tight_layout(rect=[0, 0, 0.86, 0.96])
    p = OUT / "umap_4rep_by_dominant_ct.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  saved: {p}")

    # ---- (3)(4) cluster 별 centroid 최근접 패치 grid + CSV ----
    sl_handles = {s: openslide.OpenSlide(find_svs(a.wsi_dir, s)) for s in SLIDES}
    csv_rows = []
    for rank, ct in enumerate(grid_cts, 1):
        mask = np.where(dom_ct == ct)[0]
        cen = dn[mask].mean(0)
        d = np.linalg.norm(dn[mask] - cen, axis=1)
        order = mask[np.argsort(d)]                 # 거리 오름차순 (대표성 순)
        sel = order[:min(MAX_PATCH, len(order))]
        k = len(sel)
        g = int(np.ceil(np.sqrt(k)))               # 정방형 한 변
        canvas = np.full((g * DISP_PX, g * DISP_PX, 3), 255, np.uint8)
        for gi, idx in enumerate(sel):
            s = slide_arr[idx]
            sp = spots_by_slide[s]
            row = sp[sp["spot_id"] == sid_arr[idx]].iloc[0]
            tl_x, tl_y, ts = int(row["tile_x_topleft"]), int(row["tile_y_topleft"]), int(row["tile_size"])
            patch = sl_handles[s].read_region((tl_x, tl_y), 0, (ts, ts)).convert("RGB")
            patch = patch.resize((DISP_PX, DISP_PX), Image.BILINEAR)
            r, cc = gi // g, gi % g
            canvas[r * DISP_PX:(r + 1) * DISP_PX, cc * DISP_PX:(cc + 1) * DISP_PX] = np.asarray(patch)
            csv_rows.append({
                "res": a.res_label, "cluster": ct, "cluster_rank": rank,
                "grid_index": gi, "grid_row": r, "grid_col": cc,
                "slide": short(s), "spot_id": sid_arr[idx],
                "x": int(x_arr[idx]), "y": int(y_arr[idx]),
                "dist_to_centroid": float(np.sort(d)[gi]),
            })
        out_png = OUT / f"cluster_{rank:02d}_{safe(ct)}.png"
        Image.fromarray(canvas).save(out_png)
        print(f"  saved: {out_png}  ({k} patches, {g}x{g})")
    for s in sl_handles.values():
        s.close()

    csv_df = pd.DataFrame(csv_rows)
    csv_path = OUT / f"dino_clusters_{a.res_label}.csv"
    csv_df.to_csv(csv_path, index=False)
    print(f"  saved: {csv_path}  ({len(csv_df)} patch rows, {len(grid_cts)} clusters)")
    print("done.")


if __name__ == "__main__":
    main()
