"""Visium human lung GT 셋 — GT vs Hist2Cell vs DINO 비교 (HEX 도입 전 baseline).

라벨·기준 = 실제 cell2location GT (TCGA 와 달리 Hist2Cell 예측 아님 → circular 아님).

샘플별로:
  reps = GT(log1p 80-d) / Hist2Cell prediction(log1p 80-d) / DINO(768-d), 모두 per-dim z-score.
  label = GT dominant cell type = argmax(GT 80-d).

산출 (--out-dir):
  knn_purity_by_GT.csv        — rep × sample, GT-dominant-type kNN purity(k=10)
  hist2cell_vs_GT.csv         — sample 별 dominant 일치율 + mean per-celltype Pearson r + pooled r
  umap_by_GT_dominant.png     — 샘플(행) × [GT/Hist2Cell/DINO](열) UMAP, color=GT dominant type
  summary.md (별도 작성)

Usage:
  .venv/bin/python lung_pilot/visium_gt_compare.py --out-dir lung_pilot/visium_gt/compare
"""
import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parent.parent
GT = ROOT / "lung_pilot" / "visium_gt" / "gt"
DINO = ROOT / "lung_pilot" / "visium_gt" / "dino_output"
H2C = ROOT / "lung_pilot" / "visium_gt" / "hist2cell_output"
SAMPLES = ["WSA_LngSP10193347", "WSA_LngSP8759313", "WSA_LngSP9258468"]
COLS = ["GT", "Hist2Cell", "DINO"]
UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)


def zscore(X):
    sd = X.std(0); sd[sd == 0] = 1.0
    return (X - X.mean(0)) / sd


def knn_purity(X, lab, k=10):
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    lab = np.asarray(lab)
    return float((lab[idx[:, 1:]] == lab[:, None]).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(ROOT / "lung_pilot" / "visium_gt" / "compare"))
    a = ap.parse_args()
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    EMB = OUT / "embeddings"; EMB.mkdir(exist_ok=True)
    ct_names = np.load(GT / "cell_types.npy", allow_pickle=True)

    data = {}
    all_dom = []
    for s in SAMPLES:
        gt = np.load(GT / f"{s}_celltype_gt.npy")          # [N,80] GT abundance
        pred = np.load(H2C / s / "predictions.npy")        # [N,80] Hist2Cell
        dino = np.load(DINO / s / "features_dinov2.npy")   # [N,768]
        dom = ct_names[gt.argmax(1)]                       # GT dominant type
        data[s] = dict(gt=gt, pred=pred, dino=dino, dom=dom)
        all_dom += dom.tolist()
        print(f"{s}: N={len(gt)} gt{gt.shape} pred{pred.shape} dino{dino.shape}")

    cnt = Counter(all_dom)
    major = [c for c, n in cnt.most_common() if n >= 40]
    cmap = plt.get_cmap("tab20", max(20, len(major)))
    cmap_d = {c: cmap(i % cmap.N) for i, c in enumerate(major)}
    OTHER = (0.7, 0.7, 0.7, 1.0)

    def reps(s):
        d = data[s]
        return {"GT": zscore(np.log1p(d["gt"])),
                "Hist2Cell": zscore(np.log1p(d["pred"])),
                "DINO": zscore(d["dino"])}

    # --- kNN purity by GT dominant ---
    prows = []
    for s in SAMPLES:
        R = reps(s)
        for c in COLS:
            prows.append({"sample": s, "rep": c,
                          "knn_purity_GTdom_k10": round(knn_purity(R[c], data[s]["dom"]), 4),
                          "note": "label=argmax(GT); GT rep semi-circular" if c == "GT" else ""})
    pdf = pd.DataFrame(prows); pdf.to_csv(OUT / "knn_purity_by_GT.csv", index=False)
    print("\n== kNN purity by GT dominant type (k=10) =="); print(pdf.to_string(index=False))

    # --- Hist2Cell vs GT accuracy ---
    arows = []
    for s in SAMPLES:
        gt, pred = data[s]["gt"], data[s]["pred"]
        dom_acc = float((gt.argmax(1) == pred.argmax(1)).mean())
        # per-celltype Pearson (GT std>0)
        rs = []
        for j in range(gt.shape[1]):
            if gt[:, j].std() > 1e-8 and pred[:, j].std() > 1e-8:
                rs.append(np.corrcoef(gt[:, j], pred[:, j])[0, 1])
        pooled = np.corrcoef(gt.flatten(), pred.flatten())[0, 1]
        arows.append({"sample": s, "dominant_match_acc": round(dom_acc, 4),
                      "mean_per_celltype_pearson": round(float(np.nanmean(rs)), 4),
                      "n_celltypes_scored": len(rs),
                      "pooled_pearson": round(float(pooled), 4)})
    adf = pd.DataFrame(arows); adf.to_csv(OUT / "hist2cell_vs_GT.csv", index=False)
    print("\n== Hist2Cell vs GT =="); print(adf.to_string(index=False))

    # --- UMAP: 샘플(행) × rep(열), color=GT dominant ---
    fig, axes = plt.subplots(len(SAMPLES), len(COLS), figsize=(6 * len(COLS), 6 * len(SAMPLES)))
    for r, s in enumerate(SAMPLES):
        R = reps(s); dom = data[s]["dom"]
        pc = [cmap_d.get(c, OTHER) for c in dom]
        for c, col in enumerate(COLS):
            ax = axes[r, c]
            cache = EMB / f"{s}_{col}_umap2d.npy"
            z = np.load(cache) if cache.exists() else umap.UMAP(**UMAP_KW).fit_transform(R[col])
            if not cache.exists():
                np.save(cache, z)
            ax.scatter(z[:, 0], z[:, 1], c=pc, s=5, alpha=0.6, linewidths=0)
            pr = pdf[(pdf["sample"] == s) & (pdf.rep == col)]["knn_purity_GTdom_k10"].iloc[0]
            ax.set_title((col if r == 0 else "") + f"\npurity={pr}", fontsize=12)
            ax.set_xticks([]); ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(f"{s.replace('WSA_LngSP','')}\n({len(dom)})", fontsize=11)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=cmap_d[c], markersize=8, label=c)
               for c in major]
    handles.append(plt.Line2D([], [], marker="o", linestyle="", color=OTHER, markersize=8, label="other"))
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5), fontsize=9,
               title="GT dominant cell type")
    fig.suptitle("Visium human lung GT — GT vs Hist2Cell vs DINO  "
                 "(color = 실제 GT dominant cell type, per-dim z-scored, title=kNN purity k10)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 0.9, 0.98])
    p = OUT / "umap_by_GT_dominant.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"\nsaved: {p}\ndone.")


if __name__ == "__main__":
    main()
