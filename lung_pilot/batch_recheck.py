"""Batch effect re-check — 정확한 chance baseline + balanced subsample UMAP.

배경: compare_hist2cell_vs_dinov2.py 의 1-NN purity 보고서가 chance
baseline 을 단순 `1/3 = 0.33` 으로만 표시해 prediction 0.775 가 "강한
batch" 처럼 보였다. 실제 chance 는 slide-size 가중 ∑ pᵢ² 이고, 4390-BS1
이 spot 의 69% 라 chance ≈ 0.53. 또한 cross-slide UMAP 의 시각적 4390
dominance 는 단순 spot 수 5배 효과가 크다.

본 스크립트가 하는 일:
  1) size-weighted chance + per-rep excess (= (purity-chance)/(1-chance))
     계산해 compare_output/metrics_corrected.csv 저장.
  2) 각 슬라이드 1,871 spot (가장 작은 TS1 기준) 으로 random sample 한
     balanced cross-slide UMAP 을 4 representation 별로 fit 해
     umap_output/cross_slide_balanced.png 저장. 시각적으로 spot-수
     효과를 제거한 진짜 batch 모양.

Usage:
    .venv/bin/python lung_pilot/batch_recheck.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import umap

ROOT = Path(__file__).resolve().parent.parent
INFER_DIR = ROOT / "lung_pilot" / "inference_output"
DINO_DIR = ROOT / "lung_pilot" / "dino_output"
UMAP_OUT = ROOT / "lung_pilot" / "umap_output"
CMP_OUT = ROOT / "lung_pilot" / "compare_output"

SLIDES = [
    "TCGA-05-4245-01A-01-BS1",
    "TCGA-05-4245-01A-01-TS1",
    "TCGA-05-4390-01A-01-BS1",
]

REPS = [
    ("prediction_log1p", lambda d: np.log1p(d["pred"]), 80),
    ("features_fused",   lambda d: d["ff"],            256),
    ("features_resnet",  lambda d: d["fr"],            512),
    ("features_dinov2",  lambda d: d["dn"],            768),
]

UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)
SEED = 42


def load_slide(s):
    df = pd.read_csv(INFER_DIR / s / "predictions.csv")
    return dict(
        df=df,
        pred=np.load(INFER_DIR / s / "predictions.npy"),
        ff=np.load(INFER_DIR / s / "features_fused.npy"),
        fr=np.load(INFER_DIR / s / "features_resnet.npy"),
        dn=np.load(DINO_DIR / s / "features_dinov2.npy"),
    )


def slide_1nn_purity(X, slide_lbl):
    nn = NearestNeighbors(n_neighbors=2, algorithm="auto", n_jobs=-1)
    nn.fit(X)
    _, idx = nn.kneighbors(X)
    return float((slide_lbl[idx[:, 1]] == slide_lbl).mean())


def short_slide(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def main():
    rng = np.random.default_rng(SEED)
    payloads = {s: load_slide(s) for s in SLIDES}
    sizes = {s: len(payloads[s]["df"]) for s in SLIDES}
    n_total = sum(sizes.values())
    p_share = {s: sizes[s] / n_total for s in SLIDES}
    chance_weighted = float(sum(p ** 2 for p in p_share.values()))
    chance_uniform = 1.0 / len(SLIDES)
    print(f"slide sizes: {sizes}  total={n_total}")
    print(f"  p_share: {p_share}")
    print(f"  size-weighted chance (sum p_i^2): {chance_weighted:.4f}")
    print(f"  uniform chance (1/k):            {chance_uniform:.4f}")

    # ---------- (1) reload existing purity from metrics.csv, add excess ----------
    metrics_csv = CMP_OUT / "metrics.csv"
    metrics = pd.read_csv(metrics_csv)
    purity_rows = metrics[metrics["metric"] == "slide_1nn_purity"][["representation", "value"]]
    print("\nslide 1-NN purity (from metrics.csv):")
    print(purity_rows.to_string(index=False))

    purity_rows = purity_rows.copy()
    purity_rows["chance_uniform"] = chance_uniform
    purity_rows["chance_weighted"] = chance_weighted
    purity_rows["excess_over_weighted"] = (
        (purity_rows["value"] - chance_weighted) / (1 - chance_weighted)
    )
    purity_rows.rename(columns={"value": "purity"}, inplace=True)
    out_csv = CMP_OUT / "metrics_corrected.csv"
    purity_rows.to_csv(out_csv, index=False)
    print(f"\nsaved: {out_csv}")
    print(purity_rows.to_string(index=False))

    # ---------- (2) balanced subsample cross-slide UMAP ----------
    n_per = min(sizes.values())
    print(f"\nbalanced subsample: {n_per} spots per slide  (total {n_per * len(SLIDES)})")
    sub_idx = {
        s: np.sort(rng.choice(sizes[s], size=n_per, replace=False)) for s in SLIDES
    }

    fig, axes = plt.subplots(1, len(REPS), figsize=(5.5 * len(REPS), 6.5))
    slide_colors = dict(zip(SLIDES, plt.get_cmap("tab10").colors[: len(SLIDES)]))
    for ax, (rep_name, fn, dim) in zip(axes, REPS):
        print(f"  fitting balanced UMAP for {rep_name} ({dim}-d)...")
        Xs, lbl = [], []
        for s in SLIDES:
            X_full = fn(payloads[s])
            Xs.append(X_full[sub_idx[s]])
            lbl.extend([s] * n_per)
        X_all = np.vstack(Xs).astype(np.float32)
        lbl = np.array(lbl)
        z = umap.UMAP(**UMAP_KW).fit_transform(X_all)
        # purity on balanced
        p_bal = slide_1nn_purity(X_all, lbl)
        for s in SLIDES:
            m = lbl == s
            ax.scatter(z[m, 0], z[m, 1], c=[slide_colors[s]],
                       s=5, alpha=0.55, linewidths=0, label=short_slide(s))
        ax.set_title(f"{rep_name}  ({dim}-d)\nbalanced 1-NN purity = {p_bal:.3f}",
                     fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")

    handles = [plt.Line2D([], [], marker="o", linestyle="",
                          color=slide_colors[s], markersize=8, label=short_slide(s))
               for s in SLIDES]
    fig.legend(handles=handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=9, title="slide")
    fig.suptitle(
        f"Balanced subsample cross-slide UMAP — {n_per} spots per slide ({n_per*len(SLIDES)} total)   "
        f"(equal-size sampling removes spot-count dominance; chance 1-NN purity = 1/3 = {chance_uniform:.3f})",
        fontsize=11.5,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 0.95])
    out_png = UMAP_OUT / "cross_slide_balanced.png"
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()
