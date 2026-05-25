"""Slide-별 TOP10 cell type — UMAP abundance overlay + 통계.

Breast 분석 (inference/analysis_spatial/1_085_12/cell_typing/analyze_cell_typing.py
`plot_spatial_top10`) 의 *Slide-wide mean abundance 상위 10 cell type* 정의를
lung_pilot 3 슬라이드에 적용.

산출 (lung_pilot/top10_output/):
  top10_stats.csv            — 슬라이드 × 80 cell type 의 mean/std abundance
                                + 슬라이드별 rank
  top10_union.csv            — 3 슬라이드 TOP10 의 union/intersection 표
  top10_<slide>.png          — 슬라이드별 TOP10 cell type abundance overlay
                                (2×5, prediction_log1p UMAP 좌표 위 viridis)
  summary.md                  — 통계 표 + PNG 캡션

UMAP 좌표: lung_pilot/umap_output/embeddings/<slide>_prediction_log1p_umap2d.npy
(umap_subtype_grid.py 가 만들어둔 cache 재사용 — 기존 PNG 와 동일 좌표 보장).

Usage:
    .venv/bin/python lung_pilot/top10_umap.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
INFER_DIR = ROOT / "lung_pilot" / "inference_output"
EMB_DIR = ROOT / "lung_pilot" / "umap_output" / "embeddings"
OUT_DIR = ROOT / "lung_pilot" / "top10_output"

SLIDES = [
    "TCGA-05-4245-01A-01-BS1",
    "TCGA-05-4245-01A-01-TS1",
    "TCGA-05-4390-01A-01-BS1",
]

TOP_K = 10


def short_slide(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def load_slide(s):
    df = pd.read_csv(INFER_DIR / s / "predictions.csv")
    ct_cols = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
    pred = np.load(INFER_DIR / s / "predictions.npy")
    emb = np.load(EMB_DIR / f"{s}_prediction_log1p_umap2d.npy")
    return dict(df=df, ct_cols=ct_cols, pred=pred, emb=emb)


def plot_top10_overlay(slide, payload, top10_types, out_path):
    pred = payload["pred"]
    ct_cols = payload["ct_cols"]
    emb = payload["emb"]
    type_to_idx = {c: i for i, c in enumerate(ct_cols)}

    n_row, n_col = 2, 5
    fig, axes = plt.subplots(n_row, n_col, figsize=(4.0 * n_col, 4.0 * n_row))
    for rank, ct in enumerate(top10_types):
        ax = axes[rank // n_col, rank % n_col]
        vals = pred[:, type_to_idx[ct]]
        order = np.argsort(vals)  # low → high; high on top
        sc = ax.scatter(emb[order, 0], emb[order, 1],
                        c=vals[order], cmap="viridis", s=4, alpha=0.85,
                        linewidths=0, vmin=0, vmax=float(vals.max()))
        ax.set_title(
            f"#{rank + 1}  {ct}\nmean={vals.mean():.3f}  max={vals.max():.2f}  "
            f"frac_pos(>0.1)={float((vals > 0.1).mean()):.2f}",
            fontsize=9,
        )
        ax.set_xticks([]); ax.set_yticks([])
        cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle(
        f"{short_slide(slide)}  —  TOP{TOP_K} cell types by slide-wide mean abundance"
        f"\nViridis overlay on prediction_log1p UMAP (same coords as per_slide PNGs). {len(payload['df'])} spots.",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. compute per-slide mean abundance + TOP10
    payloads = {s: load_slide(s) for s in SLIDES}
    ct_cols = payloads[SLIDES[0]]["ct_cols"]

    stats_rows = []
    top10_per_slide = {}
    for s in SLIDES:
        pred = payloads[s]["pred"]
        mean = pred.mean(axis=0)
        std = pred.std(axis=0)
        frac_pos = (pred > 0.1).mean(axis=0)
        order = np.argsort(-mean)  # desc
        rank = np.empty_like(order); rank[order] = np.arange(1, len(order) + 1)
        for ct_i, ct in enumerate(ct_cols):
            stats_rows.append({
                "slide": short_slide(s),
                "cell_type": ct,
                "mean_abundance": float(mean[ct_i]),
                "std": float(std[ct_i]),
                "frac_pos_over_0.1": float(frac_pos[ct_i]),
                "rank_in_slide": int(rank[ct_i]),
            })
        top10_per_slide[s] = [ct_cols[i] for i in order[:TOP_K]]
        print(f"\n{short_slide(s)} TOP{TOP_K}:")
        for r, ct in enumerate(top10_per_slide[s], 1):
            ct_i = ct_cols.index(ct)
            print(f"  #{r:2d}  {ct:35s}  mean={mean[ct_i]:.4f}  "
                  f"frac>0.1={frac_pos[ct_i]:.3f}")

    stats_df = pd.DataFrame(stats_rows)
    stats_csv = OUT_DIR / "top10_stats.csv"
    stats_df.to_csv(stats_csv, index=False)
    print(f"\nsaved: {stats_csv}  ({len(stats_df)} rows = {len(SLIDES)} × 80)")

    # 2. union / intersection table
    sets = {s: set(top10_per_slide[s]) for s in SLIDES}
    union = sorted(set().union(*sets.values()))
    intersection = sorted(set.intersection(*sets.values()))
    union_rows = []
    for ct in union:
        row = {"cell_type": ct}
        for s in SLIDES:
            in_top = ct in sets[s]
            slide_stats = stats_df[(stats_df["slide"] == short_slide(s)) &
                                   (stats_df["cell_type"] == ct)].iloc[0]
            row[f"{short_slide(s)}_in_top10"] = bool(in_top)
            row[f"{short_slide(s)}_rank"] = int(slide_stats["rank_in_slide"])
            row[f"{short_slide(s)}_mean"] = float(slide_stats["mean_abundance"])
        union_rows.append(row)
    union_df = pd.DataFrame(union_rows)
    union_df["n_slides_in_top10"] = (
        union_df[[f"{short_slide(s)}_in_top10" for s in SLIDES]].sum(axis=1)
    )
    union_df = union_df.sort_values(
        ["n_slides_in_top10",
         f"{short_slide(SLIDES[0])}_mean"],
        ascending=[False, False],
    ).reset_index(drop=True)
    union_csv = OUT_DIR / "top10_union.csv"
    union_df.to_csv(union_csv, index=False)
    print(f"saved: {union_csv}")
    print(f"\nUnion of TOP{TOP_K} across 3 slides: {len(union)} types")
    print(f"Intersection: {len(intersection)} types: {intersection}")

    # 3. per-slide overlay PNG
    for s in SLIDES:
        out_png = OUT_DIR / f"top10_{s}.png"
        plot_top10_overlay(s, payloads[s], top10_per_slide[s], out_png)
        print(f"saved: {out_png}")

    print("done.")


if __name__ == "__main__":
    main()
