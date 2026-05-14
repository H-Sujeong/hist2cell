"""Per-cell-type top pair extraction.

사용자 요청: "per-cell-type 별로 강한 페어 어떤 게 있나 보여줘"

joint H (94, 80) × P (94, 4168) 의 모든 페어 Pearson r 계산 후, 각
cell type 별로 양/음 top-K 페어 정리.  glob BH-FDR 보정 후 통과 표시.

산출:
 per_celltype_top_pairs.csv     — 80 cell type x top 5 (pos + neg) = 800 행
 per_celltype_max_r.csv          — 80 cell type 의 max +r / max -r 요약
 per_celltype_top_overview.png   — 3 panel:
   A) 80 cell type 의 max +r 막대 (cell type 별 best positive 페어)
   B) 80 cell type 의 max |-r| 막대 (best negative)
   C) 전체 joint top 30 페어 (label = "cell_type :: gene")
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from _proof_ver2_lib import build_roi_signatures, load_proteomics_matrix, align_modalities
from build_joint_cca import slide_cfg

TOP_K = 5  # per cell type, both directions


def build_joint():
    out = {}
    for sid in ["slide1", "slide2"]:
        cfg = slide_cfg(sid)
        sig_df, cell_cols = build_roi_signatures(cfg)
        log2_f, slide_cols = load_proteomics_matrix(cfg)
        common, H, P, _, genes = align_modalities(sig_df, log2_f, slide_cols, cell_cols)
        out[sid] = {"H": H, "P": P, "genes": genes, "cell_cols": cell_cols}
    common_genes = sorted(set(out["slide1"]["genes"]) & set(out["slide2"]["genes"]))
    for sid in ["slide1", "slide2"]:
        idx = [out[sid]["genes"].index(g) for g in common_genes]
        out[sid]["P_common"] = out[sid]["P"][:, idx]
    H = np.vstack([StandardScaler().fit_transform(out["slide1"]["H"]),
                    StandardScaler().fit_transform(out["slide2"]["H"])])
    P = np.vstack([StandardScaler().fit_transform(out["slide1"]["P_common"]),
                    StandardScaler().fit_transform(out["slide2"]["P_common"])])
    return H, P, out["slide1"]["cell_cols"], common_genes


def all_pair_r(H, P):
    n = H.shape[0]
    Hc = H - H.mean(axis=0, keepdims=True)
    Pc = P - P.mean(axis=0, keepdims=True)
    h_std = H.std(axis=0, ddof=1); p_std = P.std(axis=0, ddof=1)
    h_std = np.where(h_std == 0, 1, h_std); p_std = np.where(p_std == 0, 1, p_std)
    r = (Hc.T @ Pc) / ((n - 1) * h_std[:, None] * p_std[None, :])
    t_stat = r * np.sqrt((n - 2) / np.clip(1 - r**2, 1e-12, None))
    p = 2 * (1 - student_t.cdf(np.abs(t_stat), df=n - 2))
    return r, p


def main():
    print("[*] computing joint H,P + all-pair r …")
    H, P, cell_cols, genes = build_joint()
    print(f"    H={H.shape}, P={P.shape}, total pairs={H.shape[1]*P.shape[1]:,}")
    r, p = all_pair_r(H, P)

    # global BH-FDR
    p_bh = multipletests(p.flatten(), method="fdr_bh")[1].reshape(r.shape)

    rows_top = []
    rows_summary = []
    for i, ct in enumerate(cell_cols):
        r_row = r[i]; p_row = p[i]; pbh_row = p_bh[i]
        # top K positive
        order_pos = np.argsort(-r_row)[:TOP_K]
        # top K negative
        order_neg = np.argsort(r_row)[:TOP_K]
        for rank, j in enumerate(order_pos, start=1):
            rows_top.append({
                "cell_type": ct, "direction": "pos", "rank": rank,
                "gene": genes[j], "r": float(r_row[j]),
                "p": float(p_row[j]), "p_bh_global": float(pbh_row[j]),
                "bh_pass": bool(pbh_row[j] < 0.05),
            })
        for rank, j in enumerate(order_neg, start=1):
            rows_top.append({
                "cell_type": ct, "direction": "neg", "rank": rank,
                "gene": genes[j], "r": float(r_row[j]),
                "p": float(p_row[j]), "p_bh_global": float(pbh_row[j]),
                "bh_pass": bool(pbh_row[j] < 0.05),
            })
        # summary per cell type
        max_pos_j = order_pos[0]; max_neg_j = order_neg[0]
        n_bh_in_row = int((pbh_row < 0.05).sum())
        rows_summary.append({
            "cell_type": ct,
            "max_r_pos": float(r_row[max_pos_j]),
            "max_r_pos_gene": genes[max_pos_j],
            "max_r_neg": float(r_row[max_neg_j]),
            "max_r_neg_gene": genes[max_neg_j],
            "n_genes_BH05": n_bh_in_row,
            "pct_genes_BH05": 100 * n_bh_in_row / len(genes),
        })
    top_df = pd.DataFrame(rows_top)
    sum_df = pd.DataFrame(rows_summary)
    top_df.to_csv(HERE / "per_celltype_top_pairs.csv", index=False)
    sum_df.to_csv(HERE / "per_celltype_max_r.csv", index=False)

    print("\n[*] top 15 strongest cell-type pairs (positive direction):")
    print(top_df[top_df.direction == "pos"]
            .nlargest(15, "r")[["cell_type","gene","r","p_bh_global","bh_pass"]]
            .to_string(index=False))
    print("\n[*] top 15 strongest cell-type pairs (negative direction):")
    print(top_df[top_df.direction == "neg"]
            .nsmallest(15, "r")[["cell_type","gene","r","p_bh_global","bh_pass"]]
            .to_string(index=False))

    # ── plotting ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.30,
                           left=0.10, right=0.96, top=0.95, bottom=0.05,
                           width_ratios=[1.2, 1.0])
    axA = fig.add_subplot(gs[0, 0])  # cell type max +r 막대
    axB = fig.add_subplot(gs[1, 0])  # cell type max -r 막대
    axC = fig.add_subplot(gs[:, 1])  # 전체 top 30 페어

    # A: max +r per cell type, sorted by max +r desc
    sumA = sum_df.sort_values("max_r_pos", ascending=True)  # bottom-up bar
    y = np.arange(len(sumA))
    colors = ["#d62728" if r > 0.5 else ("#ff7f0e" if r > 0.3 else "#999999")
              for r in sumA["max_r_pos"]]
    axA.barh(y, sumA["max_r_pos"], color=colors, edgecolor="black", linewidth=0.4)
    for i, (val, gene) in enumerate(zip(sumA["max_r_pos"], sumA["max_r_pos_gene"])):
        axA.text(val + 0.005, i, gene, fontsize=6.5, va="center")
    axA.set_yticks(y)
    axA.set_yticklabels(sumA["cell_type"], fontsize=6.5)
    axA.axvline(0.3, c="#ff7f0e", lw=0.7, linestyle="--", alpha=0.6)
    axA.axvline(0.5, c="#d62728", lw=0.7, linestyle="--", alpha=0.6)
    axA.set_xlabel("max positive r  (best gene pair for the cell type)")
    axA.set_title(
        "A) Per-cell-type max +r  (label = top hit gene)\n"
        "red bar = strong (>0.5), orange = moderate (0.3-0.5), grey = weak",
        fontsize=10)
    axA.set_xlim(0, max(0.8, sumA["max_r_pos"].max() + 0.15))
    axA.grid(axis="x", alpha=0.2)

    # B: max -r per cell type
    sumB = sum_df.sort_values("max_r_neg", ascending=False)  # most negative bottom
    y = np.arange(len(sumB))
    colors = ["#1f77b4" if r < -0.5 else ("#9ecae1" if r < -0.3 else "#cccccc")
              for r in sumB["max_r_neg"]]
    axB.barh(y, sumB["max_r_neg"], color=colors, edgecolor="black", linewidth=0.4)
    for i, (val, gene) in enumerate(zip(sumB["max_r_neg"], sumB["max_r_neg_gene"])):
        axB.text(val - 0.005, i, gene, fontsize=6.5, va="center", ha="right")
    axB.set_yticks(y)
    axB.set_yticklabels(sumB["cell_type"], fontsize=6.5)
    axB.axvline(-0.3, c="#9ecae1", lw=0.7, linestyle="--", alpha=0.6)
    axB.axvline(-0.5, c="#1f77b4", lw=0.7, linestyle="--", alpha=0.6)
    axB.set_xlabel("max negative r  (best inverse-correlated gene)")
    axB.set_title("B) Per-cell-type max -r  (label = top inverse gene)\n"
                   "blue bar = strong inverse (<-0.5), light blue = moderate inverse",
                   fontsize=10)
    axB.set_xlim(min(-0.8, sumB["max_r_neg"].min() - 0.15), 0)
    axB.grid(axis="x", alpha=0.2)

    # C: 전체 top 30 페어 (양/음 합쳐, |r| 기준)
    top30 = (top_df.assign(abs_r=top_df["r"].abs())
                    .nlargest(30, "abs_r")
                    .iloc[::-1])  # smallest abs_r at top → bottom bigger
    y = np.arange(len(top30))
    colors = ["#d62728" if r > 0 else "#1f77b4" for r in top30["r"]]
    axC.barh(y, top30["r"], color=colors, edgecolor="black", linewidth=0.4)
    labels = [f"{ct} :: {g}" for ct, g in zip(top30["cell_type"], top30["gene"])]
    axC.set_yticks(y)
    axC.set_yticklabels(labels, fontsize=7.5, family="monospace")
    for i, (r_val, bh) in enumerate(zip(top30["r"], top30["p_bh_global"])):
        star = " ★" if bh < 0.001 else (" ✦" if bh < 0.01 else (" +" if bh < 0.05 else "")).replace("✦", "*").replace("★", "**")
        x_offset = 0.01 if r_val > 0 else -0.01
        ha = "left" if r_val > 0 else "right"
        axC.text(r_val + x_offset, i, f"{r_val:+.3f}{star}",
                  fontsize=7, va="center", ha=ha)
    axC.axvline(0, c="k", lw=0.5)
    axC.set_xlabel("Pearson r")
    axC.set_title("C) Top 30 strongest joint pairs by |r|\n"
                   "red = positive, blue = negative.  ** BH<0.001, * BH<0.01, + BH<0.05",
                   fontsize=10)
    axC.grid(axis="x", alpha=0.2)
    axC.set_xlim(top30["r"].min() - 0.12, top30["r"].max() + 0.12)

    fig.suptitle(
        "Per-cell-type top pairs  (joint H x P, n=94 ROI, 4168 common genes)",
        fontsize=13, fontweight="bold", y=0.985)
    fig.savefig(HERE / "per_celltype_top_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[✓] per_celltype_top_overview.png")
    print(f"    output: per_celltype_top_pairs.csv (800 rows) + per_celltype_max_r.csv (80 rows)")


if __name__ == "__main__":
    main()
