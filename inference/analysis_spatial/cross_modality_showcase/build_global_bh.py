"""Global BH-FDR test — all (cell type × gene) pairs.

가장 직접적인 답: "Hist2Cell phenotyping 과 ROI proteomics 가
통계적으로 유의한가?"  80 × G_common 페어 전체 위의 Pearson r 분포
+ 글로벌 BH-FDR 보정 + permutation 영가설.

산출:
 global_pair_correlations.csv     (every pair: r, p, p_bh)
 global_pair_summary.csv          (통과 페어 수 by significance level)
 global_pair_distribution.png     (r histogram, BH 통과 분포)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, t as student_t
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, RANDOM_SEED,
)
from build_joint_cca import slide_cfg

N_PERM = 200  # null distribution; 200 × 333k = 67M r calcs, plenty


def build_joint():
    out = {}
    for sid in ["slide1", "slide2"]:
        cfg = slide_cfg(sid)
        sig_df, cell_cols = build_roi_signatures(cfg)
        log2_f, slide_cols = load_proteomics_matrix(cfg)
        common, H, P, sig_aligned, genes = align_modalities(
            sig_df, log2_f, slide_cols, cell_cols)
        out[sid] = {"H": H, "P": P, "genes": genes, "cell_cols": cell_cols}
    g1 = set(out["slide1"]["genes"]); g2 = set(out["slide2"]["genes"])
    common_genes = sorted(g1 & g2)
    for sid in ["slide1", "slide2"]:
        idx = [out[sid]["genes"].index(g) for g in common_genes]
        out[sid]["P_common"] = out[sid]["P"][:, idx]
    H1z = StandardScaler().fit_transform(out["slide1"]["H"])
    H2z = StandardScaler().fit_transform(out["slide2"]["H"])
    P1z = StandardScaler().fit_transform(out["slide1"]["P_common"])
    P2z = StandardScaler().fit_transform(out["slide2"]["P_common"])
    H = np.vstack([H1z, H2z])
    P = np.vstack([P1z, P2z])
    return H, P, out["slide1"]["cell_cols"], common_genes


def all_pair_r(H, P):
    """Return (r_matrix [n_cells, n_genes], p_matrix)."""
    n = H.shape[0]
    Hc = H - H.mean(axis=0, keepdims=True)
    Pc = P - P.mean(axis=0, keepdims=True)
    h_std = H.std(axis=0, ddof=1)
    p_std = P.std(axis=0, ddof=1)
    # avoid div by zero
    h_std = np.where(h_std == 0, 1, h_std)
    p_std = np.where(p_std == 0, 1, p_std)
    r = (Hc.T @ Pc) / ((n - 1) * h_std[:, None] * p_std[None, :])
    # t-statistic and two-sided p
    t_stat = r * np.sqrt((n - 2) / np.clip(1 - r**2, 1e-12, None))
    p = 2 * (1 - student_t.cdf(np.abs(t_stat), df=n - 2))
    return r, p


def main():
    print("[*] building joint matrices …")
    H, P, cell_cols, genes = build_joint()
    n_cells, n_genes = H.shape[1], P.shape[1]
    print(f"    joint H = {H.shape}, joint P = {P.shape}")
    print(f"    total pairs = {n_cells} x {n_genes} = {n_cells * n_genes:,}")

    print("[*] computing all-pair Pearson r + parametric p …")
    r, p = all_pair_r(H, P)
    n_pairs = n_cells * n_genes

    # 글로벌 BH-FDR
    print("[*] global BH-FDR …")
    p_flat = p.flatten()
    r_flat = r.flatten()
    _, p_bh, _, _ = multipletests(p_flat, method="fdr_bh")
    pass_5  = int((p_bh < 0.05).sum())
    pass_1  = int((p_bh < 0.01).sum())
    pass_01 = int((p_bh < 0.001).sum())
    print(f"    BH<0.05 = {pass_5:,}/{n_pairs:,} ({100*pass_5/n_pairs:.1f}%)")
    print(f"    BH<0.01 = {pass_1:,}/{n_pairs:,} ({100*pass_1/n_pairs:.1f}%)")
    print(f"    BH<.001 = {pass_01:,}/{n_pairs:,} ({100*pass_01/n_pairs:.1f}%)")

    # permutation null on number of BH<0.05 pairs (ROI shuffle on proteomics)
    print(f"[*] permutation null on BH-passing count ({N_PERM} iters) …")
    rng = np.random.default_rng(RANDOM_SEED)
    null_pass_5 = np.zeros(N_PERM, dtype=int)
    for k in range(N_PERM):
        perm = rng.permutation(H.shape[0])
        r_null, p_null = all_pair_r(H, P[perm])
        _, p_bh_null, _, _ = multipletests(p_null.flatten(), method="fdr_bh")
        null_pass_5[k] = int((p_bh_null < 0.05).sum())
        if (k + 1) % 50 == 0:
            print(f"   perm {k+1}/{N_PERM}, BH<0.05 count median so far = "
                  f"{int(np.median(null_pass_5[:k+1]))}")
    null_mean = float(null_pass_5.mean())
    null_95_hi = float(np.percentile(null_pass_5, 97.5))
    p_count = float((null_pass_5 >= pass_5).mean())
    print(f"    null mean BH<0.05 = {null_mean:.0f}, 95% upper = {null_95_hi:.0f}")
    print(f"    permutation p = {p_count}")

    # save table
    pair_df = pd.DataFrame({
        "cell_type": np.repeat(cell_cols, n_genes),
        "gene":      np.tile(genes, n_cells),
        "r":         r_flat,
        "p":         p_flat,
        "p_bh":      p_bh,
    })
    pair_df_top = pair_df.sort_values("p_bh").head(50000)
    pair_df_top.to_csv(HERE / "global_pair_correlations.csv", index=False)

    summary = pd.DataFrame([{
        "total_pairs": n_pairs,
        "BH_0.05_count": pass_5,
        "BH_0.05_pct": 100 * pass_5 / n_pairs,
        "BH_0.01_count": pass_1,
        "BH_0.001_count": pass_01,
        "null_BH_0.05_mean": null_mean,
        "null_BH_0.05_95_upper": null_95_hi,
        "permutation_p_on_count": p_count,
    }])
    summary.to_csv(HERE / "global_pair_summary.csv", index=False)

    # ── plot ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    # A: histogram of r
    axes[0].hist(r_flat, bins=120, color="#444444", alpha=0.7, edgecolor="white")
    axes[0].axvline(0, c="k", lw=0.6, alpha=0.4)
    axes[0].set_xlabel("Pearson r  (per cell-type x gene pair)")
    axes[0].set_ylabel("number of pairs")
    axes[0].set_title(
        f"A) Distribution of all {n_pairs:,} pair-wise Pearson r\n"
        f"mean r = {r_flat.mean():+.3f}, |r|>0.3: {int((np.abs(r_flat)>0.3).sum()):,} pairs",
        fontsize=10)

    # B: BH-passing counts (observed vs null)
    axes[1].hist(null_pass_5, bins=30, color="#bbbbbb", alpha=0.8,
                  edgecolor="white", label=f"null ({N_PERM} ROI-shuffles)")
    axes[1].axvline(pass_5, c="#d62728", lw=2,
                     label=f"observed BH<0.05 = {pass_5:,}")
    axes[1].axvline(null_95_hi, c="black", lw=0.7, linestyle="--",
                     label=f"null 97.5% = {null_95_hi:.0f}")
    axes[1].set_xlabel("# pairs passing global BH-FDR<0.05")
    axes[1].set_ylabel("permutation count")
    axes[1].set_title(
        f"B) Global BH-FDR<0.05 passing count vs null\n"
        f"observed {pass_5:,} pairs ({100*pass_5/n_pairs:.1f}% of total),  "
        f"permutation p = {p_count}",
        fontsize=10)
    axes[1].legend(loc="upper right", fontsize=8, frameon=False)

    # C: bar of significance thresholds
    levels = ["BH<0.05", "BH<0.01", "BH<0.001"]
    counts = [pass_5, pass_1, pass_01]
    pcts = [100*c/n_pairs for c in counts]
    bars = axes[2].bar(levels, counts, color=["#d62728","#ff7f0e","#bcbd22"],
                        edgecolor="black")
    for bar, cnt, pct in zip(bars, counts, pcts):
        axes[2].text(bar.get_x() + bar.get_width()/2,
                      bar.get_height(),
                      f"{cnt:,}\n({pct:.1f}%)",
                      ha="center", va="bottom", fontsize=9)
    axes[2].set_ylabel("# passing pairs")
    axes[2].set_title("C) Passing pair counts by BH-FDR threshold",
                       fontsize=10)
    axes[2].set_ylim(0, max(counts) * 1.18)
    axes[2].grid(axis="y", alpha=0.2)

    fig.suptitle(
        "Global pairwise BH-FDR test  -  Hist2Cell phenotyping (80 cell types) "
        f"x  ROI proteomics ({n_genes:,} genes)  on  n=94 joint ROIs",
        fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(HERE / "global_pair_distribution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[✓] global_pair_distribution.png")
    print(f"\n=== ANSWER ===")
    print(f"Of {n_pairs:,} (cell type x gene) pairs, "
          f"{pass_5:,} pass global BH-FDR<0.05.")
    print(f"Null permutation: under random shuffling, expected ~{null_mean:.0f} "
          f"(95% upper {null_95_hi:.0f}).")
    print(f"Observed >> null  ->  permutation p = {p_count}")


if __name__ == "__main__":
    main()
