"""Effect-size assessment — significant vs strong.

사용자 비판: "유의성은 있는데 연관성이 있다고 하긴 수치들이 낮다."
정확한 지적이므로 *통계적 유의성* 과 *실질적 효과 크기* 를 분리해
보여준다.

산출:
 effect_size_summary.csv — |r| 분포 통계 + bin 별 페어 수
 effect_size_overview.png — 3 panel:
   A) 전체 |r| histogram + percentile 마커 + Cohen 기준선
   B) bin 별 페어 수 (|r| 임계 ladder)
   C) BH<0.05 통과 페어 *자체* 의 r 분포 — "유의한 페어들도 weak 가 다수?"
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
    return r.flatten(), p.flatten()


BINS = [
    (0.00, 0.10, "trivial",    "#dddddd"),
    (0.10, 0.20, "very weak",  "#bbbbbb"),
    (0.20, 0.30, "weak",       "#999999"),
    (0.30, 0.50, "moderate",   "#ff7f0e"),
    (0.50, 0.70, "strong",     "#d62728"),
    (0.70, 1.01, "very strong","#8b0000"),
]


def main():
    print("[*] computing all-pair r …")
    H, P, cell_cols, genes = build_joint()
    r, p = all_pair_r(H, P)
    n_pairs = len(r)
    _, p_bh, _, _ = multipletests(p, method="fdr_bh")

    abs_r = np.abs(r)
    pct = np.percentile(abs_r, [50, 75, 90, 95, 99])
    print(f"    |r| median={pct[0]:.3f}  Q3={pct[1]:.3f}  P90={pct[2]:.3f}  "
          f"P95={pct[3]:.3f}  P99={pct[4]:.3f}  max={abs_r.max():.3f}")
    print(f"    mean r = {r.mean():+.4f}  (signed — small positive bias)")

    # bin counts (전체 + BH 통과)
    is_bh = p_bh < 0.05
    rows = []
    for lo, hi, label, _ in BINS:
        mask = (abs_r >= lo) & (abs_r < hi)
        cnt_all = int(mask.sum())
        cnt_bh = int((mask & is_bh).sum())
        rows.append({"bin": f"|r| {lo:.2f}~{hi:.2f}", "label": label,
                     "n_pairs_all": cnt_all,
                     "pct_all": 100 * cnt_all / n_pairs,
                     "n_pairs_BH05": cnt_bh,
                     "pct_BH05_in_bin": 100 * cnt_bh / max(cnt_all, 1),
                     "pct_of_total_BH05": 100 * cnt_bh / max(is_bh.sum(), 1)})
    binsum = pd.DataFrame(rows)
    print("\n[*] bin distribution:")
    print(binsum.to_string(index=False))

    summary = pd.DataFrame([{
        "n_pairs_total": n_pairs,
        "mean_r": float(r.mean()),
        "median_abs_r": float(pct[0]),
        "P75_abs_r": float(pct[1]),
        "P90_abs_r": float(pct[2]),
        "P95_abs_r": float(pct[3]),
        "P99_abs_r": float(pct[4]),
        "max_abs_r": float(abs_r.max()),
        "n_BH_pass": int(is_bh.sum()),
        "n_BH_pass_abs_r_ge_0.3": int(((abs_r >= 0.3) & is_bh).sum()),
        "n_BH_pass_abs_r_ge_0.5": int(((abs_r >= 0.5) & is_bh).sum()),
        "n_BH_pass_abs_r_ge_0.7": int(((abs_r >= 0.7) & is_bh).sum()),
    }])
    summary.to_csv(HERE / "effect_size_summary.csv", index=False)
    binsum.to_csv(HERE / "effect_size_bins.csv", index=False)

    # ── plotting ──────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # A: |r| histogram + Cohen 기준선
    axA = axes[0]
    axA.hist(abs_r, bins=80, color="#666666", alpha=0.75, edgecolor="white")
    for thr, lbl, c in [(0.1, "trivial-weak", "#bbbbbb"),
                         (0.3, "weak-moderate", "#ff7f0e"),
                         (0.5, "moderate-strong", "#d62728"),
                         (0.7, "strong-very_strong", "#8b0000")]:
        axA.axvline(thr, color=c, lw=1.4, linestyle="--",
                     label=f"|r|={thr}  ({lbl})")
    axA.axvline(pct[0], color="black", lw=1.2,
                 label=f"median |r| = {pct[0]:.3f}")
    axA.set_xlabel("|r|  (per cell-type x gene pair)")
    axA.set_ylabel("number of pairs")
    axA.set_title(
        f"A) Distribution of |r| across all {n_pairs:,} pairs\n"
        f"median |r| = {pct[0]:.3f},  P90 = {pct[2]:.3f},  max = {abs_r.max():.3f}\n"
        f"(most pairs are weak — significance != strong association)",
        fontsize=10)
    axA.legend(loc="upper right", fontsize=7, frameon=False)

    # B: bin 별 페어 수 (전체 vs BH 통과)
    axB = axes[1]
    labels = [r["label"] for r in rows]
    counts_all = [r["n_pairs_all"] for r in rows]
    counts_bh  = [r["n_pairs_BH05"] for r in rows]
    x = np.arange(len(labels))
    w = 0.4
    bars_all = axB.bar(x - w/2, counts_all, w, label="all pairs",
                        color="#888888", edgecolor="black")
    bars_bh  = axB.bar(x + w/2, counts_bh, w, label="BH-FDR<0.05",
                        color="#d62728", edgecolor="black")
    axB.set_xticks(x)
    axB.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    axB.set_yscale("log")
    axB.set_ylabel("number of pairs (log scale)")
    for bar, cnt in zip(bars_all, counts_all):
        axB.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.05,
                  f"{cnt:,}", ha="center", va="bottom", fontsize=7)
    for bar, cnt in zip(bars_bh, counts_bh):
        if cnt > 0:
            axB.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.05,
                      f"{cnt:,}", ha="center", va="bottom", fontsize=7, color="#a02020")
    axB.set_title(
        "B) Pair count by |r| bin  -  all pairs vs BH-FDR<0.05\n"
        "(moderate/strong/very-strong counts drop sharply on log scale)",
        fontsize=10)
    axB.legend(loc="upper right", fontsize=8)
    axB.grid(axis="y", alpha=0.2)

    # C: BH-pass 페어들의 r 분포 자체
    axC = axes[2]
    r_bh = r[is_bh]
    axC.hist(r_bh, bins=60, color="#d62728", alpha=0.7, edgecolor="white")
    axC.axvline(0, c="k", lw=0.5)
    for thr in [-0.5, -0.3, 0.3, 0.5]:
        axC.axvline(thr, c="black", lw=0.8, linestyle=":")
    axC.set_xlabel("Pearson r  (BH<0.05 passing pairs only)")
    axC.set_ylabel("number of pairs")
    n_bh = is_bh.sum()
    n_bh_strong = int(((abs_r >= 0.5) & is_bh).sum())
    n_bh_mod = int(((abs_r >= 0.3) & (abs_r < 0.5) & is_bh).sum())
    n_bh_weak = int(((abs_r < 0.3) & is_bh).sum())
    axC.set_title(
        f"C) r distribution of BH<0.05 passing pairs (n={n_bh:,})\n"
        f"|r|<0.3 (weak): {n_bh_weak:,}  |  "
        f"|r| 0.3-0.5 (moderate): {n_bh_mod:,}  |  "
        f"|r|>=0.5 (strong): {n_bh_strong:,}\n"
        "(even BH-significant pairs are mostly weak-to-moderate)",
        fontsize=10)

    fig.suptitle(
        "Effect-size assessment  -  statistical significance vs association strength",
        fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(HERE / "effect_size_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[✓] effect_size_overview.png")
    print(f"\n=== HONEST ANSWER ===")
    print(f"전체 페어 {n_pairs:,} 중:")
    print(f"  median |r| = {pct[0]:.3f}, P90 |r| = {pct[2]:.3f}")
    print(f"  |r| >= 0.3 (moderate+): {int((abs_r >= 0.3).sum()):,} ({100*(abs_r >= 0.3).mean():.1f}%)")
    print(f"  |r| >= 0.5 (strong+):   {int((abs_r >= 0.5).sum()):,} ({100*(abs_r >= 0.5).mean():.1f}%)")
    print(f"  |r| >= 0.7 (very strong):{int((abs_r >= 0.7).sum()):,} ({100*(abs_r >= 0.7).mean():.1f}%)")
    print(f"BH<0.05 통과 {n_bh:,} 중:")
    print(f"  |r|<0.3 (weak): {n_bh_weak:,} ({100*n_bh_weak/n_bh:.1f}%)")
    print(f"  |r| 0.3-0.5 (moderate): {n_bh_mod:,} ({100*n_bh_mod/n_bh:.1f}%)")
    print(f"  |r|>=0.5 (strong): {n_bh_strong:,} ({100*n_bh_strong/n_bh:.1f}%)")


if __name__ == "__main__":
    main()
