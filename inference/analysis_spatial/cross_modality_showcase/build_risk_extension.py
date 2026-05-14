"""Risk-score extension to the cross-modality showcase.

핵심 가설: CCA axis 1 (= 두 modality 의 main coupling 축) 은
실은 *ROI 의 risk score gradient* 와 정렬되어 있다.  즉
"두 modality 가 합의해서 잡은 latent dimension = risk axis".

세 패널:
 A. Hist2Cell axis 1 (sign-aligned) vs risk score, n=94 ROI 합본
 B. Proteomics axis 1 (sign-aligned) vs risk score, n=94 ROI 합본
 C. Risk-score boxplot per section group (5 라벨) — section 라벨이
    실제 risk gradient 의 어느 위치에 해당하는지 sanity check

산출:
 risk_axis_showcase.png
 risk_axis_per_roi.csv  (94 ROI × {slide, tube_id, section, group,
                                   H_axis1, P_axis1, risk})
 risk_axis_correlations.csv  (slide-별 + 합본의 Pearson/Spearman r,
                              null permutation p)
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # inference/analysis_spatial

SECTION_GROUP = {
    "a": "Tumor-h", "b": "Tumor-l", "c": "Tcell-h", "d": "Tcell-l", "t": "Tumor-ctrl",
    "e": "Tumor-h", "f": "Tumor-l", "g": "Tcell-h", "h": "Tcell-l", "v": "Tumor-ctrl",
}
GROUP_COLOR = {
    "Tumor-h":    "#d62728",
    "Tumor-l":    "#1f77b4",
    "Tumor-ctrl": "#7f7f7f",
    "Tcell-h":    "#2ca02c",
    "Tcell-l":    "#9467bd",
}
GROUP_ORDER = ["Tumor-h", "Tumor-l", "Tumor-ctrl", "Tcell-h", "Tcell-l"]
SLIDE_MARKER = {"slide1": "o", "slide2": "^"}

N_PERM = 1000
SEED = 42


def load_risk(slide_id: str) -> dict:
    folder = "1_085_12" if slide_id == "slide1" else "1_152_19"
    pkl = ROOT / folder / f"{folder}_ROI_groups_risk_scores.pkl"
    with open(pkl, "rb") as f:
        return pickle.load(f)


def perm_p(x, y, observed_r, n_perm=N_PERM, seed=SEED):
    """Two-sided permutation p for Pearson r."""
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm)
    for k in range(n_perm):
        ys = rng.permutation(y)
        null[k] = pearsonr(x, ys)[0]
    return float(np.mean(np.abs(null) >= abs(observed_r))), null


def main():
    out_dir = HERE
    print("[*] loading axis 1 paired (sign-aligned) + risk scores …")
    axis = pd.read_csv(out_dir / "axis1_paired_combined.csv")
    risk1 = load_risk("slide1")
    risk2 = load_risk("slide2")
    axis["risk"] = axis.apply(
        lambda r: risk1[r["tube_id"]] if r["slide"] == "slide1" else risk2[r["tube_id"]],
        axis=1,
    )
    axis.to_csv(out_dir / "risk_axis_per_roi.csv", index=False)

    # ── per-slide + combined correlations ────────────────────────────
    rows = []
    for label, df in [
        ("slide1", axis[axis["slide"] == "slide1"]),
        ("slide2", axis[axis["slide"] == "slide2"]),
        ("combined (n=94)", axis),
    ]:
        for mod_label, col in [("Hist2Cell axis 1", "H_axis1"),
                                ("Proteomics axis 1", "P_axis1")]:
            x = df[col].values
            y = df["risk"].values
            r_p, p_p_param = pearsonr(x, y)
            r_s, p_s_param = spearmanr(x, y)
            p_perm_two, _ = perm_p(x, y, r_p)
            rows.append({
                "scope": label, "modality": mod_label,
                "n": len(df),
                "pearson_r": r_p, "pearson_p_param": p_p_param,
                "pearson_p_perm_2sided": p_perm_two,
                "spearman_rho": r_s, "spearman_p_param": p_s_param,
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "risk_axis_correlations.csv", index=False)
    print("[*] correlations:")
    print(summary[["scope", "modality", "n", "pearson_r", "pearson_p_perm_2sided",
                   "spearman_rho"]].to_string(index=False))

    # extract main numbers for plot titles
    def get(scope, mod):
        return summary[(summary.scope == scope) & (summary.modality == mod)].iloc[0]
    h1 = get("slide1", "Hist2Cell axis 1")
    h2 = get("slide2", "Hist2Cell axis 1")
    hc = get("combined (n=94)", "Hist2Cell axis 1")
    p1 = get("slide1", "Proteomics axis 1")
    p2 = get("slide2", "Proteomics axis 1")
    pc = get("combined (n=94)", "Proteomics axis 1")

    # ── plotting ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(19.5, 6.4))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[1.0, 1.0, 0.95], wspace=0.32,
        left=0.055, right=0.985, top=0.85, bottom=0.12,
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    # ── Panel A: H axis 1 vs risk ────────────────────────────────────
    for slide_id in ["slide1", "slide2"]:
        for g in GROUP_ORDER:
            sub = axis[(axis.slide == slide_id) & (axis.group == g)]
            if sub.empty: continue
            axA.scatter(sub["risk"], sub["H_axis1"],
                        marker=SLIDE_MARKER[slide_id],
                        c=GROUP_COLOR[g], s=72 if slide_id == "slide1" else 64,
                        edgecolors="black", linewidths=0.6, alpha=0.85,
                        label=f"{slide_id} {g}" if slide_id == "slide1" else None)
    # combined regression line
    x = axis["risk"].values; y = axis["H_axis1"].values
    coef = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    axA.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
    axA.axhline(0, c="k", lw=0.4, alpha=0.3)
    axA.set_xlabel("ROI risk score")
    axA.set_ylabel("Hist2Cell CCA axis-1 score   (+ = epithelial/glandular)")
    axA.set_title(
        "A) Hist2Cell axis 1  vs  risk score\n"
        f"slide1 r = {h1.pearson_r:+.3f} (p_perm={h1.pearson_p_perm_2sided:.3f})   "
        f"slide2 r = {h2.pearson_r:+.3f} (p_perm={h2.pearson_p_perm_2sided:.3f})\n"
        f"combined (n=94)  Pearson r = {hc.pearson_r:+.3f}   "
        f"Spearman rho = {hc.spearman_rho:+.3f}",
        fontsize=10)
    axA.legend(loc="upper right", fontsize=7, frameon=False, ncol=1)
    axA.grid(alpha=0.2)

    # ── Panel B: P axis 1 vs risk ────────────────────────────────────
    for slide_id in ["slide1", "slide2"]:
        for g in GROUP_ORDER:
            sub = axis[(axis.slide == slide_id) & (axis.group == g)]
            if sub.empty: continue
            axB.scatter(sub["risk"], sub["P_axis1"],
                        marker=SLIDE_MARKER[slide_id],
                        c=GROUP_COLOR[g], s=72 if slide_id == "slide1" else 64,
                        edgecolors="black", linewidths=0.6, alpha=0.85)
    x = axis["risk"].values; y = axis["P_axis1"].values
    coef = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    axB.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
    axB.axhline(0, c="k", lw=0.4, alpha=0.3)
    axB.set_xlabel("ROI risk score")
    axB.set_ylabel("Proteomics CCA axis-1 score")
    axB.set_title(
        "B) Proteomics axis 1  vs  risk score\n"
        f"slide1 r = {p1.pearson_r:+.3f} (p_perm={p1.pearson_p_perm_2sided:.3f})   "
        f"slide2 r = {p2.pearson_r:+.3f} (p_perm={p2.pearson_p_perm_2sided:.3f})\n"
        f"combined (n=94)  Pearson r = {pc.pearson_r:+.3f}   "
        f"Spearman rho = {pc.spearman_rho:+.3f}",
        fontsize=10)
    axB.grid(alpha=0.2)

    # ── Panel C: risk score boxplot per group ─────────────────────────
    box_data = []
    box_labels = []
    box_colors = []
    for g in GROUP_ORDER:
        sub = axis[axis.group == g]
        if sub.empty: continue
        box_data.append(sub["risk"].values)
        box_labels.append(f"{g}\n(n={len(sub)})")
        box_colors.append(GROUP_COLOR[g])
    bp = axC.boxplot(box_data, labels=box_labels, patch_artist=True,
                      widths=0.55, showfliers=False)
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("black")
    for median in bp["medians"]:
        median.set_color("black"); median.set_linewidth(1.4)
    # overlay individual ROI dots
    for i, (data, c, g) in enumerate(zip(box_data, box_colors, GROUP_ORDER), start=1):
        sub = axis[axis.group == g]
        for slide_id in ["slide1", "slide2"]:
            sub_s = sub[sub.slide == slide_id]
            if sub_s.empty: continue
            xs = np.random.default_rng(SEED + i).normal(i, 0.08, size=len(sub_s))
            axC.scatter(xs, sub_s["risk"].values,
                        marker=SLIDE_MARKER[slide_id], s=36, c=c,
                        edgecolors="black", linewidths=0.4, alpha=0.85)
    axC.axhline(axis["risk"].mean(), color="k", lw=0.6, linestyle=":", alpha=0.5,
                 label=f"overall mean = {axis['risk'].mean():.2f}")
    axC.set_ylabel("ROI risk score")
    # Kruskal-Wallis test across 5 groups
    from scipy.stats import kruskal
    H_kw, p_kw = kruskal(*box_data)
    axC.set_title(
        "C) Risk-score distribution per section group\n"
        f"Kruskal-Wallis H = {H_kw:.2f}   p = {p_kw:.2e}\n"
        "section group label is a 5-bin discretization of underlying risk",
        fontsize=10)
    axC.legend(loc="upper right", fontsize=8, frameon=False)
    axC.grid(axis="y", alpha=0.2)

    fig.suptitle(
        "CCA axis 1  =  Risk-score axis  -  the two-modality coupling tracks an actual biological gradient",
        fontsize=12, fontweight="bold", y=0.97)
    out_png = out_dir / "risk_axis_showcase.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[✓] {out_png}")
    print(f"    combined H axis 1 vs risk  r = {hc.pearson_r:+.3f}  (Spearman {hc.spearman_rho:+.3f})")
    print(f"    combined P axis 1 vs risk  r = {pc.pearson_r:+.3f}  (Spearman {pc.spearman_rho:+.3f})")
    print(f"    Kruskal-Wallis across 5 groups:  H = {H_kw:.2f}  p = {p_kw:.2e}")


if __name__ == "__main__":
    main()
