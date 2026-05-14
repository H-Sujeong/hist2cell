"""Risk-score grid — axis 1/2/3 × {H2C, Pro} × {all ROIs, Tumor-only}.

본 확장의 두 가지 질문:
 (1) Tumor section 만 떼면 axis 1 ↔ risk 가 살아나는가?
     (compartment 분산 빼면 within-Tumor risk 가 axis 1 으로 보이는가)
 (2) axis 2 / axis 3 가 risk axis 인가?
     (slide1 axis 2 가 a vs b 를 잡는다는 보고와 일치하는가)

산출:
 risk_axis_grid.csv          (slide × axis × modality × subset 의 r/p)
 risk_axis_grid_summary.png  (히트맵 — |Pearson r| 그리드 한눈에)
 risk_best_axis_scatter.png  (가장 강한 신호 두 개 골라 paired scatter)
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
ROOT = HERE.parent

SECTION_GROUP = {
    "a":"Tumor-h","b":"Tumor-l","c":"Tcell-h","d":"Tcell-l","t":"Tumor-ctrl",
    "e":"Tumor-h","f":"Tumor-l","g":"Tcell-h","h":"Tcell-l","v":"Tumor-ctrl",
}
GROUP_COLOR = {
    "Tumor-h":"#d62728","Tumor-l":"#1f77b4","Tumor-ctrl":"#7f7f7f",
    "Tcell-h":"#2ca02c","Tcell-l":"#9467bd",
}
GROUP_ORDER = ["Tumor-h","Tumor-l","Tumor-ctrl","Tcell-h","Tcell-l"]
TUMOR_GROUPS = {"Tumor-h","Tumor-l","Tumor-ctrl"}
SLIDE_MARKER = {"slide1":"o","slide2":"^"}
N_PERM = 1000
SEED = 42


def load_risk(slide_id: str) -> dict:
    folder = "1_085_12" if slide_id == "slide1" else "1_152_19"
    with open(ROOT / folder / f"{folder}_ROI_groups_risk_scores.pkl", "rb") as f:
        return pickle.load(f)


def load_cca_scores(slide_id: str) -> pd.DataFrame:
    folder = "1_085_12" if slide_id == "slide1" else "1_152_19"
    df = pd.read_csv(ROOT / folder / "proof_ver2" / "cca_scores_per_roi.csv")
    df["slide"] = slide_id
    df["group"] = df["section"].map(SECTION_GROUP)
    return df


def perm_p_pearson(x, y, observed_r, n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm)
    for k in range(n_perm):
        ys = rng.permutation(y)
        null[k] = pearsonr(x, ys)[0]
    return float(np.mean(np.abs(null) >= abs(observed_r)))


def sign_align(scores: pd.Series, risk: pd.Series) -> float:
    """Return sign (+1 / -1) that makes Pearson r non-negative.
    Used so combined-slide scatter shows the natural orientation."""
    r, _ = pearsonr(scores, risk)
    return +1.0 if r >= 0 else -1.0


def main():
    out_dir = HERE
    print("[*] loading CCA scores + risk scores …")
    d1 = load_cca_scores("slide1"); r1 = load_risk("slide1")
    d2 = load_cca_scores("slide2"); r2 = load_risk("slide2")
    d1["risk"] = d1["tube_id"].map(r1)
    d2["risk"] = d2["tube_id"].map(r2)
    all_df = pd.concat([d1, d2], ignore_index=True)

    rows = []
    for axis in [1, 2, 3]:
        for mod_label, col in [("Hist2Cell", f"H2C_canon{axis}"),
                                ("Proteomics", f"P_canon{axis}")]:
            for subset_label, mask in [
                ("all_ROI",     all_df["group"].notna()),
                ("Tumor_only",  all_df["group"].isin(TUMOR_GROUPS)),
            ]:
                for slide_label, sub in [
                    ("slide1",   all_df[mask & (all_df["slide"]=="slide1")]),
                    ("slide2",   all_df[mask & (all_df["slide"]=="slide2")]),
                    ("combined", all_df[mask]),
                ]:
                    x = sub[col].values
                    y = sub["risk"].values
                    if slide_label == "combined":
                        # 두 슬라이드 합본은 부호 정렬 — 슬라이드별로 risk 와의 r 부호에 맞춰 flip
                        s1 = all_df[mask & (all_df["slide"]=="slide1")]
                        s2 = all_df[mask & (all_df["slide"]=="slide2")]
                        sign1 = sign_align(s1[col], s1["risk"])
                        sign2 = sign_align(s2[col], s2["risk"])
                        x = np.concatenate([sign1*s1[col].values, sign2*s2[col].values])
                        y = np.concatenate([s1["risk"].values, s2["risk"].values])
                    if len(x) < 5:
                        continue
                    r_p, _ = pearsonr(x, y)
                    r_s, _ = spearmanr(x, y)
                    p_perm = perm_p_pearson(x, y, r_p)
                    rows.append({
                        "axis": axis, "modality": mod_label,
                        "subset": subset_label, "slide": slide_label,
                        "n": len(x),
                        "pearson_r": r_p, "pearson_p_perm": p_perm,
                        "spearman_rho": r_s,
                    })
    grid = pd.DataFrame(rows)
    grid.to_csv(out_dir / "risk_axis_grid.csv", index=False)
    print("[*] grid table:")
    print(grid.to_string(index=False))

    # ── heatmap of |Pearson r| for combined-slide entries only ──────
    combined = grid[grid["slide"]=="combined"].copy()
    pivot = combined.pivot(index=["axis", "modality"], columns="subset",
                            values="pearson_r")
    pivot_p = combined.pivot(index=["axis", "modality"], columns="subset",
                              values="pearson_p_perm")

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.6))
    data = pivot.values
    im = ax.imshow(data, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns.tolist(), fontsize=10)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels([f"axis {a} {m}" for a, m in pivot.index.tolist()],
                        fontsize=10)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            p = pivot_p.values[i, j]
            txt_color = "white" if abs(val) > 0.35 else "black"
            star = " *" if p < 0.05 else ""
            ax.text(j, i, f"{val:+.3f}{star}", ha="center", va="center",
                     fontsize=10, color=txt_color)
    fig.colorbar(im, ax=ax, label="combined Pearson r (sign-aligned)")
    ax.set_title("Risk-score correlation grid — combined slides (n_all=94, n_Tumor=65)\n"
                  "* = permutation p < 0.05",
                  fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "risk_axis_grid_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[✓] risk_axis_grid_summary.png")

    # ── pick the strongest combined signal (largest |r| with p<0.05) ─
    best = combined.assign(abs_r=combined["pearson_r"].abs()) \
                   .sort_values("abs_r", ascending=False)
    print("\n[*] top 5 signals (combined, sign-aligned):")
    print(best[["axis","modality","subset","n","pearson_r","pearson_p_perm","spearman_rho"]]
          .head(5).to_string(index=False))

    # ── scatter for the *best Tumor-only* finding (one per modality)
    pick = best[(best.subset=="Tumor_only")].head(2)
    print("\n[*] plotting Tumor-only best:", pick[["axis","modality"]].values.tolist())

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6))
    plot_specs = pick[["axis","modality"]].values.tolist()
    if len(plot_specs) < 2:
        plot_specs = best.head(2)[["axis","modality"]].values.tolist()

    for ax_obj, (axis_id, mod_label) in zip(axes, plot_specs):
        col = f"H2C_canon{axis_id}" if mod_label == "Hist2Cell" else f"P_canon{axis_id}"
        s1 = all_df[(all_df.slide=="slide1") & (all_df.group.isin(TUMOR_GROUPS))]
        s2 = all_df[(all_df.slide=="slide2") & (all_df.group.isin(TUMOR_GROUPS))]
        sign1 = sign_align(s1[col], s1["risk"])
        sign2 = sign_align(s2[col], s2["risk"])
        for slide_id, sub, sgn in [("slide1", s1, sign1), ("slide2", s2, sign2)]:
            for g in GROUP_ORDER:
                if g not in TUMOR_GROUPS: continue
                sub_g = sub[sub.group==g]
                if sub_g.empty: continue
                ax_obj.scatter(sub_g["risk"], sgn*sub_g[col].values,
                                marker=SLIDE_MARKER[slide_id],
                                c=GROUP_COLOR[g],
                                s=72 if slide_id=="slide1" else 64,
                                edgecolors="black", linewidths=0.6, alpha=0.85)
        x = np.concatenate([sign1*s1[col].values, sign2*s2[col].values])
        y = np.concatenate([s1["risk"].values, s2["risk"].values])
        r_p, _ = pearsonr(y, x)
        coef = np.polyfit(y, x, 1)
        xs = np.linspace(y.min(), y.max(), 100)
        ax_obj.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
        ax_obj.set_xlabel("ROI risk score")
        ax_obj.set_ylabel(f"axis {axis_id} score  ({mod_label}, sign-aligned)")
        ax_obj.set_title(
            f"axis {axis_id}  {mod_label}   |   Tumor only  (n={len(x)})\n"
            f"combined Pearson r = {r_p:+.3f}",
            fontsize=10)
        ax_obj.grid(alpha=0.2)
    # build a single shared legend
    handles = []
    for g in [gx for gx in GROUP_ORDER if gx in TUMOR_GROUPS]:
        handles.append(plt.scatter([],[],marker="o",c=GROUP_COLOR[g],
                                    s=64, edgecolors="black", linewidths=0.6,
                                    label=g))
    handles.append(plt.scatter([],[],marker="o",c="white",edgecolors="black",
                                linewidths=0.6, s=64, label="slide1"))
    handles.append(plt.scatter([],[],marker="^",c="white",edgecolors="black",
                                linewidths=0.6, s=64, label="slide2"))
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9,
                frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        "Tumor-only subset  -  where does risk gradient hide?  (top 2 strongest axis x modality findings)",
        fontsize=12, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(out_dir / "risk_best_axis_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[✓] risk_best_axis_scatter.png")


if __name__ == "__main__":
    main()
