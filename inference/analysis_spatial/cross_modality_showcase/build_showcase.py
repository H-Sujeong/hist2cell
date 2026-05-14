"""Cross-modality showcase — 3 panel single-page figure.

목적: Hist2Cell × proteomics 의 *양의 상관관계* 가 두 슬라이드에서
독립적으로 재현된다는 사실을 한 그림에서 명쾌히 보여주기.

세 패널:
 A. Mantel scatter — ROI×ROI 거리 (PCA10 Euclidean) 두 슬라이드 합본
 B. CCA axis 1 paired scatter — 두 슬라이드 부호 정렬 합본 (n=94 ROI)
 C. 사전 등록 8 가설 forest plot — slide1·slide2 의 effect size + match/mismatch

산출:
 cross_modality_showcase.png (3 panel)
 mantel_combined.csv          (slide-별 + 합본 Mantel r/p)
 axis1_paired_combined.csv    (94 ROI × {slide, section, H_axis1, P_axis1, sign_flipped})
 forest_hypotheses.csv        (16 행 = 8 hyp × 2 slide)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # inference/analysis_spatial
sys.path.insert(0, str(ROOT))

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, N_PCS, RANDOM_SEED,
)

SECTION_LABEL_S1 = {"a":"Tumor-h", "b":"Tumor-l", "c":"Tcell-h", "d":"Tcell-l", "t":"Tumor-ctrl"}
SECTION_LABEL_S2 = {"e":"Tumor-h", "f":"Tumor-l", "g":"Tcell-h", "h":"Tcell-l", "v":"Tumor-ctrl"}
SECTION_LABEL = {**SECTION_LABEL_S1, **SECTION_LABEL_S2}

# 두 슬라이드 *같은 의미의 section* 에 같은 색을 줘서 합본 scatter 에서 한눈에 매칭됨
GROUP_COLOR = {
    "Tumor-h":    "#d62728",  # 빨강
    "Tumor-l":    "#1f77b4",  # 파랑
    "Tcell-h":    "#2ca02c",  # 초록
    "Tcell-l":    "#9467bd",  # 보라
    "Tumor-ctrl": "#7f7f7f",  # 회색
}
GROUP_ORDER = ["Tumor-h", "Tumor-l", "Tumor-ctrl", "Tcell-h", "Tcell-l"]

SLIDE_MARKER = {"slide1": "o", "slide2": "^"}

N_PERM = 1000
rng_global = np.random.default_rng(RANDOM_SEED)


# ───────────────────────── slide config ───────────────────────────────

def slide_cfg(slide_id: str) -> SlideConfig:
    if slide_id == "slide1":
        sec_lbl = {"a":"Tumor-h","b":"Tumor-l","c":"Tcell-h","d":"Tcell-l","t":"Tumor-ctrl"}
        sec_col = {"a":"#d62728","b":"#1f77b4","c":"#2ca02c","d":"#9467bd","t":"#7f7f7f"}
        return SlideConfig(
            name="slide1 (1_085_12)",
            pred_csv=Path("/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv"),
            roi_pkl=ROOT / "1_085_12" / "1_085_12_ROI_groups.pkl",
            npy=ROOT / "1_085_12" / "meteo_1_085_12_coords.npy",
            section_label=sec_lbl, section_color=sec_col,
            sample_section_prefixes="abcdt",
            out_dir=ROOT / "1_085_12" / "proof_ver2",
        )
    else:
        sec_lbl = {"e":"Tumor-h","f":"Tumor-l","g":"Tcell-h","h":"Tcell-l","v":"Tumor-ctrl"}
        sec_col = {"e":"#d62728","f":"#1f77b4","g":"#2ca02c","h":"#9467bd","v":"#7f7f7f"}
        return SlideConfig(
            name="slide2 (1_152_19)",
            pred_csv=Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv"),
            roi_pkl=ROOT / "1_152_19" / "1_152_19_ROI_groups.pkl",
            npy=ROOT / "1_152_19" / "meteo_1_152_19_coords.npy",
            section_label=sec_lbl, section_color=sec_col,
            sample_section_prefixes="efghv",
            out_dir=ROOT / "1_152_19" / "proof_ver2",
        )


# ───────────────────────── per-slide distance matrices ────────────────

def slide_distance_pairs(slide_id: str):
    """Return per-ROI section labels and the two upper-triangle distance
    vectors (Hist2Cell, Proteomics) computed in PCA-10 z-scored Euclidean
    space — same reduction used by the CCA pipeline."""
    cfg = slide_cfg(slide_id)
    sig_df, cell_cols = build_roi_signatures(cfg)
    log2_f, slide_cols = load_proteomics_matrix(cfg)
    common, H, P, sig_aligned, gene_index = align_modalities(
        sig_df, log2_f, slide_cols, cell_cols)
    sections = [t[0] for t in common]
    group_lbl = [SECTION_LABEL[s] for s in sections]

    Hs = StandardScaler().fit_transform(H)
    Ps = StandardScaler().fit_transform(P)
    n_pcs = min(N_PCS, H.shape[0] - 1, H.shape[1], P.shape[1])
    H_pc = PCA(n_components=n_pcs, random_state=RANDOM_SEED).fit_transform(Hs)
    P_pc = PCA(n_components=n_pcs, random_state=RANDOM_SEED).fit_transform(Ps)

    d_h = pdist(H_pc, metric="euclidean")
    d_p = pdist(P_pc, metric="euclidean")
    n = H.shape[0]
    # pair labels
    pair_same_section = []
    pair_same_group   = []
    for i in range(n):
        for j in range(i + 1, n):
            pair_same_section.append(sections[i] == sections[j])
            pair_same_group.append(group_lbl[i] == group_lbl[j])
    return {
        "n": n, "sections": sections, "groups": group_lbl,
        "d_h": d_h, "d_p": d_p,
        "pair_same_section": np.array(pair_same_section),
        "pair_same_group":   np.array(pair_same_group),
        "common_ids":  common,
    }


def mantel_perm_p(d_h, d_p, n, n_perm=N_PERM, seed=RANDOM_SEED):
    """Permutation test on the rows/columns of one distance matrix —
    correctly preserves the ROI structure."""
    rng = np.random.default_rng(seed)
    Dh = squareform(d_h)
    Dp = squareform(d_p)
    obs = pearsonr(d_h, d_p)[0]
    null = np.zeros(n_perm)
    for it in range(n_perm):
        perm = rng.permutation(n)
        Dp_perm = Dp[perm][:, perm]
        d_p_perm = squareform(Dp_perm, checks=False)
        null[it] = pearsonr(d_h, d_p_perm)[0]
    p = float((np.abs(null) >= abs(obs)).mean())
    return obs, null, p


# ───────────────────────── axis 1 sign alignment ──────────────────────

def load_axis1(slide_id: str):
    """Load per-ROI CCA axis 1 scores, returning a DataFrame with
    columns {slide, tube_id, section, group, H_axis1, P_axis1}.

    Sign convention: + direction = epithelial/glandular module
    (slide1 axis1 is already that orientation; slide2 needs sign flip).
    """
    if slide_id == "slide1":
        csv = ROOT / "1_085_12" / "proof_ver2" / "cca_scores_per_roi.csv"
        sign = +1
    else:
        csv = ROOT / "1_152_19" / "proof_ver2" / "cca_scores_per_roi.csv"
        sign = -1
    df = pd.read_csv(csv)
    df["slide"] = slide_id
    df["H_axis1"] = sign * df["H2C_canon1"]
    df["P_axis1"] = sign * df["P_canon1"]
    df["group"] = df["section"].map(SECTION_LABEL)
    return df[["slide", "tube_id", "section", "group", "H_axis1", "P_axis1"]]


# ───────────────────────── forest plot data ───────────────────────────

def load_hypotheses_forest():
    s1 = pd.read_csv(ROOT / "1_085_12" / "cell_typing" / "proteomics_marker_hypotheses.csv")
    s2 = pd.read_csv(ROOT / "1_152_19" / "cell_typing" / "marker_hypotheses.csv")
    s1["slide"] = "slide1"; s2["slide"] = "slide2"
    df = pd.concat([s1, s2], ignore_index=True)
    df["hypothesis"] = df["protein_marker"] + " ↔ " + df["hist2cell_type"]
    df = df[["slide", "hypothesis", "protein_marker", "hist2cell_type",
             "predicted_direction", "observed_direction",
             "matches_hypothesis", "delta", "p_bh"]]
    # 통일된 가설 순서 (8개)
    hyp_order = (
        s1[["protein_marker", "hist2cell_type"]]
        .drop_duplicates()
        .assign(hyp=lambda d: d["protein_marker"] + " ↔ " + d["hist2cell_type"])
        ["hyp"].tolist()
    )
    df["hyp_order"] = df["hypothesis"].map({h: i for i, h in enumerate(hyp_order)})
    return df.sort_values(["hyp_order", "slide"]).reset_index(drop=True), hyp_order


# ───────────────────────── 3-panel figure ─────────────────────────────

def main():
    out_dir = HERE
    out_dir.mkdir(exist_ok=True)

    print("[*] loading slide1 / slide2 distance matrices …")
    d1 = slide_distance_pairs("slide1")
    d2 = slide_distance_pairs("slide2")

    print("[*] Mantel per slide (1000 perm) …")
    obs_r1, null_r1, p_r1 = mantel_perm_p(d1["d_h"], d1["d_p"], d1["n"])
    obs_r2, null_r2, p_r2 = mantel_perm_p(d2["d_h"], d2["d_p"], d2["n"])

    # 합본 Mantel — slide-별 z-score 후 vertical stack (스케일 다름 보정)
    zh = np.concatenate([
        (d1["d_h"] - d1["d_h"].mean()) / d1["d_h"].std(),
        (d2["d_h"] - d2["d_h"].mean()) / d2["d_h"].std(),
    ])
    zp = np.concatenate([
        (d1["d_p"] - d1["d_p"].mean()) / d1["d_p"].std(),
        (d2["d_p"] - d2["d_p"].mean()) / d2["d_p"].std(),
    ])
    same_group = np.concatenate([d1["pair_same_group"], d2["pair_same_group"]])
    slide_tag  = np.array(["slide1"] * len(d1["d_h"]) + ["slide2"] * len(d2["d_h"]))
    r_combined, _ = pearsonr(zh, zp)
    rho_combined, _ = spearmanr(zh, zp)

    pd.DataFrame({
        "slide": ["slide1", "slide2", "combined (z-score stacked)"],
        "n_pairs": [len(d1["d_h"]), len(d2["d_h"]), len(zh)],
        "pearson_r": [obs_r1, obs_r2, r_combined],
        "pearson_p_perm": [p_r1, p_r2, np.nan],
        "spearman_rho": [
            spearmanr(d1["d_h"], d1["d_p"])[0],
            spearmanr(d2["d_h"], d2["d_p"])[0],
            rho_combined,
        ],
    }).to_csv(out_dir / "mantel_combined.csv", index=False)

    print("[*] axis 1 paired scores (sign-aligned) …")
    a1 = load_axis1("slide1")
    a2 = load_axis1("slide2")
    axis = pd.concat([a1, a2], ignore_index=True)
    axis.to_csv(out_dir / "axis1_paired_combined.csv", index=False)
    r_axis_s1 = pearsonr(a1["H_axis1"], a1["P_axis1"])[0]
    r_axis_s2 = pearsonr(a2["H_axis1"], a2["P_axis1"])[0]
    r_axis_all = pearsonr(axis["H_axis1"], axis["P_axis1"])[0]

    print("[*] hypothesis forest data …")
    forest, hyp_order = load_hypotheses_forest()
    forest.to_csv(out_dir / "forest_hypotheses.csv", index=False)

    # ───────── plotting ─────────
    print("[*] building 3-panel figure …")
    fig = plt.figure(figsize=(19.5, 6.4))
    gs = fig.add_gridspec(
        1, 3,
        width_ratios=[1.0, 1.0, 1.25],
        wspace=0.36,
        left=0.055, right=0.985, top=0.85, bottom=0.11,
    )

    # ── Panel A: Mantel scatter combined
    axA = fig.add_subplot(gs[0, 0])
    diff_mask = ~same_group
    axA.scatter(zh[diff_mask], zp[diff_mask], s=10, c="#cccccc",
                alpha=0.4, edgecolors="none", label="different-group ROI pairs")
    axA.scatter(zh[same_group], zp[same_group], s=18, c="#ff7f0e",
                alpha=0.65, edgecolors="none", label="same-group ROI pairs")
    coef = np.polyfit(zh, zp, 1)
    xline = np.linspace(zh.min(), zh.max(), 100)
    axA.plot(xline, np.polyval(coef, xline), "k-", lw=1.6, label=f"regression (slope={coef[0]:.2f})")
    axA.axhline(0, c="k", lw=0.4, alpha=0.3)
    axA.axvline(0, c="k", lw=0.4, alpha=0.3)
    axA.set_xlabel("Hist2Cell ROI-pair distance (slide-z-scored)")
    axA.set_ylabel("Proteomics ROI-pair distance (slide-z-scored)")
    axA.set_title(
        "A) Mantel combined  -  ROI x ROI distance, both modalities\n"
        f"slide1 r=+{obs_r1:.3f} (p={p_r1:.3f})   slide2 r=+{obs_r2:.3f} (p={p_r2:.3f})\n"
        f"combined (n={len(zh)} pairs)  Pearson r=+{r_combined:.3f}   Spearman rho=+{rho_combined:.3f}",
        fontsize=10)
    axA.legend(loc="upper left", fontsize=8, frameon=False)
    axA.grid(alpha=0.2)

    # ── Panel B: axis 1 paired scatter (sign-aligned, n=94)
    axB = fig.add_subplot(gs[0, 1])
    for slide_id in ["slide1", "slide2"]:
        for g in GROUP_ORDER:
            sub = axis[(axis["slide"] == slide_id) & (axis["group"] == g)]
            if sub.empty:
                continue
            axB.scatter(sub["H_axis1"], sub["P_axis1"],
                        s=72 if slide_id == "slide1" else 64,
                        c=GROUP_COLOR[g],
                        marker=SLIDE_MARKER[slide_id],
                        edgecolors="black", linewidths=0.6,
                        alpha=0.85,
                        label=f"{slide_id} {g}" if slide_id == "slide1" else None)
    lo = min(axis["H_axis1"].min(), axis["P_axis1"].min()) - 0.3
    hi = max(axis["H_axis1"].max(), axis["P_axis1"].max()) + 0.3
    axB.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.5, label="y=x")
    axB.set_xlim(lo, hi); axB.set_ylim(lo, hi)
    axB.set_xlabel("Hist2Cell CCA axis-1 score   (+ = epithelial/glandular)")
    axB.set_ylabel("Proteomics CCA axis-1 score")
    axB.set_title(
        "B) CCA axis 1 paired scores  -  two slides sign-aligned, combined\n"
        f"slide1 Pearson r=+{r_axis_s1:.3f}   slide2 Pearson r=+{r_axis_s2:.3f}\n"
        f"combined n=94 ROI  Pearson r=+{r_axis_all:.3f}   (o = slide1,  ^ = slide2)",
        fontsize=10)
    leg = axB.legend(loc="upper left", fontsize=7, frameon=False, ncol=1)
    axB.grid(alpha=0.2)

    # ── Panel C: 8/8 forest — both slides
    axC = fig.add_subplot(gs[0, 2])
    n_hyp = len(hyp_order)
    y_positions = np.arange(n_hyp)[::-1]
    for _, row in forest.iterrows():
        y = y_positions[row["hyp_order"]]
        if row["slide"] == "slide1":
            y_off = +0.18; mk = "o"
        else:
            y_off = -0.18; mk = "^"
        clr = "#2ca02c" if row["matches_hypothesis"] else "#d62728"
        # BH significance star (ASCII only — Korean fonts not installed)
        star = "*" if row["p_bh"] < 0.01 else ("+" if row["p_bh"] < 0.05 else "")
        axC.scatter(row["delta"], y + y_off, marker=mk, s=110, c=clr,
                    edgecolors="black", linewidths=0.7, zorder=3)
        if star:
            axC.text(row["delta"], y + y_off, star, ha="center", va="center",
                     fontsize=7, color="white", zorder=4, fontweight="bold")
    axC.axvline(0, color="k", lw=0.8, alpha=0.5)
    axC.set_yticks(y_positions)
    short_lbl = []
    for h in hyp_order:
        marker, ctype = h.split(" ↔ ")
        marker = marker.replace(" (mitosis)", "").replace(" (smooth muscle)", "")
        short_lbl.append(f"{marker} -> {ctype}")
    axC.yaxis.tick_right()
    axC.yaxis.set_label_position("right")
    axC.set_yticklabels(short_lbl, fontsize=8.5, family="monospace")
    axC.tick_params(axis="y", which="major", pad=4)
    axC.set_xlabel("Hist2Cell delta  (High-risk - Low-risk Tumor)")
    axC.set_title(
        "C) 8 pre-registered hypotheses  -  slide1 8/8 match, slide2 5/8 match\n"
        "o = slide1,  ^ = slide2.   green = predicted direction match,  red = opposite\n"
        "*  BH-FDR < 0.01,   +  BH-FDR < 0.05",
        fontsize=10)
    axC.grid(axis="x", alpha=0.2)
    axC.set_xlim(-0.4, max(forest["delta"].max() + 0.2, 1.8))

    fig.suptitle(
        "Hist2Cell  x  Proteomics  -  positive cross-modality coupling reproduced across both slides",
        fontsize=12, fontweight="bold", y=0.97)
    out_png = out_dir / "cross_modality_showcase.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[✓] {out_png}")
    print(f"    합본 Mantel Pearson r = +{r_combined:.3f}")
    print(f"    합본 axis 1 Pearson r = +{r_axis_all:.3f}  (n=94 ROI)")
    print(f"    forest: slide1 = {forest[forest.slide=='slide1'].matches_hypothesis.sum()}/8 일치, "
          f"slide2 = {forest[forest.slide=='slide2'].matches_hypothesis.sum()}/8 일치")


if __name__ == "__main__":
    main()
