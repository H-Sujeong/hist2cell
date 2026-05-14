"""Negative-pair scatter — 음의 r 의 의미 시각화.

ROI 단위 (n=94) 의 (cell type abundance, gene log2 intensity) 산점도로
*음의 상관 = 공간 mutual exclusion* 의 직관 보여주기.

4 panel:
 A. B_plasma_IgA  vs  PTPRC (CD45)         r ≈ -0.71  — 면역 lineage 분리
 B. B_plasma_IgA  vs  LCP1                 r ≈ -0.71  — T-cell marker 분리
 C. SMG_Serous     vs  HBA1 (hemoglobin α) r ≈ -0.68  — 상피 vs 적혈구
 D. Mesothelia    vs  HBA1                 r ≈ -0.66  — 중피 vs 적혈구
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from _proof_ver2_lib import build_roi_signatures, load_proteomics_matrix, align_modalities
from build_joint_cca import slide_cfg

GROUP_COLOR = {
    "Tumor-h":"#d62728","Tumor-l":"#1f77b4","Tumor-ctrl":"#7f7f7f",
    "Tcell-h":"#2ca02c","Tcell-l":"#9467bd",
}
SECTION_GROUP = {
    "a":"Tumor-h","b":"Tumor-l","c":"Tcell-h","d":"Tcell-l","t":"Tumor-ctrl",
    "e":"Tumor-h","f":"Tumor-l","g":"Tcell-h","h":"Tcell-l","v":"Tumor-ctrl",
}
SLIDE_MARKER = {"slide1":"o","slide2":"^"}

PAIRS = [
    ("B_plasma_IgA", "PTPRC", "B_plasma_IgA  vsPTPRC (CD45, pan-leukocyte)"),
    ("B_plasma_IgA", "LCP1",  "B_plasma_IgA  vsLCP1 (lymphocyte cytosolic)"),
    ("SMG_Serous",   "HBA1",  "SMG_Serous   vsHBA1 (hemoglobin alpha)"),
    ("Mesothelia",   "HBA1",  "Mesothelia   vsHBA1 (hemoglobin alpha)"),
]


def build_joint():
    out = {}
    for sid in ["slide1", "slide2"]:
        cfg = slide_cfg(sid)
        sig_df, cell_cols = build_roi_signatures(cfg)
        log2_f, slide_cols = load_proteomics_matrix(cfg)
        common, H, P, sig_aligned, genes = align_modalities(
            sig_df, log2_f, slide_cols, cell_cols)
        out[sid] = {"common": common, "H": H, "P": P, "genes": genes,
                    "cell_cols": cell_cols,
                    "sections": [t[0] for t in common]}
    common_genes = sorted(set(out["slide1"]["genes"]) & set(out["slide2"]["genes"]))
    for sid in ["slide1", "slide2"]:
        idx = [out[sid]["genes"].index(g) for g in common_genes]
        out[sid]["P_common"] = out[sid]["P"][:, idx]
    H = np.vstack([StandardScaler().fit_transform(out["slide1"]["H"]),
                    StandardScaler().fit_transform(out["slide2"]["H"])])
    P = np.vstack([StandardScaler().fit_transform(out["slide1"]["P_common"]),
                    StandardScaler().fit_transform(out["slide2"]["P_common"])])
    slides = ["slide1"]*len(out["slide1"]["common"]) + ["slide2"]*len(out["slide2"]["common"])
    sections = out["slide1"]["sections"] + out["slide2"]["sections"]
    return H, P, out["slide1"]["cell_cols"], common_genes, slides, sections


def main():
    H, P, cell_cols, genes, slides, sections = build_joint()
    groups = [SECTION_GROUP[s] for s in sections]

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    axes = axes.flatten()
    for ax, (ct, gn, title) in zip(axes, PAIRS):
        i = cell_cols.index(ct)
        j = genes.index(gn)
        x = H[:, i]
        y = P[:, j]
        r, _ = pearsonr(x, y)
        for sid in ["slide1", "slide2"]:
            for g in ["Tumor-h","Tumor-l","Tumor-ctrl","Tcell-h","Tcell-l"]:
                mask = [(s == sid) and (gr == g) for s, gr in zip(slides, groups)]
                if not any(mask): continue
                idx = np.where(mask)[0]
                ax.scatter(x[idx], y[idx],
                            marker=SLIDE_MARKER[sid],
                            c=GROUP_COLOR[g],
                            s=80 if sid=="slide1" else 70,
                            edgecolors="black", linewidths=0.6, alpha=0.85,
                            label=f"{g}" if sid=="slide1" else None)
        # regression line
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
        ax.set_xlabel(f"{ct}  (Hist2Cell ROI abundance, z-scored)")
        ax.set_ylabel(f"{gn}  (log2 intensity, z-scored)")
        ax.set_title(f"{title}\nPearson r = {r:+.3f}  (n=94 ROI)", fontsize=10)
        ax.axhline(0, c="k", lw=0.4, alpha=0.3)
        ax.axvline(0, c="k", lw=0.4, alpha=0.3)
        ax.grid(alpha=0.2)
    # single legend
    handles = []
    for g in ["Tumor-h","Tumor-l","Tumor-ctrl","Tcell-h","Tcell-l"]:
        handles.append(plt.scatter([],[],c=GROUP_COLOR[g], s=70,
                                    edgecolors="black", linewidths=0.6, label=g))
    handles.append(plt.scatter([],[],c="white", marker="o", s=70,
                                edgecolors="black", linewidths=0.6, label="slide1"))
    handles.append(plt.scatter([],[],c="white", marker="^", s=70,
                                edgecolors="black", linewidths=0.6, label="slide2"))
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=9,
                frameon=False, bbox_to_anchor=(0.5, -0.005))

    fig.suptitle(
        "Negative-pair examples  -  spatial mutual exclusion across ROIs (n=94)\n"
        "(high cell-type abundance ROI → low gene intensity ROI, and vice versa)",
        fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.035, 1, 0.97])
    fig.savefig(HERE / "negative_pair_examples.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[✓] negative_pair_examples.png")


if __name__ == "__main__":
    main()
