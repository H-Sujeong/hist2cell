"""Per-ROI / per-section interpretation of CCA canonical axes.

Produces 2 figures + 1 CSV per slide:
  cca_dumbbell_axis1.png       per-ROI Hist2Cell vs Proteomics canonical 1
  cca_section_means.png        section means on canonical 1/2/3, both modalities
  cca_scores_per_roi.csv       tube_id, section, H2C_canon{1,2,3}, P_canon{1,2,3}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, run_cca,
)


SECTIONS = {
    "1_085_12": {
        "labels": {"a": "High-risk Tumor", "b": "Low-risk Tumor",
                   "c": "High-risk T-cell", "d": "Low-risk T-cell",
                   "t": "Middle-risk Tumor (ctrl)"},
        "colors": {"a": "#d62728", "b": "#1f77b4", "c": "#2ca02c",
                   "d": "#9467bd", "t": "#7f7f7f"},
        "order": ["a", "b", "c", "d", "t"],
        "prefix": "abcdt",
        "pred_csv": "/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv",
        "axis1_pos_module": "epithelial / glandular  (KRT7/8/18, B_plasma_IgA, AT1/AT2, SMG_Duct)",
        "axis1_neg_module": "immune / vascular / stromal  (PTPRC=CD45, CD74, COL1A1, Muscle_smooth)",
    },
    "1_152_19": {
        "labels": {"e": "High-risk Tumor", "f": "Low-risk Tumor",
                   "g": "High-risk T-cell", "h": "Low-risk T-cell",
                   "v": "Middle-risk Tumor (ctrl)"},
        "colors": {"e": "#d62728", "f": "#1f77b4", "g": "#2ca02c",
                   "h": "#9467bd", "v": "#7f7f7f"},
        "order": ["e", "f", "g", "h", "v"],
        "prefix": "efghv",
        "pred_csv": "/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv",
        "axis1_pos_module": "vascular / smooth muscle / blood  (Muscle_smooth_*, hemoglobins, collagen)",
        "axis1_neg_module": "glandular / secretory / metabolic  (SMG_Serous/Duct, KRT7, PHGDH/PYCR1)",
    },
}


def plot_dumbbell_axis1(common, sections, Hc, Pc, slide_key, scfg, out_path):
    df = pd.DataFrame({
        "tube_id": common,
        "section": sections,
        "h": Hc[:, 0],
        "p": Pc[:, 0],
    })
    df["sec_order"] = df["section"].map(lambda s: scfg["order"].index(s))
    df = df.sort_values(["sec_order", "h"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9, max(7, 0.18 * len(df))))
    y = np.arange(len(df))
    for _, r in df.iterrows():
        col = scfg["colors"][r["section"]]
        ax.plot([r["h"], r["p"]], [y[df.index[df["tube_id"] == r["tube_id"]][0]]] * 2,
                color="#888888", linewidth=0.6, alpha=0.55, zorder=1)
    for i, r in df.iterrows():
        col = scfg["colors"][r["section"]]
        ax.scatter(r["h"], i, marker="o", s=70, color=col, edgecolor="black",
                   linewidth=0.4, zorder=2)
        ax.scatter(r["p"], i, marker="^", s=70, color=col, edgecolor="black",
                   linewidth=0.4, zorder=2, alpha=0.7)

    ax.axvline(0, color="black", linewidth=0.6, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.section} · {r.tube_id}" for _, r in df.iterrows()],
                       fontsize=6)
    ax.set_xlabel("Canonical axis 1 score", fontsize=10)
    ax.set_title(f"{slide_key} — per-ROI position on CCA axis 1\n"
                 f"○ = Hist2Cell score,  △ = Proteomics score "
                 f"(grey line = same ROI)", fontsize=11)
    handles = [plt.Rectangle((0, 0), 1, 1, color=scfg["colors"][s],
                              label=scfg["labels"][s]) for s in scfg["order"]]
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    fig.text(0.02, 0.070, f"axis 1 −  : {scfg['axis1_neg_module']}",
              fontsize=8, color="#1f77b4", ha="left", va="center")
    fig.text(0.02, 0.050, f"axis 1 +  : {scfg['axis1_pos_module']}",
              fontsize=8, color="#d62728", ha="left", va="center")
    fig.legend(handles=handles, loc="lower left",
                bbox_to_anchor=(0.02, 0.005),
                ncol=len(scfg["order"]), fontsize=7.5, framealpha=0.9)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_section_means(common, sections, Hc, Pc, train_rs, slide_key, scfg, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=False)
    width = 0.36
    x_base = np.arange(len(scfg["order"]))
    for axis_idx, ax in enumerate(axes):
        means_h = []
        std_h = []
        means_p = []
        std_p = []
        for s in scfg["order"]:
            mask = np.array(sections) == s
            means_h.append(Hc[mask, axis_idx].mean() if mask.sum() else 0)
            std_h.append(Hc[mask, axis_idx].std() if mask.sum() > 1 else 0)
            means_p.append(Pc[mask, axis_idx].mean() if mask.sum() else 0)
            std_p.append(Pc[mask, axis_idx].std() if mask.sum() > 1 else 0)
        bars_h = ax.bar(x_base - width/2, means_h, width, yerr=std_h,
                         color=[scfg["colors"][s] for s in scfg["order"]],
                         edgecolor="black", linewidth=0.5,
                         capsize=3, label="Hist2Cell")
        bars_p = ax.bar(x_base + width/2, means_p, width, yerr=std_p,
                         color=[scfg["colors"][s] for s in scfg["order"]],
                         edgecolor="black", linewidth=0.5, alpha=0.55,
                         hatch="//", capsize=3, label="Proteomics")
        ax.axhline(0, color="black", linewidth=0.5, alpha=0.4)
        ax.set_xticks(x_base)
        ax.set_xticklabels([scfg["labels"][s].replace(" ", "\n", 1)
                            for s in scfg["order"]], fontsize=7)
        ax.set_ylabel("section-mean canonical score", fontsize=8)
        ax.set_title(f"Canonical axis {axis_idx+1}  (r = {train_rs[axis_idx]:+.3f})",
                     fontsize=10)
        if axis_idx == 0:
            ax.legend(loc="best", fontsize=8)
    fig.suptitle(f"{slide_key} — section means on each CCA canonical axis "
                 f"(solid = Hist2Cell, hatched = Proteomics, bars = ±1 SD)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def run_one(slide_key):
    scfg = SECTIONS[slide_key]
    slide_dir = HERE / slide_key
    out_dir = slide_dir / "proof_ver2"
    cfg = SlideConfig(
        name=f"slide {slide_key}",
        pred_csv=Path(scfg["pred_csv"]),
        roi_pkl=slide_dir / f"{slide_key}_ROI_groups.pkl",
        npy=slide_dir / f"meteo_{slide_key}_coords.npy",
        section_label=scfg["labels"], section_color=scfg["colors"],
        sample_section_prefixes=scfg["prefix"], out_dir=out_dir,
    )
    sig_df, cell_cols = build_roi_signatures(cfg)
    log2_f, slide_cols = load_proteomics_matrix(cfg)
    common, H, P, sig_aligned, gene_index = align_modalities(
        sig_df, log2_f, slide_cols, cell_cols)
    sections = list(sig_aligned["section"])
    print(f"\n=== {slide_key} (N={len(common)}) ===")

    cca_out = run_cca(H, P)
    Hc, Pc, train_rs = cca_out["Hc"], cca_out["Pc"], cca_out["train_rs"]
    print(f"  train r per axis: {[f'{r:+.3f}' for r in train_rs]}")

    scores = pd.DataFrame({
        "tube_id": common, "section": sections,
        "H2C_canon1": Hc[:, 0], "P_canon1": Pc[:, 0],
        "H2C_canon2": Hc[:, 1], "P_canon2": Pc[:, 1],
        "H2C_canon3": Hc[:, 2], "P_canon3": Pc[:, 2],
    })
    scores.to_csv(out_dir / "cca_scores_per_roi.csv", index=False)

    plot_dumbbell_axis1(common, sections, Hc, Pc, slide_key, scfg,
                         out_dir / "cca_dumbbell_axis1.png")
    plot_section_means(common, sections, Hc, Pc, train_rs, slide_key, scfg,
                        out_dir / "cca_section_means.png")

    # quick console summary: section means on axis 1
    print(f"  axis 1 section means (H2C / Pro):")
    for s in scfg["order"]:
        mask = np.array(sections) == s
        if mask.sum() == 0: continue
        print(f"    {s} ({scfg['labels'][s]:<26}): "
              f"H2C {Hc[mask, 0].mean():+.2f}  /  Pro {Pc[mask, 0].mean():+.2f}  "
              f"(n={mask.sum()})")
    print(f"  saved → cca_dumbbell_axis1.png, cca_section_means.png, "
          f"cca_scores_per_roi.csv")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
