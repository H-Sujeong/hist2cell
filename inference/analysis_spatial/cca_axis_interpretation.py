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
        "axis_modules": {
            1: ("epithelial / glandular  (KRT7/8/18, B_plasma_IgA, AT1/AT2, SMG_Duct)",
                "immune / vascular / stromal  (PTPRC=CD45, CD74, COL1A1, Muscle_smooth)"),
            2: ("ciliated alveolar epithelium / capillary  (Ciliated, AT1/2, Fibro_alveolar, HBA1/B, CNN1)",
                "adventitial fibroblast / glandular secretory  (Fibro_adventitial, SMG_Serous, SDC1, PTX3, CALML5)"),
            3: ("plasma B-cell / mucous-secretory glandular  (SMG_Duct, B_plasma_IgA, Goblet, S100P/A6, PTX3)",
                "alveolar epithelium / basal stress keratin  (AT1/2, Muscle_smooth, KRT6A/B/72/75, TAGLN)"),
        },
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
        "axis_modules": {
            1: ("vascular / smooth muscle / blood  (Muscle_smooth_*, hemoglobins, collagen)",
                "glandular / secretory / metabolic  (SMG_Serous/Duct, KRT7, PHGDH/PYCR1)"),
            2: ("cornified squamous / plasma-B  (Chondrocyte, B_plasmablast, SPRR2G/1A, CEACAM1, BPIFA2)",
                "alveolar / airway smooth muscle / immunoglobulin  (AT2, Muscle_airway, Ciliated, IGKV*, PRR4)"),
            3: ("keratinized epithelium / mucinous-S100  (Chondrocyte, Ciliated, S100A7/A7A, MUCL1, CYP1B1)",
                "stromal fibroblast / ECM collagen  (Fibro_alveolar, COL12A1, COL1A2, BGN, THBS1, TAGLN)"),
        },
    },
}


def plot_dumbbell(common, sections, Hc, Pc, slide_key, scfg, out_path,
                  axis_idx=0):
    """Dumbbell plot for any canonical axis (axis_idx 0/1/2)."""
    axis_num = axis_idx + 1
    pos_mod, neg_mod = scfg["axis_modules"][axis_num]
    df = pd.DataFrame({
        "tube_id": common,
        "section": sections,
        "h": Hc[:, axis_idx],
        "p": Pc[:, axis_idx],
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
    ax.set_xlabel(f"Canonical axis {axis_num} score", fontsize=10)
    ax.set_title(f"{slide_key} — per-ROI position on CCA axis {axis_num}\n"
                 f"○ = Hist2Cell score,  △ = Proteomics score "
                 f"(grey line = same ROI)", fontsize=11)
    handles = [plt.Rectangle((0, 0), 1, 1, color=scfg["colors"][s],
                              label=scfg["labels"][s]) for s in scfg["order"]]
    fig.tight_layout(rect=[0, 0.09, 1, 1])
    fig.text(0.02, 0.070, f"axis {axis_num} −  : {neg_mod}",
              fontsize=8, color="#1f77b4", ha="left", va="center")
    fig.text(0.02, 0.050, f"axis {axis_num} +  : {pos_mod}",
              fontsize=8, color="#d62728", ha="left", va="center")
    fig.legend(handles=handles, loc="lower left",
                bbox_to_anchor=(0.02, 0.005),
                ncol=len(scfg["order"]), fontsize=7.5, framealpha=0.9)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_loadings(h_load, p_load, cell_cols, gene_index, axis_idx, slide_key,
                  out_path, top_n=12):
    h_axis = pd.Series(h_load[:, axis_idx], index=cell_cols)
    p_axis = pd.Series(p_load[:, axis_idx], index=gene_index)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for ax, (data, title, color) in zip(
            axes.flat,
            [(h_axis.nlargest(top_n),  f"Hist2Cell top + loaders (axis {axis_idx+1})", "#d62728"),
             (h_axis.nsmallest(top_n), f"Hist2Cell top − loaders (axis {axis_idx+1})", "#1f77b4"),
             (p_axis.nlargest(top_n),  f"Proteomics top + loaders (axis {axis_idx+1})", "#d62728"),
             (p_axis.nsmallest(top_n), f"Proteomics top − loaders (axis {axis_idx+1})", "#1f77b4")]):
        data = data.sort_values()
        ax.barh(range(len(data)), data.values, color=color, alpha=0.75,
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(data)))
        ax.set_yticklabels(data.index, fontsize=8)
        ax.set_xlabel("loading on canonical axis")
        ax.set_title(title, fontsize=10)
        ax.axvline(0, color="black", linewidth=0.4)
    fig.suptitle(f"CCA loadings — {slide_key} canonical axis {axis_idx+1}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def axis_intra_inter(scores_1d, sections, n_perm=1000, seed=0):
    """One canonical axis at a time: are within-section |Δ score| pairs
    smaller than between-section pairs? Returns intra/inter ratio +
    permutation-null comparison."""
    sec = np.array(sections)
    n = len(sec)
    s = np.asarray(scores_1d, dtype=float)
    diff = np.abs(s[:, None] - s[None, :])
    same = sec[:, None] == sec[None, :]
    eye = np.eye(n, dtype=bool)
    intra = diff[same & ~eye]
    inter = diff[~same & ~eye]
    intra_mean = float(intra.mean())
    inter_mean = float(inter.mean())
    ratio = float(intra_mean / inter_mean)
    rng = np.random.default_rng(seed)
    null_ratio = np.zeros(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(sec)
        same_p = perm[:, None] == perm[None, :]
        same_p &= ~eye
        intra_p = diff[same_p]
        inter_p = diff[~same_p & ~eye]
        null_ratio[k] = intra_p.mean() / inter_p.mean()
    p_perm = float(np.mean(null_ratio <= ratio))
    return {
        "intra_mean": intra_mean,
        "inter_mean": inter_mean,
        "ratio": ratio,
        "null_ratio_mean": float(null_ratio.mean()),
        "null_ratio_95lo": float(np.percentile(null_ratio, 2.5)),
        "null_ratio_95hi": float(np.percentile(null_ratio, 97.5)),
        "p_perm": p_perm,
    }


def plot_axis_intra_inter(rows, slide_key, out_path):
    """rows: list of dicts with axis, modality, ratio, null_ratio_mean,
    null_ratio_95lo, null_ratio_95hi."""
    df = pd.DataFrame(rows)
    axes_uniq = sorted(df["axis"].unique())
    mods = ["Hist2Cell", "Proteomics"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    width = 0.35
    x_base = np.arange(len(axes_uniq))
    for m_idx, mod in enumerate(mods):
        sub = df[df["modality"] == mod].set_index("axis").loc[axes_uniq]
        offset = (m_idx - 0.5) * width
        bar = ax.bar(x_base + offset, sub["ratio"], width,
                     color=("#d62728" if mod == "Hist2Cell" else "#1f77b4"),
                     edgecolor="black", linewidth=0.5,
                     alpha=0.85 if mod == "Hist2Cell" else 0.55,
                     label=mod, hatch=("" if mod == "Hist2Cell" else "//"))
        # null 95% CI as error bars going up
        null_lo = sub["null_ratio_95lo"].values
        null_hi = sub["null_ratio_95hi"].values
        for xi, (lo, hi) in enumerate(zip(null_lo, null_hi)):
            ax.plot([x_base[xi] + offset - width*0.4,
                     x_base[xi] + offset + width*0.4],
                    [lo, lo], color="gray", linewidth=1)
            ax.plot([x_base[xi] + offset - width*0.4,
                     x_base[xi] + offset + width*0.4],
                    [hi, hi], color="gray", linewidth=1)
            ax.plot([x_base[xi] + offset]*2, [lo, hi],
                    color="gray", linewidth=1, alpha=0.7)
    ax.axhline(1.0, color="black", linewidth=0.5, alpha=0.4)
    ax.set_xticks(x_base)
    ax.set_xticklabels([f"axis {a}" for a in axes_uniq], fontsize=10)
    ax.set_ylabel("intra-section / inter-section  mean |Δ score|", fontsize=9)
    ax.set_title(f"{slide_key} — per-axis intra/inter section separation\n"
                 f"(bar = observed, grey bracket = permutation null 95% CI;  "
                 f"smaller = section structure clearer)", fontsize=10)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def print_top_loaders(h_load, p_load, cell_cols, gene_index, axis_idx, top_n=8):
    h_axis = pd.Series(h_load[:, axis_idx], index=cell_cols)
    p_axis = pd.Series(p_load[:, axis_idx], index=gene_index)
    print(f"  axis {axis_idx+1} top + Hist2Cell:  "
          f"{', '.join(f'{c}({v:+.2f})' for c, v in h_axis.nlargest(top_n).items())}")
    print(f"  axis {axis_idx+1} top − Hist2Cell:  "
          f"{', '.join(f'{c}({v:+.2f})' for c, v in h_axis.nsmallest(top_n).items())}")
    print(f"  axis {axis_idx+1} top + Proteomics: "
          f"{', '.join(f'{c}({v:+.3f})' for c, v in p_axis.nlargest(top_n).items())}")
    print(f"  axis {axis_idx+1} top − Proteomics: "
          f"{', '.join(f'{c}({v:+.3f})' for c, v in p_axis.nsmallest(top_n).items())}")


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

    for axis_idx in range(3):
        plot_dumbbell(common, sections, Hc, Pc, slide_key, scfg,
                       out_dir / f"cca_dumbbell_axis{axis_idx+1}.png",
                       axis_idx=axis_idx)
        # only generate loadings png for axis 2/3 (axis 1 already exists)
        if axis_idx in (1, 2):
            plot_loadings(cca_out["h_loadings"], cca_out["p_loadings"],
                          cell_cols, gene_index, axis_idx, slide_key,
                          out_dir / f"cca_loadings_axis{axis_idx+1}.png")
    plot_section_means(common, sections, Hc, Pc, train_rs, slide_key, scfg,
                        out_dir / "cca_section_means.png")

    # console: per-section means + top loaders for each axis (helps decide labels)
    for axis_idx in range(3):
        print(f"  --- axis {axis_idx+1}  (r = {train_rs[axis_idx]:+.3f}) ---")
        print(f"  section means (H2C / Pro):")
        for s in scfg["order"]:
            mask = np.array(sections) == s
            if mask.sum() == 0: continue
            print(f"    {s} ({scfg['labels'][s]:<26}): "
                  f"H2C {Hc[mask, axis_idx].mean():+.2f}  /  Pro {Pc[mask, axis_idx].mean():+.2f}  "
                  f"(n={mask.sum()})")
        print_top_loaders(cca_out["h_loadings"], cca_out["p_loadings"],
                          cell_cols, gene_index, axis_idx)

    print(f"  --- per-axis intra/inter check (1000-perm null) ---")
    intra_rows = []
    for axis_idx in range(3):
        for mod_name, scores in [("Hist2Cell", Hc[:, axis_idx]),
                                  ("Proteomics", Pc[:, axis_idx])]:
            stats = axis_intra_inter(scores, sections, seed=42)
            row = {"slide": slide_key, "axis": axis_idx+1, "modality": mod_name, **stats}
            intra_rows.append(row)
            print(f"    axis {axis_idx+1}  {mod_name:<10}  "
                  f"intra={stats['intra_mean']:.2f}  inter={stats['inter_mean']:.2f}  "
                  f"ratio={stats['ratio']:.3f}  "
                  f"null 95%=[{stats['null_ratio_95lo']:.3f}, {stats['null_ratio_95hi']:.3f}]  "
                  f"p={stats['p_perm']:.4f}")
    pd.DataFrame(intra_rows).to_csv(out_dir / "cca_axis_intra_inter.csv", index=False)
    plot_axis_intra_inter(intra_rows, slide_key,
                           out_dir / "cca_axis_intra_inter.png")

    print(f"  saved → cca_dumbbell_axis{{1,2,3}}.png, "
          f"cca_loadings_axis{{2,3}}.png, cca_section_means.png, "
          f"cca_scores_per_roi.csv, cca_axis_intra_inter.{{csv,png}}")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
