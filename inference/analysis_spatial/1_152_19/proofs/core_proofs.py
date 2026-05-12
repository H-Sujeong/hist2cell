"""slide2 (1_152_19) — focused 2-claim analysis (mirrors slide1/proofs/).

Claim 1 — Tumor risk-axis direction agreement (correlation magnitude
          NOT claimed — same convention as slide1)
Claim 2 — per-ROI top cell types + per-section group composition

Inputs (relative to this folder)
  ../cell_typing/roi_signatures.csv         48 × (80 + 3 scores)
  ../cell_typing/marker_hypotheses.csv      Hist2Cell-side direction check
  ../proteomics/marker_hypothesis_check.csv proteomics-side direction check
  ../../report.gg_matrix (1).tsv            log2 source
  ../../../analysis/cell_type_groups.csv

Outputs
  cross_modality_correlations.csv  exploratory panel × subset r table
  cross_modality_scatter.png       Tumor (e+f+v) per-panel scatter
  roi_top_celltypes.csv            48 × top1..top5 + lineage group
  roi_top_celltypes_heatmap.png    ROI × union-of-top z-score heatmap
  section_group_composition.csv    5 section × 10 lineage % share
  section_group_composition.png    stacked-bar
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


HERE  = Path(__file__).resolve().parent
ROI   = HERE.parent / "cell_typing" / "roi_signatures.csv"
GG    = Path("/home/sjhong/hist2cell/inference/analysis_spatial/report.gg_matrix (1).tsv")
GRP   = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")

SECTION_LABEL = {
    "e": "High-risk Tumor", "f": "Low-risk Tumor",
    "g": "High-risk T-cell", "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)", "w": "Middle-risk T-cell (ctrl)",
}
SECTION_COLOR = {"e":"#d62728","f":"#1f77b4","g":"#2ca02c","h":"#9467bd","v":"#7f7f7f","w":"#bcbd22"}
SECTION_ORDER = ["e", "f", "g", "h", "v"]

PANELS = [
    ("Smooth muscle",
     ["MYH11", "TAGLN", "CNN1", "MYLK"],
     {"group": "Stromal-muscle"}),
    ("Fibroblast",
     ["COL1A1", "COL3A1", "DCN", "LUM", "FAP", "VIM"],
     {"group": "Stromal-fibroblast"}),
    ("Epithelial (broad proxy)",
     ["KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "EPCAM", "CDH1"],
     {"flag": "is_broad_proxy"}),
    ("Macrophage",
     ["CD163", "LYZ"],
     {"name_prefix": ["Macro_", "Macrophage_"]}),
    ("B cell",
     ["IGHM", "IGHG1", "JCHAIN"],
     {"name_prefix": ["B_"]}),
    ("Endothelial",
     ["PECAM1", "VWF"],
     {"group": "Vascular"}),
]

SUBSETS = [
    ("all (e+f+g+h+v)",    SECTION_ORDER),
    ("Tumor (e+f+v)",      ["e", "f", "v"]),
    ("Tumor e vs f only",  ["e", "f"]),
    ("T-cell (g+h)",       ["g", "h"]),
]


def types_for_panel(spec, groups):
    if "flag" in spec:
        return groups[groups[spec["flag"]] == 1]["cell_type"].tolist()
    if "group" in spec:
        return groups[groups["group"] == spec["group"]]["cell_type"].tolist()
    if "name_prefix" in spec:
        return [c for c in groups["cell_type"]
                if any(c.startswith(p) for p in spec["name_prefix"])]
    return []


def load_inputs():
    sig = pd.read_csv(ROI)
    groups = pd.read_csv(GRP)
    raw = pd.read_csv(GG, sep="\t")
    slide2_cols = [c for c in raw.columns[3:] if c[0] in "efghv"]
    raw[slide2_cols] = raw[slide2_cols].replace(0, np.nan)
    log2 = np.log2(raw[slide2_cols])
    log2.index = raw["Genes"].astype(str)
    return sig, groups, log2, slide2_cols


def cross_modality_correlations(sig, groups, log2, slide2_cols):
    common = [t for t in sig["tube_id"] if t in slide2_cols]
    sig_by_tube = sig.set_index("tube_id")
    rows, paired_data = [], {}
    for panel_name, markers, spec in PANELS:
        h2c_types = types_for_panel(spec, groups)
        present = [g for g in markers if g in log2.index]
        if not present or not h2c_types:
            continue
        prot = log2.loc[present, common].mean(axis=0, skipna=True)
        h2c  = sig_by_tube.loc[common, h2c_types].sum(axis=1)
        full = pd.DataFrame({"tube": common,
                              "section": [t[0] for t in common],
                              "prot": prot.values, "h2c": h2c.values}).dropna()
        paired_data[panel_name] = full
        for subset_label, sec_set in SUBSETS:
            sub = full[full["section"].isin(sec_set)]
            if len(sub) < 4:
                continue
            rp, pp = pearsonr(sub["prot"], sub["h2c"])
            rs, ps = spearmanr(sub["prot"], sub["h2c"])
            rows.append({
                "panel": panel_name, "subset": subset_label,
                "n_ROI": len(sub),
                "n_markers": len(present),
                "n_h2c_types": len(h2c_types),
                "pearson_r": round(float(rp), 3),
                "pearson_p": float(pp),
                "spearman_r": round(float(rs), 3),
                "spearman_p": float(ps),
            })
    return pd.DataFrame(rows), paired_data


def plot_cross_modality_tumor_subset(paired_data, out_path):
    panels = list(paired_data.items())
    n = len(panels); ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows))
    axes = list(axes.flat)
    for ax, (panel, df) in zip(axes, panels):
        sub = df[df["section"].isin(["e", "f", "v"])]
        if len(sub) < 4:
            ax.axis("off"); continue
        rp, pp = pearsonr(sub["prot"], sub["h2c"])
        for s in ["e", "f", "v"]:
            ss = sub[sub.section == s]
            if len(ss):
                ax.scatter(ss["prot"], ss["h2c"], s=55,
                           c=SECTION_COLOR[s], edgecolor="black",
                           linewidth=0.4, label=SECTION_LABEL[s], alpha=0.85)
        xs = np.linspace(sub["prot"].min(), sub["prot"].max(), 50)
        slope, intercept = np.polyfit(sub["prot"], sub["h2c"], 1)
        ax.plot(xs, intercept + slope*xs, c="black", linewidth=0.8, alpha=0.5)
        ax.set_title(f"{panel}\nTumor (e+f+v) r={rp:+.3f}  p={pp:.2e}  n={len(sub)}",
                     fontsize=10)
        ax.set_xlabel("proteomics mean log2 (markers)", fontsize=8)
        ax.set_ylabel("Hist2Cell matching-types sum", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes[n:]:
        ax.axis("off")
    handles = [plt.Line2D([0],[0], marker="o", color="w",
                           markerfacecolor=SECTION_COLOR[s], markersize=8,
                           label=SECTION_LABEL[s]) for s in ["e","f","v"]]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), fontsize=9)
    fig.suptitle("Exploratory cross-modality correlation per panel "
                 "(Tumor compartment only) — slide2", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def per_roi_top_celltypes(sig, groups, top_n=5):
    cell_cols = groups["cell_type"].tolist()
    g_map = dict(zip(groups["cell_type"], groups["group"]))
    rows = []
    for _, r in sig.iterrows():
        order = r[cell_cols].sort_values(ascending=False)
        top = order.head(top_n)
        rec = {"tube_id": r.tube_id, "section": r.section,
               "section_label": SECTION_LABEL.get(r.section, "?"),
               "n_spots": int(r.n_spots)}
        for i, (ct, v) in enumerate(top.items()):
            rec[f"top{i+1}"] = f"{ct} ({v:.2f})"
            rec[f"top{i+1}_group"] = g_map.get(ct, "?")
        rows.append(rec)
    return pd.DataFrame(rows)


def section_group_composition(sig, groups, top_n=5):
    cell_cols = groups["cell_type"].tolist()
    g_map = dict(zip(groups["cell_type"], groups["group"]))
    all_groups = sorted(set(groups["group"]))
    rows = []
    for sec in SECTION_ORDER:
        sub = sig[sig["section"] == sec]
        if len(sub) == 0: continue
        counts = {g: 0 for g in all_groups}
        total = 0
        for _, r in sub.iterrows():
            top = r[cell_cols].sort_values(ascending=False).head(top_n)
            for ct in top.index:
                counts[g_map[ct]] += 1; total += 1
        rec = {"section": sec, "section_label": SECTION_LABEL[sec],
               "n_ROIs": len(sub), "total_top_slots": total}
        for g in all_groups:
            rec[f"pct_{g}"] = round(100.0 * counts[g] / total, 1) if total else 0
        rows.append(rec)
    return pd.DataFrame(rows), all_groups


def plot_section_composition(comp_df, all_groups, out_path):
    fig, ax = plt.subplots(figsize=(11, 5))
    bottom = np.zeros(len(comp_df))
    palette = plt.cm.tab10.colors
    for i, g in enumerate(all_groups):
        vals = comp_df[f"pct_{g}"].values
        ax.bar(np.arange(len(comp_df)), vals, bottom=bottom, label=g,
               color=palette[i % len(palette)], edgecolor="white", linewidth=0.5)
        for j, v in enumerate(vals):
            if v >= 5:
                ax.text(j, bottom[j] + v/2, f"{v:.0f}",
                        ha="center", va="center", fontsize=8, color="white",
                        fontweight="bold")
        bottom += vals
    ax.set_xticks(np.arange(len(comp_df)))
    ax.set_xticklabels([SECTION_LABEL[s] for s in comp_df["section"]],
                       rotation=20, fontsize=9, ha="right")
    ax.set_ylabel("share of per-ROI top-5 slots (%)")
    ax.set_title("slide2 — per-section composition of Hist2Cell top-5 cell types",
                 fontsize=11)
    ax.set_ylim(0, 100)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8,
              frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_roi_top_heatmap(sig, groups, top_df, out_path):
    union_top = sorted({top_df[c].iloc[i].split(" (")[0]
                        for c in ["top1","top2","top3","top4","top5"]
                        for i in range(len(top_df))})
    M = sig.set_index("tube_id").loc[top_df["tube_id"], union_top]
    means = M.mean(axis=0); stds = M.std(axis=0).replace(0, 1)
    Z = (M - means) / stds

    from matplotlib.gridspec import GridSpec
    fig_w = 0.32 * len(union_top) + 5.5
    fig_h = 0.30 * len(top_df) + 2
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = GridSpec(1, 3, width_ratios=[0.4, len(union_top), 0.3],
                   wspace=0.05, figure=fig)
    ax_strip = fig.add_subplot(gs[0, 0])
    ax_heat  = fig.add_subplot(gs[0, 1])
    ax_cbar  = fig.add_subplot(gs[0, 2])

    im = ax_heat.imshow(Z.values, aspect="auto", cmap="RdBu_r",
                        vmin=-2.5, vmax=2.5)
    ax_heat.set_yticks([])
    ax_heat.set_xticks(np.arange(len(union_top)))
    ax_heat.set_xticklabels(union_top, rotation=90, fontsize=7)

    for i, s in enumerate(top_df["section"]):
        ax_strip.add_patch(plt.Rectangle((0, i-0.5), 1, 1,
                                          color=SECTION_COLOR[s]))
    ax_strip.set_xlim(0, 1)
    ax_strip.set_ylim(len(top_df) - 0.5, -0.5)
    ax_strip.set_xticks([])
    ax_strip.set_yticks(np.arange(len(top_df)))
    ax_strip.set_yticklabels(top_df["tube_id"], fontsize=7)
    ax_strip.set_ylabel("tube_id (section colour)", fontsize=8)
    for sp in ("top", "right", "bottom"):
        ax_strip.spines[sp].set_visible(False)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("z (across ROIs)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    handles = [plt.Line2D([0],[0], marker="s", color="w",
                           markerfacecolor=SECTION_COLOR[s], markersize=10,
                           label=SECTION_LABEL[s]) for s in SECTION_ORDER]
    fig.legend(handles=handles, loc="upper left",
               bbox_to_anchor=(0.84, 0.97), fontsize=8, frameon=True,
               framealpha=0.9)

    fig.suptitle("slide2 — Per-ROI top cell-type heatmap "
                 "(union of every ROI's top-5, z-score across ROIs)",
                 fontsize=11, y=0.995)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    print("[load] slide2 roi_signatures + gg_matrix + groups")
    sig, groups, log2, slide2_cols = load_inputs()
    common = [t for t in sig["tube_id"] if t in slide2_cols]
    print(f"  tubes={len(sig)}, proteomics samples={len(slide2_cols)}, "
          f"intersection={len(common)}")

    print("\n[Claim 1] cross-modality correlation per panel × subset (exploratory)")
    corr_df, paired = cross_modality_correlations(sig, groups, log2, slide2_cols)
    corr_df.to_csv(HERE / "cross_modality_correlations.csv", index=False)
    pivot = corr_df.pivot_table(index="panel", columns="subset",
                                 values="pearson_r")
    cols = [c for c in ["all (e+f+g+h+v)", "Tumor (e+f+v)",
                        "Tumor e vs f only", "T-cell (g+h)"]
            if c in pivot.columns]
    print(pivot[cols].to_string(
        float_format=lambda x: f"{x:+.3f}" if not pd.isna(x) else "  na"))
    plot_cross_modality_tumor_subset(paired,
                                     HERE / "cross_modality_scatter.png")

    print("\n[Claim 2] per-ROI top cell types + per-section composition")
    top_df = per_roi_top_celltypes(sig, groups)
    top_df.to_csv(HERE / "roi_top_celltypes.csv", index=False)
    print(top_df[["tube_id", "section_label", "top1", "top2"]].head(10).to_string(
        index=False))
    comp_df, all_groups = section_group_composition(sig, groups)
    comp_df.to_csv(HERE / "section_group_composition.csv", index=False)
    print("\nper-section group composition (% of top-5 slots):")
    print(comp_df[["section_label", "n_ROIs"] +
                  [f"pct_{g}" for g in all_groups]].to_string(index=False))
    plot_section_composition(comp_df, all_groups,
                             HERE / "section_group_composition.png")
    plot_roi_top_heatmap(sig, groups, top_df,
                         HERE / "roi_top_celltypes_heatmap.png")

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png", ".md"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
