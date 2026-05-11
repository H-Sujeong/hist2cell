"""Two-claim focused proof for slide1 (per user feedback "분석이 너무 과하다").

Claim 1 — Tumor-internal cross-modality risk-correlation is positive
  Within the Tumor compartment (sections a + b + t = 32 ROIs) the
  proteomics smooth-muscle / fibroblast panel mean (log2 intensity)
  correlates with the matching Hist2Cell lineage-group sum across ROIs
  (Pearson r ≈ +0.3 to +0.4).  The Tumor a-vs-b Wilcoxon already showed
  both modalities flag the same direction (MYH11 / TAGLN ↑ in a,
  Stromal-muscle ↑ in a; broad-proxy ↑ in a, COLGALT1 / EPPK1 / ALDH1B1
  ↑ in a).

  Limitation made explicit
  ------------------------
  When all 46 ROIs are pooled (Tumor + T-cell + control) the
  correlation flattens / inverts because the lung-trained Hist2Cell
  top-5 per ROI is uniformly epithelial-dominant — see Claim 2.  This
  is a known consequence of lung→breast cross-tissue proxying
  (EPITHELIAL_PROXY_METHODOLOGY.md) and is the same reason the
  per-cell-type top-5 list barely varies between Tumor and T-cell
  sections.

Claim 2 — Per-ROI top-expressed cell types
  For each of 47 ROIs we list the top-5 Hist2Cell cell types.  We also
  collapse to lineage groups and show that every section's top-5
  draws ~95% from {Epithelial, Stromal, Vascular} — i.e. lung
  Hist2Cell does NOT place immune types in the top-5 even for T-cell
  ROIs.  The table is informative about which lung labels the model
  emphasises in each ROI; the cross-tissue caveat means the labels
  should be read as morphology categories, not cell-type ground truth.

Inputs (relative to this folder)
  ../cell_typing/roi_signatures.csv     47 × (80 + 3 scores) — ROI sigs
  ../../report.gg_matrix (1).tsv        proteomics log2 source
  ../../../analysis/cell_type_groups.csv  groups + proxy flags

Outputs (this folder)
  cross_modality_correlations.csv   panel × subset × {Pearson, Spearman}
  cross_modality_scatter.png        per-panel scatter + r annotated
                                    (subset = Tumor compartment a+b+t)
  roi_top_celltypes.csv             47 × {tube_id, section, n_spots,
                                          top1..top5}
  section_group_composition.csv     per-section top-5 group share (%)
  section_group_composition.png     stacked bar of the same
  roi_top_celltypes_heatmap.png     ROI × union-of-top cell types
                                    (z-score per cell type)
  summary.md                        2-claim focused write-up
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
    "a": "High-risk Tumor", "b": "Low-risk Tumor",
    "c": "High-risk T-cell", "d": "Low-risk T-cell",
    "t": "Middle-risk Tumor (ctrl)",
}
SECTION_COLOR = {"a":"#d62728","b":"#1f77b4","c":"#2ca02c","d":"#9467bd","t":"#7f7f7f"}
SECTION_ORDER = ["a", "b", "c", "d", "t"]

# Hist2Cell types matching each proteomics marker panel
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
    ("all (a+b+c+d+t)",  ["a", "b", "c", "d", "t"]),
    ("Tumor (a+b+t)",    ["a", "b", "t"]),
    ("Tumor a vs b only",["a", "b"]),
    ("T-cell (c+d)",     ["c", "d"]),
]


# ---- helpers ----

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
    slide1_cols = [c for c in raw.columns[3:] if c[0] in "abcdt"]
    raw[slide1_cols] = raw[slide1_cols].replace(0, np.nan)
    log2 = np.log2(raw[slide1_cols])
    log2.index = raw["Genes"].astype(str)
    return sig, groups, log2, slide1_cols


# ---- Claim 1 ----

def cross_modality_correlations(sig, groups, log2, slide1_cols):
    cell_cols = groups["cell_type"].tolist()
    common = [t for t in sig["tube_id"] if t in slide1_cols]
    sig_by_tube = sig.set_index("tube_id")
    rows = []
    paired_data = {}     # panel → DataFrame of (prot, h2c, section, tube)
    for panel_name, markers, spec in PANELS:
        h2c_types = types_for_panel(spec, groups)
        present = [g for g in markers if g in log2.index]
        if not present or not h2c_types:
            continue
        prot = log2.loc[present, common].mean(axis=0, skipna=True)
        h2c  = sig_by_tube.loc[common, h2c_types].sum(axis=1)
        full = pd.DataFrame({
            "tube": common, "section": [t[0] for t in common],
            "prot": prot.values, "h2c": h2c.values,
        }).dropna()
        paired_data[panel_name] = full
        for subset_label, sec_set in SUBSETS:
            sub = full[full["section"].isin(sec_set)]
            if len(sub) < 4:
                continue
            rp, pp = pearsonr(sub["prot"], sub["h2c"])
            rs, ps = spearmanr(sub["prot"], sub["h2c"])
            rows.append({
                "panel": panel_name,
                "subset": subset_label,
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
    """One panel per measure, scatter restricted to Tumor (a+b+t) — the
    subset where the positive correlation is cleanest."""
    panels = list(paired_data.items())
    n = len(panels); ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows))
    axes = list(axes.flat)
    for ax, (panel, df) in zip(axes, panels):
        sub = df[df["section"].isin(["a", "b", "t"])]
        if len(sub) < 4:
            ax.axis("off"); continue
        rp, pp = pearsonr(sub["prot"], sub["h2c"])
        for s in ["a", "b", "t"]:
            ss = sub[sub.section == s]
            if len(ss):
                ax.scatter(ss["prot"], ss["h2c"], s=55,
                           c=SECTION_COLOR[s], edgecolor="black",
                           linewidth=0.4, label=SECTION_LABEL[s], alpha=0.85)
        xs = np.linspace(sub["prot"].min(), sub["prot"].max(), 50)
        slope, intercept = np.polyfit(sub["prot"], sub["h2c"], 1)
        ax.plot(xs, intercept + slope*xs, c="black", linewidth=0.8, alpha=0.5)
        ax.set_title(f"{panel}\nTumor (a+b+t) r={rp:+.3f}  p={pp:.2e}  n={len(sub)}",
                     fontsize=10)
        ax.set_xlabel("proteomics mean log2 (markers)", fontsize=8)
        ax.set_ylabel("Hist2Cell matching-types sum", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes[n:]:
        ax.axis("off")
    handles = [plt.Line2D([0],[0], marker="o", color="w",
                           markerfacecolor=SECTION_COLOR[s], markersize=8,
                           label=SECTION_LABEL[s]) for s in ["a","b","t"]]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.02), fontsize=9)
    fig.suptitle("Claim 1 — Tumor-internal cross-modality correlation\n"
                 "(proteomics marker mean vs Hist2Cell matching-types sum, "
                 "Tumor compartment only)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---- Claim 2 ----

def per_roi_top_celltypes(sig, groups, top_n=5):
    cell_cols = groups["cell_type"].tolist()
    g_map = dict(zip(groups["cell_type"], groups["group"]))
    rows = []
    for _, r in sig.iterrows():
        order = r[cell_cols].sort_values(ascending=False)
        top = order.head(top_n)
        rec = {"tube_id": r.tube_id, "section": r.section,
               "section_label": SECTION_LABEL[r.section],
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
        if len(sub) == 0:
            continue
        counts = {g: 0 for g in all_groups}
        total = 0
        for _, r in sub.iterrows():
            top = r[cell_cols].sort_values(ascending=False).head(top_n)
            for ct in top.index:
                counts[g_map[ct]] += 1
                total += 1
        rec = {"section": sec, "section_label": SECTION_LABEL[sec],
               "n_ROIs": len(sub), "total_top_slots": total}
        for g in all_groups:
            rec[f"pct_{g}"] = round(100.0 * counts[g] / total, 1) if total else 0
        rows.append(rec)
    return pd.DataFrame(rows), all_groups


def plot_section_composition(comp_df, all_groups, out_path):
    fig, ax = plt.subplots(figsize=(11, 5))
    # stacked bar per section
    bottom = np.zeros(len(comp_df))
    palette = plt.cm.tab10.colors
    for i, g in enumerate(all_groups):
        col = f"pct_{g}"
        vals = comp_df[col].values
        ax.bar(np.arange(len(comp_df)), vals, bottom=bottom,
               label=g, color=palette[i % len(palette)],
               edgecolor="white", linewidth=0.5)
        # annotate non-zero
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
    ax.set_title("Claim 2 — per-section composition of Hist2Cell top-5 cell types\n"
                 "(every section is Epithelial/Stromal/Vascular-dominant; "
                 "Immune ≈ 0% — a lung-proxy limitation, see findings)",
                 fontsize=11)
    ax.set_ylim(0, 100)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8,
              frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_roi_top_heatmap(sig, groups, top_df, out_path):
    cell_cols = groups["cell_type"].tolist()
    union_top = sorted({top_df[c].iloc[i].split(" (")[0]
                        for c in ["top1","top2","top3","top4","top5"]
                        for i in range(len(top_df))})
    M = sig.set_index("tube_id").loc[top_df["tube_id"], union_top]
    means = M.mean(axis=0)
    stds  = M.std(axis=0).replace(0, 1)
    Z = (M - means) / stds
    fig, ax = plt.subplots(figsize=(0.35 * len(union_top) + 4,
                                    0.30 * len(top_df) + 2))
    im = ax.imshow(Z.values, aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5)
    ax.set_yticks(np.arange(len(top_df)))
    ax.set_yticklabels(top_df["tube_id"], fontsize=7)
    ax.set_xticks(np.arange(len(union_top)))
    ax.set_xticklabels(union_top, rotation=90, fontsize=7)
    for i, s in enumerate(top_df["section"]):
        ax.add_patch(plt.Rectangle((-1.5, i-0.5), 1, 1,
                                    color=SECTION_COLOR[s], clip_on=False))
    ax.set_xlim(-1.6, len(union_top) - 0.5)
    ax.set_title("Per-ROI top cell-type heatmap "
                 "(union of every ROI's top-5, z-score across ROIs)\n"
                 "left strip = section colour",
                 fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.02, label="z (across ROIs)")
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                           markerfacecolor=SECTION_COLOR[s], markersize=10,
                           label=SECTION_LABEL[s]) for s in SECTION_ORDER]
    ax.legend(handles=handles, loc="upper right",
              bbox_to_anchor=(1.25, 1.0), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---- main ----

def main():
    print("[load] roi_signatures + gg_matrix + cell_type_groups")
    sig, groups, log2, slide1_cols = load_inputs()
    common = [t for t in sig["tube_id"] if t in slide1_cols]
    print(f"  tubes={len(sig)}, proteomics samples={len(slide1_cols)}, "
          f"intersection={len(common)}")

    # ── Claim 1 ──
    print("\n[Claim 1] cross-modality correlation per panel × subset")
    corr_df, paired = cross_modality_correlations(sig, groups, log2, slide1_cols)
    corr_df.to_csv(HERE / "cross_modality_correlations.csv", index=False)
    pivot = corr_df.pivot_table(index="panel", columns="subset",
                                 values="pearson_r")
    print(pivot[["all (a+b+c+d+t)", "Tumor (a+b+t)",
                 "Tumor a vs b only", "T-cell (c+d)"]].to_string(
                    float_format=lambda x: f"{x:+.3f}" if not pd.isna(x) else "  na"))
    plot_cross_modality_tumor_subset(paired,
                                     HERE / "cross_modality_scatter.png")

    # ── Claim 2 ──
    print("\n[Claim 2] per-ROI top-5 cell types")
    top_df = per_roi_top_celltypes(sig, groups)
    top_df.to_csv(HERE / "roi_top_celltypes.csv", index=False)
    print(top_df[["tube_id", "section_label", "top1", "top2"]].head(10).to_string(
        index=False))
    print(f"\n  per-section group composition of top-5 slots:")
    comp_df, all_groups = section_group_composition(sig, groups)
    comp_df.to_csv(HERE / "section_group_composition.csv", index=False)
    print(comp_df[["section_label", "n_ROIs"] +
                  [f"pct_{g}" for g in all_groups]].to_string(index=False))
    plot_section_composition(comp_df, all_groups,
                              HERE / "section_group_composition.png")
    plot_roi_top_heatmap(sig, groups, top_df,
                          HERE / "roi_top_celltypes_heatmap.png")

    print("\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png", ".md"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
