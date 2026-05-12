"""slide2 (1_152_19) — proof_ver2: data-driven cross-modality validation.

Deliberately re-derives every cross-modality claim from the slide's own
Hist2Cell × proteomics matrices, *ignoring* the collaborator's
pre-selected marker panel.  See ../../_proof_ver2_lib.py for the
algorithm (PCA → CCA + permutation null; all-pair Pearson with BH-FDR;
per-ROI cosine similarity built from discovered markers).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))   # analysis_spatial/ on sys.path

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, run_cca, permutation_null,
    discover_marker_pairs, per_roi_cosine,
)


SECTION_LABEL = {
    "e": "High-risk Tumor", "f": "Low-risk Tumor",
    "g": "High-risk T-cell", "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)",
}
SECTION_COLOR = {"e":"#d62728","f":"#1f77b4","g":"#2ca02c","h":"#9467bd","v":"#7f7f7f"}
SECTION_ORDER = ["e", "f", "g", "h", "v"]

CFG = SlideConfig(
    name="slide2 (1_152_19)",
    pred_csv=Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv"),
    roi_pkl=HERE.parent / "1_152_19_ROI_groups.pkl",
    npy=HERE.parent / "meteo_1_152_19_coords.npy",
    section_label=SECTION_LABEL,
    section_color=SECTION_COLOR,
    sample_section_prefixes="efghv",
    out_dir=HERE,
)


def plot_cca_scatter(Hc, Pc, sections, train_rs, out_path):
    n_comp = Hc.shape[1]
    fig, axes = plt.subplots(1, n_comp, figsize=(5.5*n_comp, 5))
    if n_comp == 1: axes = [axes]
    for i, ax in enumerate(axes):
        for s in SECTION_ORDER:
            mask = [j for j, sec in enumerate(sections) if sec == s]
            if mask:
                ax.scatter(Hc[mask, i], Pc[mask, i], s=80,
                           c=SECTION_COLOR[s], edgecolor="black",
                           linewidth=0.4, label=SECTION_LABEL[s], alpha=0.85)
        xs = np.linspace(Hc[:, i].min(), Hc[:, i].max(), 50)
        slope, intercept = np.polyfit(Hc[:, i], Pc[:, i], 1)
        ax.plot(xs, slope*xs + intercept, c="black", linewidth=0.8, alpha=0.5)
        ax.set_xlabel(f"Hist2Cell canonical {i+1}", fontsize=9)
        ax.set_ylabel(f"Proteomics canonical {i+1}", fontsize=9)
        ax.set_title(f"Canonical pair {i+1}: r = {train_rs[i]:+.3f}", fontsize=11)
        if i == 0:
            ax.legend(loc="best", fontsize=8)
    fig.suptitle("CCA — slide2: data-driven canonical correlations between modalities",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_permutation_null(null_top, observed, out_path):
    null_top = null_top[~np.isnan(null_top)]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(null_top, bins=40, color="#bbbbbb", edgecolor="white", alpha=0.85,
            label=f"permutation null (n={len(null_top)})")
    ax.axvline(observed, color="#d62728", linewidth=2,
               label=f"observed top canonical r = {observed:+.3f}")
    p_emp = float(np.mean(np.abs(null_top) >= np.abs(observed)))
    ax.set_xlabel("top canonical correlation r")
    ax.set_ylabel("permutation count")
    ax.set_title(f"slide2 — permutation null for top canonical r  "
                 f"(empirical p = {p_emp:.4f})", fontsize=11)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return p_emp


def plot_top_loaders(h_load, p_load, cell_cols, gene_index, axis_idx,
                     out_path, top_n=12):
    h_axis = pd.Series(h_load[:, axis_idx], index=cell_cols)
    p_axis = pd.Series(p_load[:, axis_idx], index=gene_index)
    h_top_pos = h_axis.nlargest(top_n)
    h_top_neg = h_axis.nsmallest(top_n)
    p_top_pos = p_axis.nlargest(top_n)
    p_top_neg = p_axis.nsmallest(top_n)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for ax, (data, title, color) in zip(
            axes.flat,
            [(h_top_pos, f"Hist2Cell top + loaders (axis {axis_idx+1})", "#d62728"),
             (h_top_neg, f"Hist2Cell top − loaders (axis {axis_idx+1})", "#1f77b4"),
             (p_top_pos, f"Proteomics top + loaders (axis {axis_idx+1})", "#d62728"),
             (p_top_neg, f"Proteomics top − loaders (axis {axis_idx+1})", "#1f77b4")]):
        data = data.sort_values()
        ax.barh(range(len(data)), data.values, color=color, alpha=0.75,
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(data)))
        ax.set_yticklabels(data.index, fontsize=8)
        ax.set_xlabel("loading on canonical axis")
        ax.set_title(title, fontsize=10)
        ax.axvline(0, color="black", linewidth=0.4)
    fig.suptitle(f"CCA loadings — slide2 canonical axis {axis_idx+1}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_per_roi_similarity(sims, sections, out_path):
    fig, ax = plt.subplots(figsize=(14, 5))
    pairs = sorted(enumerate(zip(sections, sims)),
                   key=lambda kv: (SECTION_ORDER.index(kv[1][0]), kv[0]))
    ordered_idx = [k for k, _ in pairs]
    ordered_sec = [sections[k] for k in ordered_idx]
    ordered_sim = [sims[k] for k in ordered_idx]
    colors = [SECTION_COLOR[s] for s in ordered_sec]
    ax.bar(range(len(ordered_sim)), ordered_sim, color=colors,
           edgecolor="black", linewidth=0.3)
    ax.axhline(np.nanmean(sims), color="black", linestyle="--",
               linewidth=0.8, label=f"mean = {np.nanmean(sims):+.3f}")
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_xticks(range(len(ordered_sec)))
    ax.set_xticklabels([f"{s}" for s in ordered_sec], rotation=0, fontsize=6)
    ax.set_ylabel("cosine similarity")
    ax.set_title("slide2 — per-ROI cross-modality similarity "
                 "(Hist2Cell vector ↔ proteomics-derived cell-type score)",
                 fontsize=11)
    handles = [plt.Rectangle((0,0),1,1, color=SECTION_COLOR[s],
                              label=SECTION_LABEL[s]) for s in SECTION_ORDER]
    ax.legend(handles=handles + ax.get_legend_handles_labels()[0], loc="best",
              fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    print(f"[{CFG.name}] build Hist2Cell ROI signatures + load proteomics")
    sig_df, cell_cols = build_roi_signatures(CFG)
    log2_f, slide_cols = load_proteomics_matrix(CFG)
    print(f"  Hist2Cell tubes: {len(sig_df)}  proteomics samples: {len(slide_cols)}")
    print(f"  proteomics genes after detect≥50% filter: {len(log2_f)}")

    common, H, P, sig_aligned, gene_index = align_modalities(
        sig_df, log2_f, slide_cols, cell_cols)
    sections = list(sig_aligned["section"])
    print(f"  ROIs with both modalities: {len(common)}  "
          f"(H matrix {H.shape}, P matrix {P.shape})")

    print("\n[CCA] PCA → CCA → train canonical correlations")
    cca_out = run_cca(H, P)
    print(f"  PCA explained var (Hist2Cell first 3): "
          f"{', '.join(f'{v*100:.1f}%' for v in cca_out['pca_h_var'][:3])}")
    print(f"  PCA explained var (Proteomics first 3): "
          f"{', '.join(f'{v*100:.1f}%' for v in cca_out['pca_p_var'][:3])}")
    print(f"  Top {len(cca_out['train_rs'])} canonical correlations: "
          f"{', '.join(f'{r:+.3f}' for r in cca_out['train_rs'])}")

    print("\n[Perm] permutation null (top canonical r)")
    null_top = permutation_null(H, P)
    p_emp = float(np.mean(np.abs(null_top[~np.isnan(null_top)])
                          >= np.abs(cca_out["train_rs"][0])))
    print(f"  observed top r = {cca_out['train_rs'][0]:+.3f}, "
          f"empirical p (two-sided permutation) = {p_emp:.4f}")
    print(f"  null distribution mean = {np.nanmean(null_top):+.3f}, "
          f"95% range = [{np.nanpercentile(null_top, 2.5):+.3f}, "
          f"{np.nanpercentile(null_top, 97.5):+.3f}]")

    cca_summary = pd.DataFrame({
        "axis": range(1, len(cca_out["train_rs"])+1),
        "train_r": cca_out["train_rs"],
    })
    cca_summary["observed_top_r"] = cca_out["train_rs"][0]
    cca_summary["null_mean"] = np.nanmean(null_top)
    cca_summary["null_95_lo"] = np.nanpercentile(null_top, 2.5)
    cca_summary["null_95_hi"] = np.nanpercentile(null_top, 97.5)
    cca_summary["permutation_p_top_axis"] = p_emp
    cca_summary.to_csv(HERE/"cca_summary.csv", index=False)

    print("\n[Discover] all-pair Pearson + BH-FDR (top 5/type each direction)")
    pairs = discover_marker_pairs(H, P, cell_cols, gene_index)
    pairs.to_csv(HERE/"discovered_marker_pairs.csv", index=False)
    sig_pos = pairs[(pairs.direction=="pos") & (pairs.p_bh < 0.05)]
    sig_neg = pairs[(pairs.direction=="neg") & (pairs.p_bh < 0.05)]
    print(f"  BH<.05 positive pairs: {len(sig_pos)}")
    print(f"  BH<.05 negative pairs: {len(sig_neg)}")
    print(f"  top 10 positive (highest r):")
    print(sig_pos.nlargest(10, "r")[["cell_type","gene","r","p","p_bh"]]
          .to_string(index=False))

    print("\n[Cosine] per-ROI cross-modality similarity")
    P_score, used_genes, score_cell_order = per_roi_cosine(
        H, P, gene_index, pairs, top_n_per_type=3)
    marker_cell_cols = [c for c in score_cell_order
                        if used_genes.get(c)]
    cols_idx = [score_cell_order.index(c) for c in marker_cell_cols]
    P_score_sub = P_score[:, cols_idx]
    H_idx = [cell_cols.index(c) for c in marker_cell_cols]
    H_sub = H[:, H_idx]
    sims = []
    for i in range(H_sub.shape[0]):
        h, p = H_sub[i], P_score_sub[i]
        if np.linalg.norm(h) == 0 or np.linalg.norm(p) == 0:
            sims.append(np.nan)
        else:
            sims.append(float(np.dot(h, p) / (np.linalg.norm(h)*np.linalg.norm(p))))
    sims = np.array(sims)
    print(f"  per-ROI cosine similarity: mean = {np.nanmean(sims):+.3f}, "
          f"range = [{np.nanmin(sims):+.3f}, {np.nanmax(sims):+.3f}]")
    pd.DataFrame({"tube_id": common, "section": sections,
                  "cosine_similarity": sims}).to_csv(
        HERE/"per_roi_cosine_similarity.csv", index=False)

    print("\n[plot] saving 5 PNGs")
    plot_cca_scatter(cca_out["Hc"], cca_out["Pc"], sections,
                     cca_out["train_rs"], HERE/"cca_scatter.png")
    plot_permutation_null(null_top, cca_out["train_rs"][0],
                          HERE/"permutation_null.png")
    plot_top_loaders(cca_out["h_loadings"], cca_out["p_loadings"],
                     cell_cols, gene_index, 0,
                     HERE/"cca_loadings_axis1.png")
    plot_per_roi_similarity(sims, sections, HERE/"per_roi_cosine.png")

    top20 = sig_pos.nlargest(20, "r")
    if len(top20) > 0:
        fig, ax = plt.subplots(figsize=(8, max(4, 0.3*len(top20))))
        ax.barh(range(len(top20)), top20["r"], color="#d62728",
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels([f"{r.gene} ↔ {r.cell_type}" for _, r in top20.iterrows()],
                            fontsize=8)
        ax.set_xlabel("Pearson r across ROIs")
        ax.invert_yaxis()
        ax.set_title("slide2 — top 20 data-driven positive marker-celltype pairs (BH<0.05)",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(HERE/"top_discovered_pairs.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

    print(f"\nDone. Outputs in {HERE}:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
