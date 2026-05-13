"""Two-figure quick check: do Hist2Cell and proteomics see the same ROI structure?

Figure A — side-by-side UMAP (one per modality), coloured by section.
Figure B — joint UMAP after CCA alignment: each ROI appears as two points
           (circle = Hist2Cell, triangle = proteomics) with a thin line
           connecting matched ROIs. Short lines = strong cross-modality
           agreement.

Outputs land in each slide's proof_ver2/ next to the existing figures.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, N_PCS, RANDOM_SEED,
)
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
import umap


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
    },
}


def fit_umap(X, n_neighbors=10, min_dist=0.3, seed=RANDOM_SEED):
    return umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                      n_components=2, random_state=seed).fit_transform(X)


def side_by_side_umap(H, P, sections, slide_name, scfg, out_path):
    Hs = StandardScaler().fit_transform(H)
    Ps = StandardScaler().fit_transform(P)
    H_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Hs)
    P_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Ps)
    H_um = fit_umap(H_pcs)
    P_um = fit_umap(P_pcs)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, emb, title in [(axes[0], H_um, "Hist2Cell UMAP"),
                            (axes[1], P_um, "Proteomics UMAP")]:
        for s in scfg["order"]:
            mask = [i for i, sec in enumerate(sections) if sec == s]
            if mask:
                ax.scatter(emb[mask, 0], emb[mask, 1], s=90,
                           c=scfg["colors"][s], edgecolor="black",
                           linewidth=0.5, label=scfg["labels"][s], alpha=0.85)
        ax.set_xlabel("UMAP 1", fontsize=9)
        ax.set_ylabel("UMAP 2", fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
    axes[0].legend(loc="best", fontsize=8)
    fig.suptitle(f"{slide_name} — same ROIs embedded by each modality "
                 f"(side-by-side UMAP)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def joint_umap_after_cca(H, P, sections, common, slide_name, scfg, out_path):
    Hs = StandardScaler().fit_transform(H)
    Ps = StandardScaler().fit_transform(P)
    H_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Hs)
    P_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Ps)
    # align with CCA
    cca = CCA(n_components=min(N_PCS, len(common)-1), max_iter=1000)
    Hc, Pc = cca.fit_transform(H_pcs, P_pcs)
    # stack: rows = Hist2Cell ROIs then Proteomics ROIs (same order)
    stacked = np.vstack([Hc, Pc])
    um = fit_umap(stacked, n_neighbors=8, min_dist=0.25)
    n = len(common)
    H_um = um[:n]; P_um = um[n:]

    # Compute average pair distance for sanity number
    dists = np.linalg.norm(H_um - P_um, axis=1)
    median_dist = float(np.median(dists))

    fig, ax = plt.subplots(figsize=(9, 7))
    # connecting lines (light grey)
    for i in range(n):
        ax.plot([H_um[i, 0], P_um[i, 0]], [H_um[i, 1], P_um[i, 1]],
                color="#888888", linewidth=0.5, alpha=0.55, zorder=1)
    # points
    for s in scfg["order"]:
        mask = [i for i, sec in enumerate(sections) if sec == s]
        if mask:
            ax.scatter(H_um[mask, 0], H_um[mask, 1], s=110,
                       c=scfg["colors"][s], edgecolor="black",
                       linewidth=0.5, alpha=0.95, marker="o",
                       label=f"{scfg['labels'][s]} (Hist2Cell)", zorder=2)
            ax.scatter(P_um[mask, 0], P_um[mask, 1], s=110,
                       c=scfg["colors"][s], edgecolor="black",
                       linewidth=0.5, alpha=0.6, marker="^",
                       zorder=2)
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{slide_name} — joint UMAP after CCA alignment\n"
                 f"○ = Hist2Cell, △ = Proteomics, grey lines = same ROI  "
                 f"| median pair dist = {median_dist:.2f}", fontsize=11)
    ax.legend(loc="best", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return median_dist


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
    print(f"[{slide_key}] N={len(common)} ROIs, H={H.shape}, P={P.shape}")

    side_by_side_umap(H, P, sections, slide_key, scfg,
                      out_dir / "umap_side_by_side.png")
    median_d = joint_umap_after_cca(H, P, sections, common, slide_key, scfg,
                                     out_dir / "umap_joint_paired.png")
    print(f"  saved → umap_side_by_side.png  umap_joint_paired.png  "
          f"(median pair dist={median_d:.2f})")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
