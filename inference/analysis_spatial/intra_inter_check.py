"""Quantitative check: do ROIs within the same section cluster tightly
in each modality? Is intra-section distance reliably smaller than
inter-section distance, in *both* Hist2Cell and proteomics?

Outputs (per slide, in proof_ver2/)
  distance_heatmap_hist2cell.png      ROI×ROI distance, section-ordered
  distance_heatmap_proteomics.png     same for proteomics
  silhouette_null.png                 observed silhouette vs 1000-perm null
  intra_inter_summary.csv             silhouette + intra/inter ratio + p
"""
import sys
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, N_PCS, RANDOM_SEED,
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

N_PERM = 1000


def compute_stats(X_pcs, sections):
    sil = float(silhouette_score(X_pcs, sections, metric="euclidean"))
    D = squareform(pdist(X_pcs))
    sec = np.array(sections)
    n = len(sec)
    same = sec[:, None] == sec[None, :]
    np.fill_diagonal(same, False)
    intra_pairs = D[same]
    inter_pairs = D[~same & ~np.eye(n, dtype=bool)]
    ratio = float(intra_pairs.mean() / inter_pairs.mean())
    return {
        "silhouette": sil,
        "intra_mean": float(intra_pairs.mean()),
        "inter_mean": float(inter_pairs.mean()),
        "intra_inter_ratio": ratio,
        "D": D,
    }


def permutation_null(X_pcs, sections, n_perm=N_PERM, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    sil_null = np.zeros(n_perm)
    ratio_null = np.zeros(n_perm)
    secs = np.array(sections)
    D = squareform(pdist(X_pcs))
    n = len(secs)
    eye_mask = np.eye(n, dtype=bool)
    for k in range(n_perm):
        permuted = rng.permutation(secs)
        try:
            sil_null[k] = silhouette_score(X_pcs, permuted, metric="euclidean")
        except ValueError:
            sil_null[k] = np.nan
        same = permuted[:, None] == permuted[None, :]
        same[eye_mask] = False
        intra_pairs = D[same]
        inter_pairs = D[~same & ~eye_mask]
        ratio_null[k] = intra_pairs.mean() / inter_pairs.mean()
    return sil_null, ratio_null


def plot_distance_heatmap(D, sections, scfg, title, out_path):
    sec_arr = np.array(sections)
    order = np.argsort([scfg["order"].index(s) for s in sec_arr])
    D_ord = D[order][:, order]
    sec_ord = sec_arr[order]

    fig = plt.figure(figsize=(10, 9))
    # Row 0 = legend strip, row 1 = top section colour strip, row 2 = main row
    # Col 0 = left section strip, col 1 = main heatmap, col 2 = colorbar
    gs = fig.add_gridspec(3, 3,
                          width_ratios=[0.035, 1, 0.045],
                          height_ratios=[0.08, 0.035, 1],
                          wspace=0.04, hspace=0.04)
    ax_legend = fig.add_subplot(gs[0, :])
    ax_top    = fig.add_subplot(gs[1, 1])
    ax_left   = fig.add_subplot(gs[2, 0])
    ax_main   = fig.add_subplot(gs[2, 1])
    ax_cb     = fig.add_subplot(gs[2, 2])

    im = ax_main.imshow(D_ord, cmap="magma_r", aspect="auto")
    ax_main.set_xticks([]); ax_main.set_yticks([])

    section_strip = np.array([[scfg["order"].index(s) for s in sec_ord]])
    cmap_strip = matplotlib.colors.ListedColormap(
        [scfg["colors"][s] for s in scfg["order"]])
    ax_top.imshow(section_strip, aspect="auto", cmap=cmap_strip,
                  vmin=0, vmax=len(scfg["order"])-1)
    ax_top.set_xticks([]); ax_top.set_yticks([])
    ax_left.imshow(section_strip.T, aspect="auto", cmap=cmap_strip,
                   vmin=0, vmax=len(scfg["order"])-1)
    ax_left.set_xticks([]); ax_left.set_yticks([])

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("Euclidean distance (PC10 space)", fontsize=9)

    handles = [plt.Rectangle((0, 0), 1, 1, color=scfg["colors"][s],
                              label=scfg["labels"][s]) for s in scfg["order"]]
    ax_legend.axis("off")
    ax_legend.legend(handles=handles, loc="center",
                      ncol=len(scfg["order"]), fontsize=9, frameon=False,
                      handletextpad=0.4, columnspacing=1.5)

    fig.suptitle(title, fontsize=12, y=0.99)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_silhouette_null(sil_h, sil_p, null_h, null_p, slide_key, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, observed, null, title in [
            (axes[0], sil_h, null_h, "Hist2Cell"),
            (axes[1], sil_p, null_p, "Proteomics")]:
        null = null[~np.isnan(null)]
        ax.hist(null, bins=40, color="#bbbbbb", edgecolor="white", alpha=0.85,
                label=f"permutation null (n={len(null)})")
        ax.axvline(observed, color="#d62728", linewidth=2,
                   label=f"observed = {observed:+.3f}")
        p_emp = float(np.mean(null >= observed))
        ax.set_xlabel("silhouette score")
        ax.set_ylabel("permutation count")
        ax.set_title(f"{title}  (empirical p = {p_emp:.4f})", fontsize=11)
        ax.legend(loc="best", fontsize=8)
    fig.suptitle(f"{slide_key} — section-label silhouette: observed vs null",
                 fontsize=12)
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

    Hs = StandardScaler().fit_transform(H)
    Ps = StandardScaler().fit_transform(P)
    H_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Hs)
    P_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Ps)

    stats_h = compute_stats(H_pcs, sections)
    stats_p = compute_stats(P_pcs, sections)
    print(f"  Hist2Cell:  silhouette = {stats_h['silhouette']:+.3f}, "
          f"intra/inter = {stats_h['intra_inter_ratio']:.3f}  "
          f"(intra mean {stats_h['intra_mean']:.2f}, "
          f"inter mean {stats_h['inter_mean']:.2f})")
    print(f"  Proteomics: silhouette = {stats_p['silhouette']:+.3f}, "
          f"intra/inter = {stats_p['intra_inter_ratio']:.3f}  "
          f"(intra mean {stats_p['intra_mean']:.2f}, "
          f"inter mean {stats_p['inter_mean']:.2f})")

    plot_distance_heatmap(
        stats_h["D"], sections, scfg,
        f"{slide_key} Hist2Cell — ROI×ROI distance (section-ordered)",
        out_dir / "distance_heatmap_hist2cell.png")
    plot_distance_heatmap(
        stats_p["D"], sections, scfg,
        f"{slide_key} Proteomics — ROI×ROI distance (section-ordered)",
        out_dir / "distance_heatmap_proteomics.png")

    print(f"  [perm] {N_PERM} shuffles …", flush=True)
    t0 = time.time()
    sil_null_h, ratio_null_h = permutation_null(H_pcs, sections)
    sil_null_p, ratio_null_p = permutation_null(P_pcs, sections)
    print(f"  perm time = {time.time()-t0:.1f}s")

    p_sil_h = float(np.mean(sil_null_h[~np.isnan(sil_null_h)]
                            >= stats_h["silhouette"]))
    p_sil_p = float(np.mean(sil_null_p[~np.isnan(sil_null_p)]
                            >= stats_p["silhouette"]))
    p_ratio_h = float(np.mean(ratio_null_h <= stats_h["intra_inter_ratio"]))
    p_ratio_p = float(np.mean(ratio_null_p <= stats_p["intra_inter_ratio"]))
    print(f"  silhouette empirical p:   Hist2Cell {p_sil_h:.4f},  "
          f"Proteomics {p_sil_p:.4f}")
    print(f"  intra/inter empirical p:  Hist2Cell {p_ratio_h:.4f},  "
          f"Proteomics {p_ratio_p:.4f}")

    plot_silhouette_null(stats_h["silhouette"], stats_p["silhouette"],
                          sil_null_h, sil_null_p, slide_key,
                          out_dir / "silhouette_null.png")

    summary = pd.DataFrame([
        {"modality": "Hist2Cell",  **{k: v for k, v in stats_h.items() if k != "D"},
         "p_silhouette": p_sil_h, "p_intra_inter": p_ratio_h,
         "null_sil_mean": float(np.nanmean(sil_null_h)),
         "null_sil_95hi": float(np.nanpercentile(sil_null_h, 97.5)),
         "null_ratio_mean": float(np.nanmean(ratio_null_h)),
         "null_ratio_95lo": float(np.nanpercentile(ratio_null_h, 2.5))},
        {"modality": "Proteomics", **{k: v for k, v in stats_p.items() if k != "D"},
         "p_silhouette": p_sil_p, "p_intra_inter": p_ratio_p,
         "null_sil_mean": float(np.nanmean(sil_null_p)),
         "null_sil_95hi": float(np.nanpercentile(sil_null_p, 97.5)),
         "null_ratio_mean": float(np.nanmean(ratio_null_p)),
         "null_ratio_95lo": float(np.nanpercentile(ratio_null_p, 2.5))},
    ])
    summary.to_csv(out_dir / "intra_inter_summary.csv", index=False)
    print(f"  saved → distance_heatmap_*.png, silhouette_null.png, "
          f"intra_inter_summary.csv")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
