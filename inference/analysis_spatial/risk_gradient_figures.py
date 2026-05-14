"""Re-render CCA scatter and Mantel scatter with per-ROI risk-score
gradient overlay.  Replaces the existing PNGs in proof_ver2/.

CCA scatter — fill color = risk score (magma colormap),
              marker shape = section, separate panel per canonical axis.
Mantel scatter — point color = |Δ risk| between the two ROIs in the pair
              (viridis_r colormap; brighter = more similar risk).
"""
import pickle
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
sys.path.insert(0, str(HERE))

from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, run_cca, N_PCS, RANDOM_SEED,
)


SECTIONS = {
    "1_085_12": {
        "labels": {"a": "High-risk Tumor", "b": "Low-risk Tumor",
                   "c": "High-risk T-cell", "d": "Low-risk T-cell",
                   "t": "Middle-risk Tumor (ctrl)"},
        "order": ["a", "b", "c", "d", "t"],
        "markers": {"a": "o", "b": "s", "c": "^", "d": "D", "t": "P"},
        "prefix": "abcdt",
        "pred_csv": "/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv",
        "colors": {"a": "#d62728", "b": "#1f77b4", "c": "#2ca02c",
                   "d": "#9467bd", "t": "#7f7f7f"},
    },
    "1_152_19": {
        "labels": {"e": "High-risk Tumor", "f": "Low-risk Tumor",
                   "g": "High-risk T-cell", "h": "Low-risk T-cell",
                   "v": "Middle-risk Tumor (ctrl)"},
        "order": ["e", "f", "g", "h", "v"],
        "markers": {"e": "o", "f": "s", "g": "^", "h": "D", "v": "P"},
        "prefix": "efghv",
        "pred_csv": "/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv",
        "colors": {"e": "#d62728", "f": "#1f77b4", "g": "#2ca02c",
                   "h": "#9467bd", "v": "#7f7f7f"},
    },
}


def load_risk_scores(slide_key):
    pkl = HERE / slide_key / f"{slide_key}_ROI_groups_risk_scores.pkl"
    with open(pkl, "rb") as f:
        return pickle.load(f)


def plot_cca_scatter_risk(Hc, Pc, common, sections, train_rs, risk, scfg,
                           slide_key, out_path):
    """3-panel scatter, color = risk gradient (magma),
    marker shape = section. Colorbar in its own gridspec column."""
    n_comp = Hc.shape[1]
    risk_arr = np.array([risk[t] for t in common])
    vmin, vmax = risk_arr.min(), risk_arr.max()
    fig = plt.figure(figsize=(6.0 * n_comp + 1.2, 6.0))
    gs = fig.add_gridspec(1, n_comp + 1,
                          width_ratios=[1]*n_comp + [0.04],
                          wspace=0.30)
    axes = [fig.add_subplot(gs[0, i]) for i in range(n_comp)]
    cax = fig.add_subplot(gs[0, n_comp])
    for i, ax in enumerate(axes):
        sc = None
        for s in scfg["order"]:
            mask = [j for j, sec in enumerate(sections) if sec == s]
            if not mask: continue
            sc = ax.scatter(Hc[mask, i], Pc[mask, i],
                            c=risk_arr[mask], cmap="magma",
                            vmin=vmin, vmax=vmax,
                            marker=scfg["markers"][s], s=110,
                            edgecolor="black", linewidth=0.4, alpha=0.95)
        xs = np.linspace(Hc[:, i].min(), Hc[:, i].max(), 50)
        slope, intercept = np.polyfit(Hc[:, i], Pc[:, i], 1)
        ax.plot(xs, slope*xs + intercept, c="black", linewidth=0.7, alpha=0.5)
        ax.set_xlabel(f"Hist2Cell canonical {i+1}", fontsize=9)
        ax.set_ylabel(f"Proteomics canonical {i+1}", fontsize=9)
        ax.set_title(f"Canonical pair {i+1}: r = {train_rs[i]:+.3f}",
                     fontsize=11)
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label("ROI risk score", fontsize=9)
    # section marker legend at bottom
    handles = [plt.scatter([], [], marker=scfg["markers"][s], s=110,
                            c="white", edgecolor="black", linewidth=0.4,
                            label=scfg["labels"][s])
               for s in scfg["order"]]
    fig.legend(handles=handles, loc="lower center",
                bbox_to_anchor=(0.5, -0.02),
                ncol=len(scfg["order"]), fontsize=8.5, frameon=True,
                framealpha=0.9)
    fig.suptitle(f"CCA — {slide_key}: canonical pairs coloured by risk score "
                 f"(marker = section)", fontsize=12, y=0.99)
    fig.subplots_adjust(bottom=0.14, top=0.92, left=0.05, right=0.95)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_mantel_risk(D_h, D_p, sections, common, risk, scfg, slide_key,
                      out_path, n_perm=1000, seed=RANDOM_SEED):
    """Mantel scatter coloured by |Δ risk| between ROI pairs."""
    n = D_h.shape[0]
    iu = np.triu_indices(n, k=1)
    v1 = D_h[iu]
    v2 = D_p[iu]
    risk_arr = np.array([risk[t] for t in common])
    drisk = np.abs(risk_arr[iu[0]] - risk_arr[iu[1]])

    r_p, _ = pearsonr(v1, v2)
    r_s, _ = spearmanr(v1, v2)

    # null
    rng = np.random.default_rng(seed)
    null_r = np.zeros(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(n)
        D_p_perm = D_p[perm][:, perm]
        null_r[k], _ = pearsonr(v1, D_p_perm[iu])
    p_perm = float(np.mean(null_r >= r_p))

    fig = plt.figure(figsize=(13, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])
    ax_sc = fig.add_subplot(gs[0, 0])
    ax_null = fig.add_subplot(gs[0, 1])

    order = np.argsort(drisk)[::-1]  # plot larger Δrisk first so smaller on top
    sc = ax_sc.scatter(v1[order], v2[order], c=drisk[order], cmap="viridis_r",
                        s=11, alpha=0.7, edgecolor="none")
    slope, intercept = np.polyfit(v1, v2, 1)
    xs = np.linspace(v1.min(), v1.max(), 50)
    ax_sc.plot(xs, slope*xs + intercept, c="black", linewidth=0.8, alpha=0.6)
    ax_sc.set_xlabel("Hist2Cell ROI×ROI distance (PC10)", fontsize=9)
    ax_sc.set_ylabel("Proteomics ROI×ROI distance (PC10)", fontsize=9)
    ax_sc.set_title(f"Paired distances  "
                    f"(Pearson r = {r_p:+.3f},  Spearman ρ = {r_s:+.3f})",
                    fontsize=10)
    cb = fig.colorbar(sc, ax=ax_sc, shrink=0.85, pad=0.02)
    cb.set_label("|Δ risk score| between ROI pair", fontsize=8.5)

    ax_null.hist(null_r, bins=40, color="#bbbbbb", edgecolor="white",
                  alpha=0.85, label=f"permutation null (n={n_perm})")
    ax_null.axvline(r_p, color="#d62728", linewidth=2,
                     label=f"observed Pearson r = {r_p:+.3f}")
    ax_null.set_xlabel("Mantel Pearson r")
    ax_null.set_ylabel("permutation count")
    ax_null.set_title(f"Permutation null  (empirical p = {p_perm:.4f})",
                      fontsize=10)
    ax_null.legend(loc="best", fontsize=8)

    fig.suptitle(f"{slide_key} — Mantel test: ROI pairs coloured by risk-score "
                 f"difference  (bright = similar risk, dark = different risk)",
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
    risk = load_risk_scores(slide_key)
    print(f"\n=== {slide_key} (N={len(common)}) ===")
    print(f"  risk score: min={min(risk.values()):.3f} "
          f"max={max(risk.values()):.3f} "
          f"mean={np.mean(list(risk.values())):.3f}")

    # CCA scatter
    cca_out = run_cca(H, P)
    plot_cca_scatter_risk(cca_out["Hc"], cca_out["Pc"], common, sections,
                           cca_out["train_rs"], risk, scfg, slide_key,
                           out_dir / "cca_scatter.png")
    print(f"  ✓ cca_scatter.png  (risk-gradient overlay)")

    # Mantel scatter
    Hs = StandardScaler().fit_transform(H)
    Ps = StandardScaler().fit_transform(P)
    H_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Hs)
    P_pcs = PCA(n_components=N_PCS, random_state=RANDOM_SEED).fit_transform(Ps)
    D_h = squareform(pdist(H_pcs))
    D_p = squareform(pdist(P_pcs))
    plot_mantel_risk(D_h, D_p, sections, common, risk, scfg, slide_key,
                      out_dir / "mantel_scatter.png")
    print(f"  ✓ mantel_scatter.png  (|Δrisk| overlay)")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
