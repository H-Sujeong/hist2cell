"""Mantel test: are the ROI×ROI distance matrices of Hist2Cell and
proteomics correlated? This is the most direct test of "the two modalities
see the same notion of which ROIs are similar / different", independent
of any section labelling.

Outputs (per slide, in proof_ver2/)
  mantel_scatter.png       paired off-diagonal distances + null comparison
  mantel_summary.csv       observed r (Pearson, Spearman) + permutation p
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr

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


def upper_tri(D):
    n = D.shape[0]
    iu = np.triu_indices(n, k=1)
    return D[iu]


def mantel(D1, D2, n_perm=N_PERM, seed=RANDOM_SEED):
    v1 = upper_tri(D1)
    v2 = upper_tri(D2)
    r_p, p_p_param = pearsonr(v1, v2)
    r_s, p_s_param = spearmanr(v1, v2)

    rng = np.random.default_rng(seed)
    n = D1.shape[0]
    null_pearson = np.zeros(n_perm)
    null_spearman = np.zeros(n_perm)
    for k in range(n_perm):
        perm = rng.permutation(n)
        D2_perm = D2[perm][:, perm]
        v2_perm = upper_tri(D2_perm)
        null_pearson[k], _ = pearsonr(v1, v2_perm)
        null_spearman[k], _ = spearmanr(v1, v2_perm)
    p_pearson_emp = float(np.mean(null_pearson >= r_p))
    p_spearman_emp = float(np.mean(null_spearman >= r_s))
    return {
        "pearson_r": float(r_p),
        "pearson_p_param": float(p_p_param),
        "pearson_p_perm": p_pearson_emp,
        "spearman_r": float(r_s),
        "spearman_p_param": float(p_s_param),
        "spearman_p_perm": p_spearman_emp,
        "null_pearson": null_pearson,
        "null_spearman": null_spearman,
    }


def plot_mantel(D1, D2, sections, res, slide_key, scfg, out_path):
    v1 = upper_tri(D1)
    v2 = upper_tri(D2)
    # color: same-section pairs vs cross-section pairs
    sec = np.array(sections)
    n = D1.shape[0]
    iu = np.triu_indices(n, k=1)
    pair_same = sec[iu[0]] == sec[iu[1]]

    fig = plt.figure(figsize=(13, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])
    ax_sc = fig.add_subplot(gs[0, 0])
    ax_null = fig.add_subplot(gs[0, 1])

    # scatter
    ax_sc.scatter(v1[~pair_same], v2[~pair_same], s=10, c="#888888",
                  alpha=0.4, edgecolor="none", label="cross-section pairs")
    ax_sc.scatter(v1[pair_same], v2[pair_same], s=14, c="#d62728",
                  alpha=0.8, edgecolor="none", label="same-section pairs")
    slope, intercept = np.polyfit(v1, v2, 1)
    xs = np.linspace(v1.min(), v1.max(), 50)
    ax_sc.plot(xs, slope*xs + intercept, c="black", linewidth=0.8, alpha=0.6)
    ax_sc.set_xlabel("Hist2Cell ROI×ROI distance (PC10)", fontsize=9)
    ax_sc.set_ylabel("Proteomics ROI×ROI distance (PC10)", fontsize=9)
    ax_sc.set_title(f"Paired distances  "
                    f"(Pearson r = {res['pearson_r']:+.3f},  "
                    f"Spearman ρ = {res['spearman_r']:+.3f})",
                    fontsize=10)
    ax_sc.legend(loc="best", fontsize=8)

    # null
    null = res["null_pearson"]
    ax_null.hist(null, bins=40, color="#bbbbbb", edgecolor="white", alpha=0.85,
                  label=f"permutation null (n={len(null)})")
    ax_null.axvline(res["pearson_r"], color="#d62728", linewidth=2,
                     label=f"observed Pearson r = {res['pearson_r']:+.3f}")
    ax_null.set_xlabel("Mantel Pearson r (Hist2Cell vs Proteomics distance)")
    ax_null.set_ylabel("permutation count")
    ax_null.set_title(f"Permutation null  "
                      f"(empirical p = {res['pearson_p_perm']:.4f})",
                      fontsize=10)
    ax_null.legend(loc="best", fontsize=8)

    fig.suptitle(f"{slide_key} — Mantel test: do the two modalities order "
                 f"ROI pairs by similarity in the same way?", fontsize=12)
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
    D_h = squareform(pdist(H_pcs))
    D_p = squareform(pdist(P_pcs))

    res = mantel(D_h, D_p)
    print(f"  Mantel Pearson  r = {res['pearson_r']:+.3f},  perm p = {res['pearson_p_perm']:.4f}")
    print(f"  Mantel Spearman ρ = {res['spearman_r']:+.3f},  perm p = {res['spearman_p_perm']:.4f}")
    print(f"  null Pearson    mean = {res['null_pearson'].mean():+.3f}, "
          f"95% range = [{np.percentile(res['null_pearson'],2.5):+.3f}, "
          f"{np.percentile(res['null_pearson'],97.5):+.3f}]")

    plot_mantel(D_h, D_p, sections, res, slide_key, scfg,
                 out_dir / "mantel_scatter.png")

    pd.DataFrame([{
        "slide": slide_key,
        "n_rois": len(common),
        "pearson_r": res["pearson_r"],
        "pearson_p_perm": res["pearson_p_perm"],
        "spearman_r": res["spearman_r"],
        "spearman_p_perm": res["spearman_p_perm"],
        "null_pearson_mean": float(res["null_pearson"].mean()),
        "null_pearson_95lo": float(np.percentile(res["null_pearson"], 2.5)),
        "null_pearson_95hi": float(np.percentile(res["null_pearson"], 97.5)),
    }]).to_csv(out_dir / "mantel_summary.csv", index=False)
    print(f"  saved → mantel_scatter.png, mantel_summary.csv")


def main():
    for k in SECTIONS:
        run_one(k)


if __name__ == "__main__":
    main()
