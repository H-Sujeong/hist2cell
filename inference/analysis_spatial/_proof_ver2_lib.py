"""Shared data-driven cross-modality validation library for proof_ver2.

Approach (deliberately ignores the collaborator's pre-selected marker
panel; rederives everything from the slide's own Hist2Cell × proteomics
matrices):

  1. CCA on PCA-reduced features → top canonical correlation between
     the two modalities, with permutation null to gauge how far the
     observed value sits above the random-pairing baseline.
  2. Data-driven marker discovery — for every (proteomics gene,
     Hist2Cell cell type) pair, compute Pearson r across the ROIs that
     have both modalities; rank by signed r with BH-FDR control.
  3. Per-ROI cross-modality similarity — for each ROI, build a
     proteomics-derived cell-type score using the markers discovered in
     step 2, then compute cosine similarity against the Hist2Cell
     vector for that same ROI.  Average across ROIs.

The library exposes one `run_proof_ver2(...)` entry point used by each
slide's `core_proofs_v2.py`.  All numeric output goes into CSVs in the
slide's proof_ver2/ folder; plots are PNGs of CCA scatter, permutation
null histogram, top-pair heatmap, and per-ROI similarity bars.

Limitations made explicit in the slide-level summary:
  - N=46-48 ROIs is small; CCA can overfit even after PCA reduction,
    hence the permutation test (any genuine canonical structure should
    sit far outside the shuffled null).
  - Discovered marker pairs are post-hoc — BH-FDR controls the FDR but
    selection within rank-ordered lists still inflates effect sizes.
  - Lung-trained Hist2Cell limitation still applies: discovered
    matches may reflect *morphology-shared* features rather than
    actual breast cell-type identity.
"""
from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


NPY_TILE      = 512
ROI_PATCH     = 1024
ROI_OFFSETS   = [(0, 0), (512, 0), (0, 512), (512, 512)]
DETECT_MIN    = 0.50          # gene must be detected in ≥50% of ROIs
N_PCS         = 10            # PCA components per modality before CCA
N_CCA_COMP    = 3             # canonical components to extract
N_PERM        = 1000          # permutation replicates
RANDOM_SEED   = 42

GROUPS_CSV = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")
GG_MATRIX  = Path("/home/sjhong/hist2cell/inference/analysis_spatial/report.gg_matrix (1).tsv")


@dataclass
class SlideConfig:
    name: str
    pred_csv: Path
    roi_pkl: Path
    npy: Path
    section_label: Dict[str, str]
    section_color: Dict[str, str]
    sample_section_prefixes: str   # e.g. "abcdt" for slide1, "efghv" for slide2
    out_dir: Path


# ── Hist2Cell ROI signature build (same logic as cell_typing/) ────────

def sort_tubes(keys):
    return sorted(keys, key=lambda t: (t[0],
                                       int(t[1:]) if t[1:].isdigit() else 0))


def per_tile_signature(npy_tiles, XY, P):
    n, m = len(npy_tiles), P.shape[1]
    sig = np.zeros((n, m), dtype=np.float32)
    n_spots = np.zeros(n, dtype=int)
    tree = cKDTree(XY)
    r = NPY_TILE * np.sqrt(2)
    for i, (x, y) in enumerate(npy_tiles):
        cand = tree.query_ball_point((x + NPY_TILE/2, y + NPY_TILE/2), r=r)
        if not cand:
            continue
        cx, cy = XY[cand, 0], XY[cand, 1]
        keep = ((cx >= x) & (cx < x + NPY_TILE) &
                (cy >= y) & (cy < y + NPY_TILE))
        idx = [cand[j] for j, k in enumerate(keep) if k]
        if idx:
            sig[i] = P[idx].mean(axis=0)
            n_spots[i] = len(idx)
    return sig, n_spots


def build_roi_signatures(cfg: SlideConfig) -> pd.DataFrame:
    preds = pd.read_csv(cfg.pred_csv)
    cell_cols = [c for c in preds.columns if c not in ("spot_id", "X", "Y")]
    P = preds[cell_cols].values.astype(np.float32)
    XY = preds[["X", "Y"]].values.astype(np.float64)
    with open(cfg.roi_pkl, "rb") as f:
        roi = pickle.load(f)
    npy = np.load(cfg.npy).astype(np.float64)
    tile_sig, n_spots_tile = per_tile_signature(npy, XY, P)
    nset = {tuple(v): i for i, v in enumerate(npy.tolist())}

    rows = []
    for tid in sort_tubes(roi.keys()):
        idx, spots_total = [], 0
        for px, py in roi[tid]:
            for dx, dy in ROI_OFFSETS:
                k = (float(px+dx), float(py+dy))
                if k in nset:
                    j = nset[k]
                    idx.append(j); spots_total += int(n_spots_tile[j])
        sig = tile_sig[idx].mean(axis=0) if idx else np.zeros(P.shape[1], dtype=np.float32)
        rec = {"tube_id": tid, "section": tid[0], "n_spots": int(spots_total),
               **{c: float(sig[i]) for i, c in enumerate(cell_cols)}}
        rows.append(rec)
    return pd.DataFrame(rows), cell_cols


# ── Proteomics matrix load + filter ──────────────────────────────────

def load_proteomics_matrix(cfg: SlideConfig) -> Tuple[pd.DataFrame, list]:
    m = pd.read_csv(GG_MATRIX, sep="\t")
    slide_cols = [c for c in m.columns[3:] if c[0] in cfg.sample_section_prefixes]
    raw = m[slide_cols].copy().replace(0, np.nan)
    log2 = np.log2(raw); log2.index = m["Genes"].astype(str)
    # quality filter
    detect = log2.notna().mean(axis=1)
    keep = detect >= DETECT_MIN
    log2_f = log2.loc[keep].copy()
    # median-impute remaining NaN
    medians = log2_f.median(axis=1)
    log2_f = log2_f.apply(lambda col: col.fillna(medians))
    return log2_f, slide_cols      # rows = genes, cols = samples


# ── Alignment + matrices ─────────────────────────────────────────────

def align_modalities(sig_df, log2_f, slide_cols, cell_cols):
    common = [t for t in sig_df["tube_id"] if t in slide_cols]
    sig_aligned = sig_df.set_index("tube_id").loc[common]
    H = sig_aligned[cell_cols].values.astype(float)         # [N, 80]
    P = log2_f[common].T.values.astype(float)               # [N, n_genes]
    return common, H, P, sig_aligned, list(log2_f.index)


# ── CCA + permutation null ───────────────────────────────────────────

def run_cca(H, P, n_pcs=N_PCS, n_components=N_CCA_COMP, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    n_pcs = min(n_pcs, H.shape[0]-1, H.shape[1], P.shape[1])
    pca_h = PCA(n_components=n_pcs, random_state=seed)
    pca_p = PCA(n_components=n_pcs, random_state=seed)
    H_pc = pca_h.fit_transform(H)
    P_pc = pca_p.fit_transform(P)
    cca = CCA(n_components=n_components, max_iter=1000)
    cca.fit(H_pc, P_pc)
    Hc, Pc = cca.transform(H_pc, P_pc)
    rs = []
    for i in range(n_components):
        rs.append(float(pearsonr(Hc[:, i], Pc[:, i])[0]))
    # loadings back to original feature space:
    #   H feature loading = pca_h.components_.T @ cca.x_loadings_ → [80, n_comp]
    h_load = pca_h.components_.T @ cca.x_loadings_
    p_load = pca_p.components_.T @ cca.y_loadings_
    return {
        "train_rs": rs,
        "pca_h_var": pca_h.explained_variance_ratio_,
        "pca_p_var": pca_p.explained_variance_ratio_,
        "h_loadings": h_load,         # [80, n_components]
        "p_loadings": p_load,         # [n_genes, n_components]
        "Hc": Hc, "Pc": Pc, "n_pcs": n_pcs,
    }


def permutation_null(H, P, n_perm=N_PERM, n_pcs=N_PCS, n_components=N_CCA_COMP,
                     seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    n = H.shape[0]
    n_pcs = min(n_pcs, n-1, H.shape[1], P.shape[1])
    null_top = np.zeros(n_perm)
    for it in range(n_perm):
        perm = rng.permutation(n)
        try:
            pca_h = PCA(n_components=n_pcs, random_state=seed)
            pca_p = PCA(n_components=n_pcs, random_state=seed)
            H_pc = pca_h.fit_transform(H)
            P_pc = pca_p.fit_transform(P[perm])
            cca = CCA(n_components=n_components, max_iter=1000)
            cca.fit(H_pc, P_pc)
            Hc, Pc = cca.transform(H_pc, P_pc)
            null_top[it] = float(pearsonr(Hc[:, 0], Pc[:, 0])[0])
        except Exception:
            null_top[it] = np.nan
    return null_top


# ── data-driven marker-celltype pair discovery ───────────────────────

def discover_marker_pairs(H, P, cell_cols, gene_index,
                          top_k_per_type=5, p_threshold=0.05):
    """For every (gene × cell type) pair compute Pearson r across N ROIs.
    Return a sorted DataFrame with BH-FDR per cell type."""
    n, m_h = H.shape
    _, m_p = P.shape
    rows = []
    for j, ctype in enumerate(cell_cols):
        h_col = H[:, j]
        if h_col.std() == 0:
            continue
        # vectorized pearson r for h_col vs each gene column
        H_c = h_col - h_col.mean()
        P_c = P - P.mean(axis=0, keepdims=True)
        # avoid div by zero
        h_std = h_col.std(ddof=1)
        p_std = P.std(axis=0, ddof=1)
        rs = (H_c @ P_c) / ((n-1) * h_std * np.where(p_std==0, 1, p_std))
        # p-values (Fisher transform approximation)
        with np.errstate(divide='ignore', invalid='ignore'):
            t_stat = rs * np.sqrt((n - 2) / np.clip(1 - rs**2, 1e-12, None))
            from scipy.stats import t as student_t
            ps = 2 * (1 - student_t.cdf(np.abs(t_stat), df=n-2))
        # take top_k by signed r (separately keep top positive and top negative)
        order = np.argsort(-rs)   # descending — top positive first
        for k in range(min(top_k_per_type, m_p)):
            gi = order[k]
            rows.append({"cell_type": ctype, "gene": gene_index[gi],
                         "r": float(rs[gi]), "p": float(ps[gi]),
                         "rank": k+1, "direction": "pos"})
        order_neg = np.argsort(rs)   # ascending — top negative first
        for k in range(min(top_k_per_type, m_p)):
            gi = order_neg[k]
            rows.append({"cell_type": ctype, "gene": gene_index[gi],
                         "r": float(rs[gi]), "p": float(ps[gi]),
                         "rank": k+1, "direction": "neg"})
    df = pd.DataFrame(rows)
    df["p_bh"] = multipletests(df["p"].clip(lower=1e-300), method="fdr_bh")[1]
    return df.sort_values(["direction", "cell_type", "rank"]).reset_index(drop=True)


# ── Per-ROI cross-modality similarity ────────────────────────────────

def per_roi_cosine(H, P, gene_index, pair_df, top_n_per_type=3):
    """Build proteomics-derived cell-type score: for each cell type,
    average log2 intensity of its top-N positively-correlated genes
    (data-driven marker set).  Then compute per-ROI cosine similarity
    against Hist2Cell vector."""
    cell_cols = list(pair_df["cell_type"].unique())
    n_cells = len(cell_cols)
    n = H.shape[0]
    # Hist2Cell sub-matrix (only those cell types we have markers for)
    h_idx_map = {c: i for i, c in enumerate(cell_cols)}
    proteomics_score = np.zeros((n, n_cells), dtype=float)
    pos_pairs = pair_df[pair_df["direction"] == "pos"]
    used_genes_per_type = {}
    for ctype, sub in pos_pairs.groupby("cell_type"):
        top = sub.nlargest(top_n_per_type, "r")
        genes = top["gene"].tolist()
        used_genes_per_type[ctype] = genes
        gene_pos = [gene_index.index(g) for g in genes if g in gene_index]
        if gene_pos:
            proteomics_score[:, h_idx_map[ctype]] = P[:, gene_pos].mean(axis=1)
    return proteomics_score, used_genes_per_type, cell_cols


def roi_cosine_similarity(H, P_score, hist_cell_cols, marker_cell_cols):
    """For each ROI compute cosine similarity between
    Hist2Cell[ROI, marker_cell_cols] and P_score[ROI, marker_cell_cols]."""
    # restrict Hist2Cell to the same cell types we have markers for
    idx = [hist_cell_cols.index(c) for c in marker_cell_cols if c in hist_cell_cols]
    H_sub = H[:, idx]
    P_sub = P_score[:, :len(idx)]
    sims = []
    for i in range(H_sub.shape[0]):
        h, p = H_sub[i], P_sub[i]
        if np.linalg.norm(h) == 0 or np.linalg.norm(p) == 0:
            sims.append(np.nan); continue
        sims.append(float(np.dot(h, p) / (np.linalg.norm(h) * np.linalg.norm(p))))
    return np.array(sims)
