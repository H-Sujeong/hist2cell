"""Hist2Cell ROI-level analysis for slide1 (1_085_12).

Two PNG families:
  section_*.png — ROI-level plots (47 tube nodes coloured by descriptive
                  section label; uses the cropped tissue mask as backdrop
                  so each tube sits in its slide-anatomical context)
  spatial_*.png — Hist2Cell spot-level heatmaps in the style of
                  analysis_filtered/ (dense scatter coloured by abundance,
                  filtered to the cropped tissue-mask X range so the plot
                  fits the slide silhouette)

Section labels (per user spec — descriptive, not prefix):
  a → High-risk Tumor
  b → Low-risk Tumor
  c → High-risk T-cell
  d → Low-risk T-cell
  t → Middle-risk Tumor (control)

Inputs
  ../1_085_12_ROI_groups.pkl           ROI tubes (47 / 181 patches)
  ../meteo_1_085_12_coords.npy         5,227 candidate 512-tile top-lefts
  ../meteo_1_085_12_coords_cropped.npy filtered npy (X-range applied)
  ../tissue_mask.png                   thumbnail tissue mask (full slide)
  ../tissue_mask_cropped.png           same mask zero-ed outside X∈[30000, 175000]
  ../../analysis/cell_type_groups.csv  cell type → group + strict/broad flags
  /home/sjhong/hist2cell/inference/slide1_085_12_v2/
       predictions.csv                  35,821 spots × 80 cell types
       slide1_085_12_coords.h5          WSI dims + tile metadata

Outputs (this folder)
  roi_signatures.csv               per-tube 47 × (4 + 80 + 3)
  roi_spot_counts.csv              per-tube n_patches / n_tiles / n_spots
  section_stats.csv                Wilcoxon a vs b / c vs d
  per_celltype_wilcoxon.csv        80-row table with BH-FDR (a vs b)
  proteomics_marker_hypotheses.csv pre-registered marker matches
  moran_within_roi.csv             80×80 Moran R, ROI subgraph
  moran_slide_wide.csv             80×80 Moran R, filtered Hist2Cell graph

  section_boxplots.png             3 scores per section (47 tubes)
  section_top10_celltypes.png      47 ROI dots × 10 top cell types over mask
  section_group_heatmaps.png       47 ROI dots × (10 groups + 2 proxies)
  section_immune_vs_epithelial.png 47 ROI dots × 3 scores
  section_subgraph.png             47-tube subgraph layout over mask

  spatial_top10_celltypes.png      Hist2Cell-spot scatter × 10 top types
  spatial_group_heatmaps.png       Hist2Cell-spot scatter × 12 group panels
  spatial_immune_vs_epithelial.png Hist2Cell-spot scatter × 3 scores
  moran_r_clustermap.png           80×80 clustermap, ROI subgraph
  moran_r_clustermap_slide.png     80×80 clustermap, slide filtered graph
"""
import pickle
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ---- paths ----

HERE      = Path(__file__).resolve().parent
PARENT    = HERE.parent
PRED_CSV  = Path("/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv")
COORDS_H5 = Path("/home/sjhong/hist2cell/inference/slide1_085_12_v2/slide1_085_12_coords.h5")
GROUPS_CSV = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")

ROI_PKL   = PARENT / "1_085_12_ROI_groups.pkl"
NPY_FULL  = PARENT / "meteo_1_085_12_coords.npy"
NPY_CROP  = PARENT / "meteo_1_085_12_coords_cropped.npy"
MASK_FULL = PARENT / "tissue_mask.png"
MASK_CROP = PARENT / "tissue_mask_cropped.png"


# ---- constants ----

NPY_TILE   = 512        # candidate-tile size in level-0 px
ROI_PATCH  = 1024       # proteomics ROI patch size (270 μm)
ROI_OFFSETS = [(0, 0), (512, 0), (0, 512), (512, 512)]
MORAN_KNN  = 12         # ROI-tube subgraph k
MORAN_KNN_SLIDE = 20    # slide-wide Hist2Cell graph k

X_KEEP = (30000, 175000)   # cropped-mask X range (user-defined)

SECTION_LABEL = {
    "a": "High-risk Tumor",
    "b": "Low-risk Tumor",
    "c": "High-risk T-cell",
    "d": "Low-risk T-cell",
    "t": "Middle-risk Tumor (ctrl)",
}
SECTION_COLOR = {
    "a": "#d62728", "b": "#1f77b4", "c": "#2ca02c",
    "d": "#9467bd", "t": "#7f7f7f",
}
SECTION_ORDER = ["a", "b", "c", "d", "t"]


# ---- I/O ----

def load_inputs():
    preds = pd.read_csv(PRED_CSV)
    cell_cols = [c for c in preds.columns if c not in ("spot_id", "X", "Y")]
    P_full = preds[cell_cols].values.astype(np.float32)
    XY_full = preds[["X", "Y"]].values.astype(np.float64)
    with open(ROI_PKL, "rb") as f:
        roi = pickle.load(f)
    npy_full = np.load(NPY_FULL).astype(np.float64)
    npy_crop = np.load(NPY_CROP).astype(np.float64)
    mask_full = np.array(Image.open(MASK_FULL))
    mask_crop = np.array(Image.open(MASK_CROP))
    groups = pd.read_csv(GROUPS_CSV)
    with h5py.File(COORDS_H5, "r") as f:
        meta = {k: v for k, v in f["metadata"].attrs.items()}
    W_lvl0 = int(meta["wsi_width"])
    H_lvl0 = int(meta["wsi_height"])
    return (preds, cell_cols, P_full, XY_full,
            roi, npy_full, npy_crop,
            mask_full, mask_crop, groups,
            W_lvl0, H_lvl0)


def sort_tubes(keys):
    return sorted(keys, key=lambda t: (t[0],
                                       int(t[1:]) if t[1:].isdigit() else 0))


# ---- Hist2Cell spot filtering (matches cropped mask) ----

def filter_spots_to_mask_range(XY, P):
    keep = (XY[:, 0] >= X_KEEP[0]) & (XY[:, 0] <= X_KEEP[1])
    return XY[keep], P[keep], keep


# ---- per-npy-tile signatures (used internally for tube aggregation) ----

def per_tile_signature(npy_tiles, XY_spots, P,
                       tile_size=NPY_TILE):
    n = len(npy_tiles); m = P.shape[1]
    sig = np.zeros((n, m), dtype=np.float32)
    n_spots = np.zeros(n, dtype=int)
    tree = cKDTree(XY_spots)
    r = tile_size * np.sqrt(2)
    for i, (x, y) in enumerate(npy_tiles):
        cand = tree.query_ball_point((x + tile_size/2, y + tile_size/2), r=r)
        if not cand:
            continue
        cx = XY_spots[cand, 0]; cy = XY_spots[cand, 1]
        keep = ((cx >= x) & (cx < x + tile_size) &
                (cy >= y) & (cy < y + tile_size))
        idx = [cand[j] for j, k in enumerate(keep) if k]
        if idx:
            sig[i] = P[idx].mean(axis=0)
            n_spots[i] = len(idx)
    return sig, n_spots


def annotate_tiles(npy_tiles, roi):
    nset = {tuple(v): i for i, v in enumerate(npy_tiles.tolist())}
    tube_per_tile = [""] * len(npy_tiles)
    for tid, patches in roi.items():
        for px, py in patches:
            for dx, dy in ROI_OFFSETS:
                k = (float(px + dx), float(py + dy))
                if k in nset:
                    tube_per_tile[nset[k]] = tid
    section_per_tile = [t[0] if t else "n" for t in tube_per_tile]
    return tube_per_tile, section_per_tile


# ---- tube-level aggregation ----

def per_tube_signature(roi, npy_tiles, tile_sig, n_spots_tile, cell_cols):
    nset = {tuple(v): i for i, v in enumerate(npy_tiles.tolist())}
    rows_sig, rows_cnt, centers = [], [], []
    for tid in sort_tubes(roi.keys()):
        idx = []
        patch_xy = []
        spots_total = 0
        for px, py in roi[tid]:
            for dx, dy in ROI_OFFSETS:
                k = (float(px + dx), float(py + dy))
                if k in nset:
                    j = nset[k]
                    idx.append(j)
                    spots_total += int(n_spots_tile[j])
            patch_xy.append((px + ROI_PATCH/2, py + ROI_PATCH/2))
        sig = tile_sig[idx].mean(axis=0) if idx else np.zeros(tile_sig.shape[1], dtype=np.float32)
        rows_sig.append({"tube_id": tid, "section": tid[0],
                         "section_label": SECTION_LABEL.get(tid[0], "?"),
                         "n_patches": len(roi[tid]),
                         "n_tiles": len(idx),
                         "n_spots": int(spots_total),
                         **{c: float(sig[i]) for i, c in enumerate(cell_cols)}})
        rows_cnt.append({"tube_id": tid, "section": tid[0],
                         "section_label": SECTION_LABEL.get(tid[0], "?"),
                         "n_patches": len(roi[tid]),
                         "n_tiles": len(idx),
                         "n_spots": int(spots_total)})
        centers.append(np.mean(patch_xy, axis=0))
    return pd.DataFrame(rows_sig), pd.DataFrame(rows_cnt), np.array(centers)


def add_score_columns(df, groups, cell_cols):
    strict = groups[groups.is_strict_proxy == 1]["cell_type"].tolist()
    broad  = groups[groups.is_broad_proxy  == 1]["cell_type"].tolist()
    immune = groups[groups.group.isin(["Immune-lymphoid", "Immune-myeloid"])]["cell_type"].tolist()
    df["score_strict_proxy"] = df[strict].sum(axis=1)
    df["score_broad_proxy"]  = df[broad].sum(axis=1)
    df["score_immune_total"] = df[immune].sum(axis=1)
    return df, strict, broad, immune


# ---- statistics ----

def mw(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return {"n_a": len(a), "n_b": len(b),
                "mean_a": float(a.mean()) if len(a) else np.nan,
                "mean_b": float(b.mean()) if len(b) else np.nan,
                "delta": np.nan, "U": np.nan, "p": np.nan}
    try:
        U, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        U, p = np.nan, 1.0
    return {"n_a": int(len(a)), "n_b": int(len(b)),
            "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "delta": float(a.mean() - b.mean()),
            "U": float(U), "p": float(p)}


def section_stats(sig_df):
    sections = sig_df["section"].tolist()
    rows = []
    pairs = [
        (f"Tumor ({SECTION_LABEL['a']} vs {SECTION_LABEL['b']})", "a", "b"),
        (f"T-cell ({SECTION_LABEL['c']} vs {SECTION_LABEL['d']})", "c", "d"),
    ]
    for label, sa, sb in pairs:
        for score in ("score_strict_proxy", "score_broad_proxy",
                      "score_immune_total"):
            v = sig_df[score].values
            a = v[[i for i, s in enumerate(sections) if s == sa]]
            b = v[[i for i, s in enumerate(sections) if s == sb]]
            rows.append({"comparison": label, "score": score, **mw(a, b)})
    return pd.DataFrame(rows)


def per_celltype_wilcoxon(sig_df, cell_cols, sa="a", sb="b"):
    sections = sig_df["section"].tolist()
    P = sig_df[cell_cols].values
    a_idx = [i for i, s in enumerate(sections) if s == sa]
    b_idx = [i for i, s in enumerate(sections) if s == sb]
    rows = []
    for j, c in enumerate(cell_cols):
        a, b = P[a_idx, j], P[b_idx, j]
        try:
            U, p = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            U, p = np.nan, 1.0
        rows.append({"cell_type": c, "mean_a": float(a.mean()),
                     "mean_b": float(b.mean()),
                     "delta": float(a.mean() - b.mean()),
                     "U": float(U), "p": float(p)})
    df = pd.DataFrame(rows)
    valid = df["p"].notna()
    df.loc[valid, "p_bh"] = multipletests(df.loc[valid, "p"],
                                          method="fdr_bh")[1]
    return df.sort_values("p").reset_index(drop=True)


def proteomics_marker_check(per_ct):
    H = [
        ("KIF20A / KIF22 / INCENP (mitosis)", "Dividing_AT2",                "a>b"),
        ("KIF20A / KIF22 / INCENP (mitosis)", "Dividing_Basal",              "a>b"),
        ("KIF20A / KIF22 / INCENP (mitosis)", "Basal",                       "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",     "Muscle_smooth_syst_arterial", "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",     "Muscle_smooth_pulmonary",     "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",     "Muscle_airway",               "a>b"),
        ("(generic active Tumor)",            "AT2",                         "a>b"),
        ("(generic active Tumor)",            "Suprabasal",                  "a>b"),
    ]
    rows = []
    for prot, ctype, predicted in H:
        row = per_ct[per_ct.cell_type == ctype]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        observed = "a>b" if r.delta > 0 else "a<b"
        rows.append({"protein_marker": prot, "hist2cell_type": ctype,
                     "predicted_direction": predicted,
                     "observed_direction": observed,
                     "matches_hypothesis": predicted == observed,
                     "delta": float(r.delta),
                     "p": float(r.p),
                     "p_bh": float(r.get("p_bh", np.nan))})
    return pd.DataFrame(rows)


# ---- Moran R ----

def moran_r(P, XY, k):
    n = len(XY); kk = min(k + 1, n)
    _, nn = cKDTree(XY).query(XY, k=kk)
    nn = nn[:, 1:]
    rows = np.repeat(np.arange(n), nn.shape[1]); cols = nn.ravel()
    data = np.ones(rows.size, dtype=np.float32)
    W = csr_matrix((data, (rows, cols)), shape=(n, n))
    W = W + W.T
    W.data[:] = 1.0
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    W = W.multiply(1.0 / rs[:, None]).tocsr()
    means = P.mean(axis=0, keepdims=True)
    stds  = P.std(axis=0, keepdims=True); stds[stds == 0] = 1.0
    Z = (P - means) / stds
    return (Z.T @ (W @ Z)) / n


def moran_long(R, cell_cols):
    rows = [{"A": cell_cols[i], "B": cell_cols[j], "R": float(R[i, j])}
            for i in range(len(cell_cols))
            for j in range(i, len(cell_cols))]
    return pd.DataFrame(rows)


# ---- plot helpers ----

def _draw_mask_bg(ax, mask, W_lvl0, H_lvl0, alpha=0.3):
    ax.imshow(mask, extent=[0, W_lvl0, H_lvl0, 0],
              cmap="Greys", alpha=alpha, vmin=0, vmax=255)


def _zoom_to_mask(ax, mask, W_lvl0, H_lvl0):
    """Zoom axes to the non-zero region of the mask in level-0 px."""
    nz = np.where(mask > 0)
    if len(nz[0]) == 0:
        ax.set_xlim(0, W_lvl0); ax.set_ylim(H_lvl0, 0); return
    sy = mask.shape[0] / H_lvl0
    sx = mask.shape[1] / W_lvl0
    x0 = nz[1].min() / sx; x1 = (nz[1].max() + 1) / sx
    y0 = nz[0].min() / sy; y1 = (nz[0].max() + 1) / sy
    pad_x = (x1 - x0) * 0.02; pad_y = (y1 - y0) * 0.02
    ax.set_xlim(x0 - pad_x, x1 + pad_x)
    ax.set_ylim(y1 + pad_y, y0 - pad_y)


def _legend_sections(ax, counts):
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                          markerfacecolor=SECTION_COLOR[s], markersize=10,
                          label=f"{SECTION_LABEL[s]} ({counts.get(s, 0)})")
               for s in SECTION_ORDER]
    ax.legend(handles=handles, loc="upper right", fontsize=8,
              framealpha=0.85)


# ---- section_* plots (ROI level) ----

def plot_section_boxplots(sig_df, out_path):
    scores = [("score_strict_proxy", "strict epithelial-proliferative proxy"),
              ("score_broad_proxy",  "broad epithelial-activity proxy"),
              ("score_immune_total", "immune total")]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (col, lbl) in zip(axes, scores):
        data = [sig_df.loc[sig_df.section == s, col].values for s in SECTION_ORDER]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        labels=[SECTION_LABEL[s].replace(" ", "\n") for s in SECTION_ORDER])
        for patch, s in zip(bp["boxes"], SECTION_ORDER):
            patch.set_facecolor(SECTION_COLOR[s]); patch.set_alpha(0.55)
        for i, s in enumerate(SECTION_ORDER):
            ys = sig_df.loc[sig_df.section == s, col].values
            xs = np.random.normal(loc=i+1, scale=0.04, size=len(ys))
            ax.scatter(xs, ys, s=14, c="black", alpha=0.55, zorder=3)
        ax.set_title(lbl, fontsize=11)
        ax.set_ylabel("per-tube ROI mean")
    fig.suptitle("Per-section ROI scores — slide1 (1_085_12)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_section_scatter_panels(sig_df, tube_centers, cell_cols, mask,
                                W_lvl0, H_lvl0, columns, title, out_path,
                                ncols=5):
    """Multi-panel ROI-dot scatter over cropped tissue mask."""
    n = len(columns)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = list(axes.flat) if n > 1 else [axes]
    counts = sig_df["section"].value_counts().to_dict()
    for ax, (col, lbl) in zip(axes, columns):
        _draw_mask_bg(ax, mask, W_lvl0, H_lvl0)
        v = sig_df[col].values if col in sig_df.columns else None
        if col == "__section__":
            for s in SECTION_ORDER:
                pts = tube_centers[[i for i, sec in enumerate(sig_df.section) if sec == s]]
                if len(pts):
                    ax.scatter(pts[:, 0], pts[:, 1], s=70, c=SECTION_COLOR[s],
                               edgecolor="black", linewidth=0.3,
                               label=SECTION_LABEL[s])
        else:
            vmax = max(v.max(), 1e-6)
            sc = ax.scatter(tube_centers[:, 0], tube_centers[:, 1],
                            c=v, s=70, cmap="viridis", vmin=0, vmax=vmax,
                            edgecolor="black", linewidth=0.3)
            plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_aspect("equal")
        _zoom_to_mask(ax, mask, W_lvl0, H_lvl0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(lbl, fontsize=10)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_section_subgraph(sig_df, tube_centers, mask, W_lvl0, H_lvl0,
                          out_path):
    fig, ax = plt.subplots(figsize=(13, 7))
    _draw_mask_bg(ax, mask, W_lvl0, H_lvl0)
    XY = tube_centers
    tree = cKDTree(XY)
    _, nn = tree.query(XY, k=min(5, len(XY)))
    seen = set()
    for i in range(len(XY)):
        for j in nn[i, 1:]:
            key = tuple(sorted((i, int(j))))
            if key in seen: continue
            seen.add(key)
            ax.plot([XY[i, 0], XY[j, 0]], [XY[i, 1], XY[j, 1]],
                    c="#888", alpha=0.35, linewidth=0.6, zorder=1)
    for i, row in sig_df.iterrows():
        s = row.section
        cx, cy = tube_centers[i]
        ax.scatter(cx, cy, s=110, c=SECTION_COLOR[s],
                   edgecolor="black", linewidth=0.5, zorder=2)
        ax.annotate(row.tube_id, (cx, cy), fontsize=6, ha="center",
                    va="center", zorder=3)
    counts = sig_df["section"].value_counts().to_dict()
    _legend_sections(ax, counts)
    ax.set_aspect("equal")
    _zoom_to_mask(ax, mask, W_lvl0, H_lvl0)
    ax.set_title("ROI tube subgraph (47 nodes, kNN k=4) — coloured by section, "
                 "tissue mask backdrop",
                 fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---- spatial_* plots (Hist2Cell-spot level, with mask backdrop) ----

def _scatter_spots(ax, XY, c, *, cmap="viridis", s=1, vmin=None, vmax=None):
    sc = ax.scatter(XY[:, 0], XY[:, 1], c=c, s=s, cmap=cmap,
                    vmin=vmin, vmax=vmax)
    plt.colorbar(sc, ax=ax, fraction=0.04)
    return sc


def plot_spatial_top10(XY, P, cell_cols, ct_stats_df, mask,
                       W_lvl0, H_lvl0, out_path):
    top = ct_stats_df.head(10)
    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    for ax, (_, row) in zip(axes.flat, top.iterrows()):
        _draw_mask_bg(ax, mask, W_lvl0, H_lvl0)
        j = name_to_idx[row["cell_type"]]
        _scatter_spots(ax, XY, P[:, j], s=1, vmax=P[:, j].max())
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W_lvl0, H_lvl0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{row['cell_type']}  μ={row['mean']:.2f}", fontsize=10)
    fig.suptitle("Top-10 cell type spatial scatter — Hist2Cell spots over cropped tissue mask",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_spatial_groups(XY, P, groups_df, cell_cols, mask,
                        W_lvl0, H_lvl0, out_path):
    glist = sorted(groups_df["group"].unique())
    extras = [
        ("Strict epithelial-proliferative proxy",
         groups_df[groups_df.is_strict_proxy == 1]["cell_type"].tolist()),
        ("Broad epithelial-activity proxy",
         groups_df[groups_df.is_broad_proxy == 1]["cell_type"].tolist()),
    ]
    panels = [(g, groups_df[groups_df["group"] == g]["cell_type"].tolist())
              for g in glist] + extras
    n = len(panels); cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows))
    axes = list(axes.flat)
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    for ax, (gname, members) in zip(axes, panels):
        _draw_mask_bg(ax, mask, W_lvl0, H_lvl0)
        idx = [name_to_idx[c] for c in members]
        gsum = P[:, idx].sum(axis=1) if idx else np.zeros(len(XY))
        _scatter_spots(ax, XY, gsum, s=1)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W_lvl0, H_lvl0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{gname} (n={len(members)})  μ={gsum.mean():.2f}",
                     fontsize=10)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Spatial group + epithelial-activity proxy heatmaps — Hist2Cell spots",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_spatial_immune_vs_epithelial(XY, P, groups_df, cell_cols, mask,
                                      W_lvl0, H_lvl0, out_path):
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    immune_idx = [name_to_idx[c] for c in groups_df[groups_df.group.isin(
        ["Immune-lymphoid", "Immune-myeloid"])]["cell_type"]]
    strict_idx = [name_to_idx[c] for c in groups_df[groups_df.is_strict_proxy==1]["cell_type"]]
    broad_idx  = [name_to_idx[c] for c in groups_df[groups_df.is_broad_proxy==1]["cell_type"]]
    panels = [
        ("immune total",                    P[:, immune_idx].sum(axis=1)),
        ("strict epithelial-proliferative", P[:, strict_idx].sum(axis=1)),
        ("broad epithelial-activity",       P[:, broad_idx].sum(axis=1)),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    for ax, (lbl, v) in zip(axes, panels):
        _draw_mask_bg(ax, mask, W_lvl0, H_lvl0)
        _scatter_spots(ax, XY, v, s=1)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W_lvl0, H_lvl0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{lbl}\nμ={v.mean():.2f}  max={v.max():.2f}", fontsize=10)
    fig.suptitle("Immune vs strict / broad epithelial-activity proxy — spatial scatter",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_moran_clustermap(R, cell_cols, title, out_path):
    df = pd.DataFrame(R, index=cell_cols, columns=cell_cols)
    cg = sns.clustermap(df, cmap="vlag", center=0, vmin=-0.3, vmax=0.3,
                        figsize=(14, 14), xticklabels=True, yticklabels=True,
                        dendrogram_ratio=(0.10, 0.10),
                        cbar_pos=(0.02, 0.92, 0.05, 0.06))
    cg.ax_heatmap.tick_params(axis="x", labelsize=6, rotation=90)
    cg.ax_heatmap.tick_params(axis="y", labelsize=6, rotation=0)
    cg.fig.suptitle(title, y=1.01, fontsize=11)
    cg.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(cg.fig)


# ---- main ----

def main():
    print(f"[load] data")
    (preds, cell_cols, P_full, XY_full,
     roi, npy_full, npy_crop,
     mask_full, mask_crop, groups,
     W_lvl0, H_lvl0) = load_inputs()
    print(f"       spots={len(preds)} cell_types={len(cell_cols)}  "
          f"WSI={W_lvl0}×{H_lvl0}  mask={mask_crop.shape}")

    # Filter Hist2Cell spots to cropped mask X-range (slide-anchor)
    XY, P, keep_mask = filter_spots_to_mask_range(XY_full, P_full)
    print(f"[filter] spots within X∈{X_KEEP}: {keep_mask.sum()}/{len(keep_mask)}")

    # Per-tile signatures (using npy_full for ROI mapping coverage)
    print("[A.1] per-tile signatures over 5,227 npy tiles")
    tile_sig, n_spots_tile = per_tile_signature(npy_full, XY_full, P_full)
    tube_per_tile, section_per_tile = annotate_tiles(npy_full, roi)
    n_roi_tiles = sum(1 for t in tube_per_tile if t)
    print(f"      ROI tiles: {n_roi_tiles} ({n_roi_tiles*100/len(npy_full):.1f}%)")

    # Per-tube signatures (Scale B)
    print("[B.1] per-tube aggregation (47 tubes)")
    sig_df, cnt_df, tube_centers = per_tube_signature(
        roi, npy_full, tile_sig, n_spots_tile, cell_cols)
    sig_df, strict, broad, immune = add_score_columns(sig_df, groups, cell_cols)
    sig_df.to_csv(HERE / "roi_signatures.csv", index=False)
    cnt_df.to_csv(HERE / "roi_spot_counts.csv", index=False)
    print(cnt_df.groupby("section_label")[["n_patches", "n_tiles",
                                            "n_spots"]].agg(
        ["count", "sum", "mean"]).round(1).to_string())

    # Statistics
    print("[stat] tube-level Wilcoxon")
    sec_df = section_stats(sig_df)
    sec_df.to_csv(HERE / "section_stats.csv", index=False)
    print(sec_df.to_string(index=False))

    print("[stat] per-cell-type a vs b (tube)")
    per_ct = per_celltype_wilcoxon(sig_df, cell_cols, "a", "b")
    per_ct.to_csv(HERE / "per_celltype_wilcoxon.csv", index=False)
    print(f"      p_bh<.05: {int((per_ct.p_bh<.05).sum())}/80")
    print(per_ct.head(10)[["cell_type","mean_a","mean_b","delta","p","p_bh"]]
          .to_string(index=False))

    print("[chek] proteomics marker hypotheses")
    pm = proteomics_marker_check(per_ct)
    pm.to_csv(HERE / "proteomics_marker_hypotheses.csv", index=False)
    print(pm[["protein_marker","hist2cell_type","predicted_direction",
              "observed_direction","matches_hypothesis","delta","p","p_bh"]]
          .to_string(index=False))

    # Moran R within ROI subgraph (47 tube nodes)
    print(f"[moran] ROI subgraph (47 tubes, k={MORAN_KNN})")
    R_roi = moran_r(sig_df[cell_cols].values, tube_centers, MORAN_KNN)
    moran_long(R_roi, cell_cols).to_csv(HERE / "moran_within_roi.csv",
                                        index=False)

    # Moran R on filtered Hist2Cell spot graph (slide-wide, cropped)
    print(f"[moran] slide-wide spot graph "
          f"({len(XY)} spots, k={MORAN_KNN_SLIDE})")
    R_slide = moran_r(P, XY, MORAN_KNN_SLIDE)
    moran_long(R_slide, cell_cols).to_csv(HERE / "moran_slide_wide.csv",
                                          index=False)

    # ---- plots ----
    print("[plot] section_*")

    # section_subgraph
    plot_section_subgraph(sig_df, tube_centers, mask_crop, W_lvl0, H_lvl0,
                          HERE / "section_subgraph.png")
    plot_section_boxplots(sig_df, HERE / "section_boxplots.png")

    # section top-10 (use ROI-mean ranking)
    means_roi = sig_df[cell_cols].mean().sort_values(ascending=False)
    top_cols = [(c, f"{c}  μ_roi={means_roi[c]:.2f}") for c in means_roi.head(10).index]
    plot_section_scatter_panels(sig_df, tube_centers, cell_cols, mask_crop,
                                W_lvl0, H_lvl0, top_cols,
                                "ROI-mean top-10 cell types — per-tube scatter",
                                HERE / "section_top10_celltypes.png",
                                ncols=5)

    # section_group_heatmaps
    glist = sorted(groups["group"].unique())
    group_cols = []
    sums = {}
    for g in glist:
        members = groups[groups.group == g]["cell_type"].tolist()
        sig_df[f"_grp_{g}"] = sig_df[members].sum(axis=1)
        group_cols.append((f"_grp_{g}",
                           f"{g} (n={len(members)})  μ={sig_df[f'_grp_{g}'].mean():.2f}"))
    group_cols.append(("score_strict_proxy", "Strict epithelial-proliferative (n=3)"))
    group_cols.append(("score_broad_proxy",  "Broad epithelial-activity (n=5)"))
    plot_section_scatter_panels(sig_df, tube_centers, cell_cols, mask_crop,
                                W_lvl0, H_lvl0, group_cols,
                                "Lineage groups + proxy scores — per-tube",
                                HERE / "section_group_heatmaps.png",
                                ncols=4)

    # section_immune_vs_epithelial (3-panel)
    plot_section_scatter_panels(sig_df, tube_centers, cell_cols, mask_crop,
                                W_lvl0, H_lvl0,
                                [("score_immune_total", "immune total"),
                                 ("score_strict_proxy", "strict epithelial-proliferative"),
                                 ("score_broad_proxy",  "broad epithelial-activity")],
                                "Immune vs strict / broad epithelial-activity proxy — per-tube",
                                HERE / "section_immune_vs_epithelial.png",
                                ncols=3)

    # ---- spatial_* (Hist2Cell-spot heatmaps over cropped mask) ----
    print("[plot] spatial_*")
    ct_stats_slide = pd.DataFrame({
        "cell_type": cell_cols,
        "mean": P.mean(axis=0),
        "max":  P.max(axis=0),
    }).sort_values("mean", ascending=False).reset_index(drop=True)
    plot_spatial_top10(XY, P, cell_cols, ct_stats_slide, mask_crop,
                       W_lvl0, H_lvl0, HERE / "spatial_top10_celltypes.png")
    plot_spatial_groups(XY, P, groups, cell_cols, mask_crop, W_lvl0, H_lvl0,
                        HERE / "spatial_group_heatmaps.png")
    plot_spatial_immune_vs_epithelial(XY, P, groups, cell_cols, mask_crop,
                                      W_lvl0, H_lvl0,
                                      HERE / "spatial_immune_vs_epithelial.png")

    # Moran clustermaps
    plot_moran_clustermap(R_roi, cell_cols,
        f"Moran's R — 47-tube ROI subgraph (k={MORAN_KNN})",
        HERE / "moran_r_clustermap.png")
    plot_moran_clustermap(R_slide, cell_cols,
        f"Moran's R — Hist2Cell spot graph (cropped X-range, {len(XY)} spots, "
        f"k={MORAN_KNN_SLIDE})",
        HERE / "moran_r_clustermap_slide.png")

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
