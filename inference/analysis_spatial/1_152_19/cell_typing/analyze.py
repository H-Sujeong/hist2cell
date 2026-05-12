"""Hist2Cell ROI-level analysis for slide2 (1_152_19) — full PNG set
matching slide1's heavy version (15 figures across cell_typing + proteomics).

Section labels (slide2):
  e → High-risk Tumor
  f → Low-risk Tumor
  g → High-risk T-cell
  h → Low-risk T-cell
  v → Middle-risk Tumor (control)
  w → Middle-risk T-cell (control)  (absent in this pkl)

Outputs (same set as slide1 cell_typing/)
  roi_signatures.csv
  roi_spot_counts.csv
  section_stats.csv
  per_celltype_wilcoxon.csv
  marker_hypotheses.csv
  moran_within_roi.csv
  moran_slide_wide.csv
  section_boxplots.png
  section_subgraph.png
  section_top10_celltypes.png
  section_group_heatmaps.png
  section_immune_vs_epithelial.png
  spatial_top10_celltypes.png
  spatial_group_heatmaps.png
  spatial_immune_vs_epithelial.png
  moran_r_clustermap.png
  moran_r_clustermap_slide.png
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


HERE      = Path(__file__).resolve().parent
PARENT    = HERE.parent
PRED_CSV  = Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv")
COORDS_H5 = Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/slide2_152_19_coords.h5")
GROUPS_CSV = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")
ROI_PKL   = PARENT / "1_152_19_ROI_groups.pkl"
NPY_FULL  = PARENT / "meteo_1_152_19_coords.npy"
MASK_CROP = PARENT / "tissue_mask_cropped.png"

NPY_TILE       = 512
ROI_PATCH      = 1024
ROI_OFFSETS    = [(0, 0), (512, 0), (0, 512), (512, 512)]
MORAN_KNN      = 12
MORAN_KNN_SLIDE = 20

# slide2 X-crop: npy bounds + small buffer
X_KEEP = (65072, 153552)

SECTION_LABEL = {
    "e": "High-risk Tumor",
    "f": "Low-risk Tumor",
    "g": "High-risk T-cell",
    "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)",
    "w": "Middle-risk T-cell (ctrl)",
}
SECTION_COLOR = {
    "e":"#d62728","f":"#1f77b4","g":"#2ca02c","h":"#9467bd","v":"#7f7f7f","w":"#bcbd22",
}
SECTION_ORDER = ["e", "f", "g", "h", "v"]

HYPOTHESES = [
    ("KIF20A / KIF22 / INCENP (mitosis)",  "Dividing_AT2",                "e>f"),
    ("KIF20A / KIF22 / INCENP (mitosis)",  "Dividing_Basal",              "e>f"),
    ("KIF20A / KIF22 / INCENP (mitosis)",  "Basal",                       "e>f"),
    ("MYH11 / TAGLN (smooth muscle)",      "Muscle_smooth_syst_arterial", "e>f"),
    ("MYH11 / TAGLN (smooth muscle)",      "Muscle_smooth_pulmonary",     "e>f"),
    ("MYH11 / TAGLN (smooth muscle)",      "Muscle_airway",               "e>f"),
    ("(generic active Tumor)",             "AT2",                         "e>f"),
    ("(generic active Tumor)",             "Suprabasal",                  "e>f"),
]


def sort_tubes(keys):
    return sorted(keys, key=lambda t: (t[0],
                                       int(t[1:]) if t[1:].isdigit() else 0))


def load_inputs():
    preds = pd.read_csv(PRED_CSV)
    cell_cols = [c for c in preds.columns if c not in ("spot_id", "X", "Y")]
    P = preds[cell_cols].values.astype(np.float32)
    XY = preds[["X", "Y"]].values.astype(np.float64)
    with open(ROI_PKL, "rb") as f:
        roi = pickle.load(f)
    npy = np.load(NPY_FULL).astype(np.float64)
    mask = np.array(Image.open(MASK_CROP))
    groups = pd.read_csv(GROUPS_CSV)
    with h5py.File(COORDS_H5) as f:
        meta = dict(f["metadata"].attrs.items())
    W = int(meta["wsi_width"]); H = int(meta["wsi_height"])
    return preds, cell_cols, P, XY, roi, npy, mask, groups, W, H


def filter_spots(XY, P):
    keep = (XY[:, 0] >= X_KEEP[0]) & (XY[:, 0] <= X_KEEP[1])
    return XY[keep], P[keep]


def per_tile_signature(npy_tiles, XY, P, tile_size=NPY_TILE):
    n = len(npy_tiles); m = P.shape[1]
    sig = np.zeros((n, m), dtype=np.float32)
    n_spots = np.zeros(n, dtype=int)
    tree = cKDTree(XY)
    r = tile_size * np.sqrt(2)
    for i, (x, y) in enumerate(npy_tiles):
        cand = tree.query_ball_point((x + tile_size/2, y + tile_size/2), r=r)
        if not cand:
            continue
        cx = XY[cand, 0]; cy = XY[cand, 1]
        keep = ((cx >= x) & (cx < x + tile_size) &
                (cy >= y) & (cy < y + tile_size))
        idx = [cand[j] for j, k in enumerate(keep) if k]
        if idx:
            sig[i] = P[idx].mean(axis=0)
            n_spots[i] = len(idx)
    return sig, n_spots


def per_tube_signature(roi, npy_tiles, tile_sig, n_spots_tile, cell_cols):
    nset = {tuple(v): i for i, v in enumerate(npy_tiles.tolist())}
    rows_sig, rows_cnt, centers = [], [], []
    for tid in sort_tubes(roi.keys()):
        idx, patch_xy, spots_total = [], [], 0
        for px, py in roi[tid]:
            for dx, dy in ROI_OFFSETS:
                k = (float(px+dx), float(py+dy))
                if k in nset:
                    j = nset[k]
                    idx.append(j); spots_total += int(n_spots_tile[j])
            patch_xy.append((px + ROI_PATCH/2, py + ROI_PATCH/2))
        sig = tile_sig[idx].mean(axis=0) if idx else np.zeros(tile_sig.shape[1], dtype=np.float32)
        rec_common = {"tube_id": tid, "section": tid[0],
                      "section_label": SECTION_LABEL.get(tid[0], "?"),
                      "n_patches": len(roi[tid]),
                      "n_tiles": len(idx),
                      "n_spots": int(spots_total)}
        rows_sig.append({**rec_common, **{c: float(sig[i]) for i, c in enumerate(cell_cols)}})
        rows_cnt.append(rec_common)
        centers.append(np.mean(patch_xy, axis=0))
    return pd.DataFrame(rows_sig), pd.DataFrame(rows_cnt), np.array(centers)


def add_scores(df, groups, cell_cols):
    strict = groups[groups.is_strict_proxy == 1]["cell_type"].tolist()
    broad  = groups[groups.is_broad_proxy  == 1]["cell_type"].tolist()
    immune = groups[groups.group.isin(["Immune-lymphoid", "Immune-myeloid"])]["cell_type"].tolist()
    df["score_strict_proxy"] = df[strict].sum(axis=1)
    df["score_broad_proxy"]  = df[broad].sum(axis=1)
    df["score_immune_total"] = df[immune].sum(axis=1)
    return df


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
    sec = sig_df["section"].tolist()
    rows = []
    for label, sa, sb in [
        (f"Tumor ({SECTION_LABEL['e']} vs {SECTION_LABEL['f']})", "e", "f"),
        (f"T-cell ({SECTION_LABEL['g']} vs {SECTION_LABEL['h']})", "g", "h"),
    ]:
        for score in ("score_strict_proxy", "score_broad_proxy", "score_immune_total"):
            v = sig_df[score].values
            a = [v[i] for i, s in enumerate(sec) if s == sa]
            b = [v[i] for i, s in enumerate(sec) if s == sb]
            rows.append({"comparison": label, "score": score, **mw(a, b)})
    return pd.DataFrame(rows)


def per_celltype_wilcoxon(sig_df, cell_cols, sa="e", sb="f"):
    sec = sig_df["section"].tolist()
    P = sig_df[cell_cols].values
    a_idx = [i for i, s in enumerate(sec) if s == sa]
    b_idx = [i for i, s in enumerate(sec) if s == sb]
    rows = []
    for j, c in enumerate(cell_cols):
        a, b = P[a_idx, j], P[b_idx, j]
        try:
            U, p = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            U, p = np.nan, 1.0
        rows.append({"cell_type": c,
                     "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                     "delta": float(a.mean() - b.mean()),
                     "U": float(U), "p": float(p)})
    df = pd.DataFrame(rows)
    valid = df["p"].notna()
    df.loc[valid, "p_bh"] = multipletests(df.loc[valid, "p"],
                                          method="fdr_bh")[1]
    return df.sort_values("p").reset_index(drop=True)


def marker_check(per_ct):
    rows = []
    for prot, ctype, predicted in HYPOTHESES:
        row = per_ct[per_ct["cell_type"] == ctype]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        observed = "e>f" if r["delta"] > 0 else "e<f"
        rows.append({"protein_marker": prot, "hist2cell_type": ctype,
                     "predicted_direction": predicted,
                     "observed_direction": observed,
                     "matches_hypothesis": predicted == observed,
                     "delta": float(r["delta"]),
                     "p": float(r["p"]),
                     "p_bh": float(r.get("p_bh", np.nan))})
    return pd.DataFrame(rows)


def moran_r(P, XY, k):
    n = len(XY); kk = min(k+1, n)
    _, nn = cKDTree(XY).query(XY, k=kk)
    nn = nn[:, 1:]
    rows = np.repeat(np.arange(n), nn.shape[1])
    cols = nn.ravel()
    data = np.ones(rows.size, dtype=np.float32)
    W = csr_matrix((data, (rows, cols)), shape=(n, n))
    W = W + W.T; W.data[:] = 1.0
    rs = np.asarray(W.sum(axis=1)).ravel(); rs[rs == 0] = 1.0
    W = W.multiply(1.0/rs[:, None]).tocsr()
    means = P.mean(axis=0, keepdims=True); stds = P.std(axis=0, keepdims=True)
    stds[stds == 0] = 1.0
    Z = (P - means) / stds
    return (Z.T @ (W @ Z)) / n


def moran_long(R, cell_cols):
    rows = [{"A": cell_cols[i], "B": cell_cols[j], "R": float(R[i, j])}
            for i in range(len(cell_cols))
            for j in range(i, len(cell_cols))]
    return pd.DataFrame(rows)


# ---- plot helpers ----

def _draw_mask_bg(ax, mask, W, H, alpha=0.3):
    ax.imshow(mask, extent=[0, W, H, 0], cmap="Greys", alpha=alpha,
              vmin=0, vmax=255)


def _zoom_to_mask(ax, mask, W, H):
    nz = np.where(mask > 0)
    if len(nz[0]) == 0:
        ax.set_xlim(0, W); ax.set_ylim(H, 0); return
    sy = mask.shape[0] / H; sx = mask.shape[1] / W
    x0 = nz[1].min()/sx; x1 = (nz[1].max()+1)/sx
    y0 = nz[0].min()/sy; y1 = (nz[0].max()+1)/sy
    pad_x = (x1-x0)*0.02; pad_y = (y1-y0)*0.02
    ax.set_xlim(x0-pad_x, x1+pad_x); ax.set_ylim(y1+pad_y, y0-pad_y)


def plot_section_boxplots(sig_df, out):
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
    fig.suptitle("Per-section ROI scores — slide2 (1_152_19)", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_section_subgraph(sig_df, tube_centers, mask, W, H, out):
    fig, ax = plt.subplots(figsize=(13, 7))
    _draw_mask_bg(ax, mask, W, H)
    XY = tube_centers
    tree = cKDTree(XY); _, nn = tree.query(XY, k=min(5, len(XY)))
    seen = set()
    for i in range(len(XY)):
        for j in nn[i, 1:]:
            key = tuple(sorted((i, int(j))))
            if key in seen: continue
            seen.add(key)
            ax.plot([XY[i,0], XY[j,0]], [XY[i,1], XY[j,1]],
                    c="#888", alpha=0.35, linewidth=0.6, zorder=1)
    for i, row in sig_df.iterrows():
        s = row.section; cx, cy = tube_centers[i]
        ax.scatter(cx, cy, s=110, c=SECTION_COLOR[s],
                   edgecolor="black", linewidth=0.5, zorder=2)
        ax.annotate(row.tube_id, (cx, cy), fontsize=6, ha="center", va="center", zorder=3)
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                           markerfacecolor=SECTION_COLOR[s], markersize=10,
                           label=f"{SECTION_LABEL[s]} ({sig_df[sig_df.section==s].shape[0]})")
               for s in SECTION_ORDER]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.85)
    ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W, H)
    ax.set_title("ROI tube subgraph (48 nodes, kNN k=4) — slide2 — tissue mask backdrop",
                 fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def _scatter_roi(ax, tube_centers, vals, title, s=70):
    sc = ax.scatter(tube_centers[:, 0], tube_centers[:, 1],
                    c=vals, s=s, cmap="viridis", edgecolor="black", linewidth=0.3)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.04)


def plot_section_panels(sig_df, tube_centers, mask, W, H, columns, title, out, ncols=5):
    n = len(columns); nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5*ncols, 4.2*nrows))
    axes = list(axes.flat) if n > 1 else [axes]
    for ax, (col, lbl) in zip(axes, columns):
        _draw_mask_bg(ax, mask, W, H)
        v = sig_df[col].values
        vmax = max(v.max(), 1e-6)
        sc = ax.scatter(tube_centers[:, 0], tube_centers[:, 1],
                        c=v, s=70, cmap="viridis", vmin=0, vmax=vmax,
                        edgecolor="black", linewidth=0.3)
        plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W, H)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(lbl, fontsize=10)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


def plot_spatial_top10(XY, P, cell_cols, ct_stats, mask, W, H, out):
    top = ct_stats.head(10)
    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    for ax, (_, row) in zip(axes.flat, top.iterrows()):
        _draw_mask_bg(ax, mask, W, H)
        j = name_to_idx[row["cell_type"]]
        sc = ax.scatter(XY[:, 0], XY[:, 1], c=P[:, j], s=1, cmap="viridis",
                        vmin=0, vmax=P[:, j].max())
        plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W, H)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{row['cell_type']}  μ={row['mean']:.2f}", fontsize=10)
    fig.suptitle("Top-10 cell type spatial scatter — slide2 Hist2Cell spots", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


def plot_spatial_groups(XY, P, groups, cell_cols, mask, W, H, out):
    glist = sorted(groups["group"].unique())
    extras = [("Strict epithelial-proliferative proxy",
               groups[groups.is_strict_proxy == 1]["cell_type"].tolist()),
              ("Broad epithelial-activity proxy",
               groups[groups.is_broad_proxy == 1]["cell_type"].tolist())]
    panels = [(g, groups[groups.group == g]["cell_type"].tolist()) for g in glist] + extras
    n = len(panels); cols = 4; rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5*cols, 4.5*rows))
    axes = list(axes.flat)
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    for ax, (gname, members) in zip(axes, panels):
        _draw_mask_bg(ax, mask, W, H)
        idx = [name_to_idx[c] for c in members]
        gsum = P[:, idx].sum(axis=1) if idx else np.zeros(len(XY))
        sc = ax.scatter(XY[:, 0], XY[:, 1], c=gsum, s=1, cmap="viridis")
        plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W, H)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{gname} (n={len(members)})  μ={gsum.mean():.2f}", fontsize=10)
    for ax in axes[n:]: ax.axis("off")
    fig.suptitle("Spatial group + proxy heatmaps — slide2 Hist2Cell spots", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


def plot_spatial_immune_vs_epithelial(XY, P, groups, cell_cols, mask, W, H, out):
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    immune_idx = [name_to_idx[c] for c in groups[groups.group.isin(
        ["Immune-lymphoid","Immune-myeloid"])]["cell_type"]]
    strict_idx = [name_to_idx[c] for c in groups[groups.is_strict_proxy==1]["cell_type"]]
    broad_idx  = [name_to_idx[c] for c in groups[groups.is_broad_proxy==1]["cell_type"]]
    panels = [
        ("immune total", P[:, immune_idx].sum(axis=1)),
        ("strict epithelial-proliferative", P[:, strict_idx].sum(axis=1)),
        ("broad epithelial-activity", P[:, broad_idx].sum(axis=1)),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    for ax, (lbl, v) in zip(axes, panels):
        _draw_mask_bg(ax, mask, W, H)
        sc = ax.scatter(XY[:, 0], XY[:, 1], c=v, s=1, cmap="viridis")
        plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_aspect("equal"); _zoom_to_mask(ax, mask, W, H)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{lbl}\nμ={v.mean():.2f}, max={v.max():.2f}", fontsize=10)
    fig.suptitle("Immune vs strict / broad epithelial-activity proxy — slide2", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


def plot_moran_clustermap(R, cell_cols, title, out):
    df = pd.DataFrame(R, index=cell_cols, columns=cell_cols)
    cg = sns.clustermap(df, cmap="vlag", center=0, vmin=-0.3, vmax=0.3,
                        figsize=(14, 14), xticklabels=True, yticklabels=True,
                        dendrogram_ratio=(0.10, 0.10),
                        cbar_pos=(0.02, 0.92, 0.05, 0.06))
    cg.ax_heatmap.tick_params(axis="x", labelsize=6, rotation=90)
    cg.ax_heatmap.tick_params(axis="y", labelsize=6, rotation=0)
    cg.fig.suptitle(title, y=1.01, fontsize=11)
    cg.savefig(out, dpi=120, bbox_inches="tight"); plt.close(cg.fig)


def main():
    print("[load]")
    preds, cell_cols, P_full, XY_full, roi, npy, mask, groups, W_lvl0, H_lvl0 = load_inputs()
    print(f"  spots={len(preds)}, npy_tiles={len(npy)}, tubes={len(roi)}")

    XY, P = filter_spots(XY_full, P_full)
    print(f"[filter] spots in X∈{X_KEEP}: {len(XY)}/{len(XY_full)}")

    print("[A] per-tile signatures")
    tile_sig, n_spots_tile = per_tile_signature(npy, XY_full, P_full)

    print("[B] per-tube signatures (48 tubes)")
    sig_df, cnt_df, tube_centers = per_tube_signature(
        roi, npy, tile_sig, n_spots_tile, cell_cols)
    sig_df = add_scores(sig_df, groups, cell_cols)
    sig_df.to_csv(HERE/"roi_signatures.csv", index=False)
    cnt_df.to_csv(HERE/"roi_spot_counts.csv", index=False)

    print("[stat] section + per-cell-type")
    sec_df = section_stats(sig_df)
    sec_df.to_csv(HERE/"section_stats.csv", index=False)
    per_ct = per_celltype_wilcoxon(sig_df, cell_cols, "e", "f")
    per_ct.to_csv(HERE/"per_celltype_wilcoxon.csv", index=False)
    mk = marker_check(per_ct); mk.to_csv(HERE/"marker_hypotheses.csv", index=False)
    print(f"  per_ct p_bh<.05: {int((per_ct.p_bh<.05).sum())}/80")

    print(f"[moran] ROI subgraph + slide-wide")
    R_roi = moran_r(sig_df[cell_cols].values, tube_centers, MORAN_KNN)
    moran_long(R_roi, cell_cols).to_csv(HERE/"moran_within_roi.csv", index=False)
    R_slide = moran_r(P, XY, MORAN_KNN_SLIDE)
    moran_long(R_slide, cell_cols).to_csv(HERE/"moran_slide_wide.csv", index=False)

    print("[plot] section_*")
    plot_section_subgraph(sig_df, tube_centers, mask, W_lvl0, H_lvl0,
                          HERE/"section_subgraph.png")
    plot_section_boxplots(sig_df, HERE/"section_boxplots.png")

    means_roi = sig_df[cell_cols].mean().sort_values(ascending=False)
    top_cols = [(c, f"{c}  μ_roi={means_roi[c]:.2f}") for c in means_roi.head(10).index]
    plot_section_panels(sig_df, tube_centers, mask, W_lvl0, H_lvl0, top_cols,
                        "Top-10 ROI-mean cell types — slide2 per-tube",
                        HERE/"section_top10_celltypes.png", ncols=5)

    glist = sorted(groups["group"].unique())
    group_cols = []
    for g in glist:
        members = groups[groups.group == g]["cell_type"].tolist()
        sig_df[f"_grp_{g}"] = sig_df[members].sum(axis=1)
        group_cols.append((f"_grp_{g}",
                           f"{g} (n={len(members)})  μ={sig_df[f'_grp_{g}'].mean():.2f}"))
    group_cols.append(("score_strict_proxy", "Strict epithelial-proliferative (n=3)"))
    group_cols.append(("score_broad_proxy",  "Broad epithelial-activity (n=5)"))
    plot_section_panels(sig_df, tube_centers, mask, W_lvl0, H_lvl0, group_cols,
                        "Lineage groups + proxy — slide2 per-tube",
                        HERE/"section_group_heatmaps.png", ncols=4)

    plot_section_panels(sig_df, tube_centers, mask, W_lvl0, H_lvl0,
                        [("score_immune_total", "immune total"),
                         ("score_strict_proxy", "strict epithelial-proliferative"),
                         ("score_broad_proxy",  "broad epithelial-activity")],
                        "Immune vs strict / broad proxy — slide2 per-tube",
                        HERE/"section_immune_vs_epithelial.png", ncols=3)

    print("[plot] spatial_*")
    ct_stats_slide = pd.DataFrame({
        "cell_type": cell_cols,
        "mean": P.mean(axis=0),
        "max": P.max(axis=0),
    }).sort_values("mean", ascending=False).reset_index(drop=True)
    plot_spatial_top10(XY, P, cell_cols, ct_stats_slide, mask, W_lvl0, H_lvl0,
                       HERE/"spatial_top10_celltypes.png")
    plot_spatial_groups(XY, P, groups, cell_cols, mask, W_lvl0, H_lvl0,
                       HERE/"spatial_group_heatmaps.png")
    plot_spatial_immune_vs_epithelial(XY, P, groups, cell_cols, mask, W_lvl0, H_lvl0,
                                      HERE/"spatial_immune_vs_epithelial.png")

    plot_moran_clustermap(R_roi, cell_cols,
        f"Moran's R — 48-tube ROI subgraph (slide2, k={MORAN_KNN})",
        HERE/"moran_r_clustermap.png")
    plot_moran_clustermap(R_slide, cell_cols,
        f"Moran's R — slide2 Hist2Cell spot graph "
        f"(cropped X-range, {len(XY)} spots, k={MORAN_KNN_SLIDE})",
        HERE/"moran_r_clustermap_slide.png")

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png", ".py"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
