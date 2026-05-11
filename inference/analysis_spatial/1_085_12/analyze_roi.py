"""ROI-level quantitative analysis for slide1 (1_085_12) — subgraph at
the proteomics 270 μm patch scale, anchored on the user-provided
coordinate set.

Inputs
  ./1_085_12_ROI_groups.pkl              dict { tube_id : list[(x, y), ...] }
                                          47 tubes / 181 patches at level-0 px
  ./meteo_1_085_12_coords.npy            (5227, 2) tilemap of all candidate
                                          patch top-lefts at 512-px grid;
                                          context only — ROI coords are a
                                          subset of this set (181/181 verified)
  /home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv
                                          35,821 spots × 80 cell types
                                          (spot center X, Y at level-0 px,
                                           tile_size 400)
  /home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv
                                          strict / broad proxy flags

Geometry
  Each ROI patch coordinate is the TOP-LEFT of a 1024×1024 level-0 area
  (= 270 μm at 0.2615 mpp = proteomics extraction patch). A Hist2Cell
  spot (X, Y) is inside the patch iff X ∈ [px, px+1024) and Y ∈ [py, py+1024).

  At 400-px Hist2Cell tile size, each 1024×1024 ROI patch contains
  ~(1024/400)² ≈ 6.5 spots on average. Per-tube (3-6 patches) gives
  ~20-40 spots → stable aggregation for ROI-level statistics.

Analysis
  1. ROI patch → Hist2Cell spot mapping by bbox containment.
  2. Per-tube signature: mean cell-type abundance over contained spots.
  3. Section comparisons (Mann-Whitney U):
     - a vs b  : high-risk vs low-risk Tumor
     - c vs d  : high-risk vs low-risk T-cell
     Three pre-registered scores: strict / broad epithelial-activity
     proxy + immune total.
  4. Per-cell-type a vs b Wilcoxon for all 80 types, BH-FDR correction.
  5. Proteomics marker hypothesis checks (slide1 high-risk Tumor
     markers from existing findings: KIF20A/22/INCENP ↔ Dividing_*/Basal,
     MYH11/TAGLN ↔ Stromal-muscle).
  6. Within-ROI subgraph Moran R: kNN (k=12) on the 181 patch CENTERS
     using cKDTree, bivariate Moran R for the 80×80 cell-type matrix.
  7. PNG: ROI placement overlay coloured by section,
          per-section boxplots of the three scores,
          per-tube subgraph (patches as nodes coloured by section,
                             edges within tubes).

Outputs (this folder)
  roi_spot_counts.csv             per-tube n_patches + n_spots aggregated
  roi_signatures.csv              per-tube cell-type mean + 3 scores
  section_stats.csv               three scores' a vs b + c vs d Wilcoxon
  per_celltype_wilcoxon.csv       80 cell type a vs b with BH-FDR
  proteomics_marker_matches.csv   pre-registered marker hypothesis checks
  moran_within_roi.csv            80×80 Moran R on ROI patch subgraph
  spatial_roi_overlay.png         slide-level scatter + ROI rectangles
  section_boxplots.png            strict / broad / immune per section
  roi_subgraph.png                47 tube subgraph layout
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
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


HERE         = Path(__file__).resolve().parent
PRED_CSV     = Path("/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv")
ROI_PKL      = HERE / "1_085_12_ROI_groups.pkl"
CAND_NPY     = HERE / "meteo_1_085_12_coords.npy"
GROUPS_CSV   = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")

PATCH_SIZE_PX = 1024     # 270 μm @ 0.2615 mpp; level-0 px
MORAN_KNN     = 12       # ROI patch subgraph is smaller than slide-wide


# ---- I/O ----

def load_inputs():
    preds = pd.read_csv(PRED_CSV)
    cell_cols = [c for c in preds.columns if c not in ("spot_id", "X", "Y")]
    P = preds[cell_cols].values.astype(np.float32)
    XY = preds[["X", "Y"]].values.astype(np.float64)
    with open(ROI_PKL, "rb") as f:
        roi = pickle.load(f)
    candidate = np.load(CAND_NPY)
    groups = pd.read_csv(GROUPS_CSV)
    return preds, cell_cols, P, XY, roi, candidate, groups


# ---- ROI → spot mapping ----

def map_patches_to_spots(roi: dict, XY: np.ndarray):
    """Return list of (tube_id, patch_idx_in_tube, [spot_idx, ...])."""
    out = []
    for tid in sort_tubes(roi.keys()):
        for k, (px, py) in enumerate(roi[tid]):
            mask = ((XY[:, 0] >= px) & (XY[:, 0] < px + PATCH_SIZE_PX) &
                    (XY[:, 1] >= py) & (XY[:, 1] < py + PATCH_SIZE_PX))
            out.append((tid, k, np.nonzero(mask)[0].tolist(),
                        (int(px), int(py))))
    return out


def sort_tubes(keys):
    return sorted(keys, key=lambda t: (t[0],
                                       int(t[1:]) if t[1:].isdigit() else 0))


# ---- per-tube signature ----

def per_tube_signature(roi: dict, XY: np.ndarray, P: np.ndarray, cell_cols):
    rows_sig = []
    rows_cnt = []
    tube_centers = []
    for tid in sort_tubes(roi.keys()):
        spots_union = set()
        patch_xy = []
        for px, py in roi[tid]:
            mask = ((XY[:, 0] >= px) & (XY[:, 0] < px + PATCH_SIZE_PX) &
                    (XY[:, 1] >= py) & (XY[:, 1] < py + PATCH_SIZE_PX))
            spots_union.update(np.nonzero(mask)[0].tolist())
            patch_xy.append((px + PATCH_SIZE_PX / 2, py + PATCH_SIZE_PX / 2))
        idx = sorted(spots_union)
        if idx:
            sig = P[idx].mean(axis=0)
        else:
            sig = np.zeros(P.shape[1], dtype=np.float32)
        rows_sig.append({"tube_id": tid, "section": tid[0],
                         "n_patches": len(roi[tid]), "n_spots": len(idx),
                         **{c: float(sig[i]) for i, c in enumerate(cell_cols)}})
        rows_cnt.append({"tube_id": tid, "section": tid[0],
                         "n_patches": len(roi[tid]), "n_spots": len(idx)})
        tube_centers.append(np.mean(patch_xy, axis=0))
    sig_df = pd.DataFrame(rows_sig)
    cnt_df = pd.DataFrame(rows_cnt)
    return sig_df, cnt_df, np.array(tube_centers)


def add_score_columns(sig_df, groups, cell_cols):
    strict = groups[groups.is_strict_proxy == 1]["cell_type"].tolist()
    broad  = groups[groups.is_broad_proxy  == 1]["cell_type"].tolist()
    immune = groups[groups.group.isin(
        ["Immune-lymphoid", "Immune-myeloid"])]["cell_type"].tolist()
    sig_df["score_strict_proxy"] = sig_df[strict].sum(axis=1)
    sig_df["score_broad_proxy"]  = sig_df[broad].sum(axis=1)
    sig_df["score_immune_total"] = sig_df[immune].sum(axis=1)
    return sig_df, strict, broad, immune


# ---- statistics ----

def mw(values, sections, sec_a, sec_b):
    a = values[[i for i, s in enumerate(sections) if s == sec_a]]
    b = values[[i for i, s in enumerate(sections) if s == sec_b]]
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
    for label, sa, sb in [("Tumor a vs b (high vs low risk)", "a", "b"),
                          ("T-cell c vs d (high vs low risk)", "c", "d")]:
        for score in ("score_strict_proxy", "score_broad_proxy",
                      "score_immune_total"):
            r = mw(sig_df[score].values, sections, sa, sb)
            rows.append({"comparison": label, "score": score, **r})
    return pd.DataFrame(rows)


def per_celltype_wilcoxon(sig_df, cell_cols, sec_a="a", sec_b="b"):
    sections = sig_df["section"].tolist()
    P = sig_df[cell_cols].values
    a_idx = [i for i, s in enumerate(sections) if s == sec_a]
    b_idx = [i for i, s in enumerate(sections) if s == sec_b]
    rows = []
    for j, ctype in enumerate(cell_cols):
        a = P[a_idx, j]
        b = P[b_idx, j]
        try:
            U, p = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            U, p = np.nan, 1.0
        rows.append({"cell_type": ctype,
                     "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                     "delta": float(a.mean() - b.mean()),
                     "U": float(U), "p": float(p)})
    df = pd.DataFrame(rows)
    valid = df["p"].notna()
    df.loc[valid, "p_bh"] = multipletests(df.loc[valid, "p"],
                                          method="fdr_bh")[1]
    return df.sort_values("p").reset_index(drop=True)


def proteomics_marker_check(per_ct):
    """Pre-registered (from existing findings.md):
       high-risk Tumor (a) > low-risk (b) for these protein/cell pairs."""
    H = [
        ("KIF20A / KIF22 / INCENP (mitosis)",  "Dividing_AT2",                  "a>b"),
        ("KIF20A / KIF22 / INCENP (mitosis)",  "Dividing_Basal",                "a>b"),
        ("KIF20A / KIF22 / INCENP (mitosis)",  "Basal",                         "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",      "Muscle_smooth_syst_arterial",   "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",      "Muscle_smooth_pulmonary",       "a>b"),
        ("MYH11 / TAGLN (smooth muscle)",      "Muscle_airway",                 "a>b"),
        ("(generic active Tumor)",             "AT2",                           "a>b"),
        ("(generic active Tumor)",             "Suprabasal",                    "a>b"),
    ]
    rows = []
    for prot, ctype, predicted in H:
        row = per_ct[per_ct["cell_type"] == ctype]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        observed = "a>b" if r["delta"] > 0 else "a<b"
        match = (predicted == observed)
        rows.append({"protein_marker": prot, "hist2cell_type": ctype,
                     "predicted_direction": predicted,
                     "observed_direction": observed,
                     "matches_hypothesis": bool(match),
                     "delta": float(r["delta"]),
                     "p": float(r["p"]),
                     "p_bh": float(r.get("p_bh", np.nan))})
    return pd.DataFrame(rows)


# ---- Moran R on ROI subgraph ----

def moran_r_subgraph(sig_df, cell_cols, tube_centers, k=MORAN_KNN):
    """Moran R built on the 47 (or 181) tube-level nodes using their
    spatial centers + kNN graph. We use TUBE centers (one node per tube,
    47 nodes) since the proteomics signal is per-tube, not per-patch."""
    P = sig_df[cell_cols].values.astype(np.float32)
    XY = tube_centers
    n = len(XY)
    kk = min(k + 1, n)
    _, nn = cKDTree(XY).query(XY, k=kk)
    nn = nn[:, 1:]
    rows = np.repeat(np.arange(n), nn.shape[1])
    cols = nn.ravel()
    data = np.ones(rows.size, dtype=np.float32)
    W = csr_matrix((data, (rows, cols)), shape=(n, n))
    W = W + W.T
    W.data[:] = 1.0
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    W = W.multiply(1.0 / rs[:, None]).tocsr()
    means = P.mean(axis=0, keepdims=True)
    stds  = P.std(axis=0, keepdims=True)
    stds[stds == 0] = 1.0
    Z = (P - means) / stds
    WZ = W @ Z
    R = (Z.T @ WZ) / n
    return R


# ---- plots ----

SECTION_COLOR = {"a": "#d62728", "b": "#1f77b4", "c": "#2ca02c",
                 "d": "#9467bd", "t": "#7f7f7f"}

def plot_overlay(XY_all, roi: dict, sig_df, out_path):
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.scatter(XY_all[:, 0], XY_all[:, 1], s=0.4, c="#dddddd", alpha=0.5)
    for s in ["b", "d", "c", "a", "t"]:
        for tid in [t for t in sort_tubes(roi.keys()) if t[0] == s]:
            for px, py in roi[tid]:
                rect = plt.Rectangle((px, py), PATCH_SIZE_PX, PATCH_SIZE_PX,
                                     linewidth=0.6, edgecolor=SECTION_COLOR[s],
                                     facecolor=SECTION_COLOR[s], alpha=0.35)
                ax.add_patch(rect)
    legend = [plt.Line2D([0], [0], marker="s", color="w",
                         markerfacecolor=SECTION_COLOR[s], markersize=10,
                         label=f"{s} ({sig_df[sig_df.section==s].shape[0]} tubes)")
              for s in ["a", "b", "c", "d", "t"]]
    ax.legend(handles=legend, loc="best", frameon=True)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title("ROI proteomics patches over Hist2Cell spots — slide1 (1_085_12)\n"
                 "(grey = Hist2Cell spots, colored squares = 1024×1024 px ROI patches)")
    ax.set_xlabel("X (level-0 px)")
    ax.set_ylabel("Y (level-0 px)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_boxplots(sig_df, out_path):
    order = ["a", "b", "c", "d", "t"]
    scores = [("score_strict_proxy", "strict epithelial-proliferative proxy"),
              ("score_broad_proxy",  "broad epithelial-activity proxy"),
              ("score_immune_total", "immune total")]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (col, lbl) in zip(axes, scores):
        data = [sig_df.loc[sig_df.section == s, col].values for s in order]
        bp = ax.boxplot(data, labels=order, patch_artist=True, widths=0.6)
        for patch, s in zip(bp["boxes"], order):
            patch.set_facecolor(SECTION_COLOR[s]); patch.set_alpha(0.55)
        # add raw points
        for i, s in enumerate(order):
            ys = sig_df.loc[sig_df.section == s, col].values
            xs = np.random.normal(loc=i+1, scale=0.04, size=len(ys))
            ax.scatter(xs, ys, s=10, c="black", alpha=0.55, zorder=3)
        ax.set_title(lbl, fontsize=11)
        ax.set_xlabel("section")
        ax.set_ylabel("per-ROI mean")
    fig.suptitle("ROI scores per section — slide1 (a/b: Tumor high/low, "
                 "c/d: T-cell high/low, t: Tumor ctrl)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _scatter_roi(ax, tube_centers, values, title, *, cmap="viridis", s=80):
    sc = ax.scatter(tube_centers[:, 0], tube_centers[:, 1],
                    c=values, s=s, cmap=cmap, edgecolor="black", linewidth=0.4)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.04)


def plot_top10_roi(sig_df, cell_cols, tube_centers, out_path):
    """10-panel scatter of 47 ROI tube centers coloured by each top-10
    cell type's per-tube mean."""
    means = sig_df[cell_cols].mean().sort_values(ascending=False)
    top10 = means.head(10)
    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    for ax, (ctype, mval) in zip(axes.flat, top10.items()):
        _scatter_roi(ax, tube_centers, sig_df[ctype].values,
                     f"{ctype}  μ_roi={mval:.2f}")
    fig.suptitle("Top-10 cell types by ROI mean — abundance across 47 ROI tubes "
                 "(slide1, ROI-level aggregation)", fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_groups_roi(sig_df, groups, cell_cols, tube_centers, out_path):
    """One panel per lineage group + 2 score panels (strict / broad)."""
    glist = sorted(groups["group"].unique())
    extras = [
        ("Strict epithelial-proliferative proxy",
         groups[groups.is_strict_proxy == 1]["cell_type"].tolist()),
        ("Broad epithelial-activity proxy",
         groups[groups.is_broad_proxy == 1]["cell_type"].tolist()),
    ]
    panels = [(g, groups[groups.group == g]["cell_type"].tolist())
              for g in glist] + extras
    n = len(panels)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows))
    axes = list(axes.flat)
    for ax, (gname, members) in zip(axes, panels):
        vals = sig_df[members].sum(axis=1).values if members else \
            np.zeros(len(sig_df))
        _scatter_roi(ax, tube_centers, vals,
                     f"{gname} (n={len(members)})  μ={vals.mean():.2f}")
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Lineage groups + epithelial-activity proxies at ROI scale "
                 "— slide1 (47 tubes)", fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_immune_vs_epithelial_roi(sig_df, tube_centers, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    panels = [
        ("score_immune_total", "immune total (36 types)"),
        ("score_strict_proxy", "strict epithelial-proliferative (3 types)"),
        ("score_broad_proxy",  "broad epithelial-activity (5 types)"),
    ]
    for ax, (col, lbl) in zip(axes, panels):
        v = sig_df[col].values
        _scatter_roi(ax, tube_centers, v,
                     f"{lbl}\nμ_roi={v.mean():.2f}, max={v.max():.2f}",
                     s=110)
    fig.suptitle("Per-ROI scores: immune vs strict vs broad epithelial-activity proxy "
                 "— slide1 (47 tubes, lung-derived proxy)", fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_moran_clustermap_roi(R, cell_cols, out_path):
    df = pd.DataFrame(R, index=cell_cols, columns=cell_cols)
    cg = sns.clustermap(df, cmap="vlag", center=0, vmin=-0.3, vmax=0.3,
                        figsize=(14, 14), xticklabels=True, yticklabels=True,
                        dendrogram_ratio=(0.10, 0.10),
                        cbar_pos=(0.02, 0.92, 0.05, 0.06))
    cg.ax_heatmap.tick_params(axis="x", labelsize=6, rotation=90)
    cg.ax_heatmap.tick_params(axis="y", labelsize=6, rotation=0)
    cg.fig.suptitle("Bivariate Moran's R — within-ROI subgraph (47 tube nodes, "
                    "kNN k=12) — slide1", y=1.01, fontsize=11)
    cg.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(cg.fig)


def plot_subgraph(sig_df, tube_centers, out_path):
    fig, ax = plt.subplots(figsize=(11, 7))
    for i, row in sig_df.iterrows():
        s = row.section
        cx, cy = tube_centers[i]
        ax.scatter(cx, cy, s=120, c=SECTION_COLOR[s], edgecolor="black",
                   linewidth=0.5, zorder=2)
        ax.annotate(row.tube_id, (cx, cy), fontsize=7, ha="center",
                    va="center", zorder=3)
    # kNN edges within k=4 for visual
    XY = tube_centers
    tree = cKDTree(XY)
    _, nn = tree.query(XY, k=min(5, len(XY)))
    seen = set()
    for i in range(len(XY)):
        for j in nn[i, 1:]:
            key = tuple(sorted((i, int(j))))
            if key in seen:
                continue
            seen.add(key)
            ax.plot([XY[i, 0], XY[j, 0]], [XY[i, 1], XY[j, 1]],
                    c="#888888", alpha=0.4, linewidth=0.5, zorder=1)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title("ROI tube subgraph (kNN k=4) — 47 nodes coloured by section, "
                 "edges = spatial nearest neighbours among tubes",
                 fontsize=11)
    ax.set_xlabel("X (level-0 px)")
    ax.set_ylabel("Y (level-0 px)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---- main ----

def main():
    print(f"[load] {PRED_CSV}")
    preds, cell_cols, P, XY, roi, candidate, groups = load_inputs()
    print(f"       n_spots={len(preds)}  n_celltypes={len(cell_cols)}")
    print(f"       ROI dict: {len(roi)} tubes, "
          f"{sum(len(v) for v in roi.values())} patches")

    print("[map ] ROI patches → Hist2Cell spots (1024 px bbox)")
    sig_df, cnt_df, tube_centers = per_tube_signature(roi, XY, P, cell_cols)
    cnt_df.to_csv(HERE / "roi_spot_counts.csv", index=False)
    print(cnt_df.groupby("section")[["n_patches","n_spots"]].agg(
        ["count","sum","mean"]).round(1).to_string())

    print("[sig ] per-tube signatures + score columns")
    sig_df, strict, broad, immune = add_score_columns(sig_df, groups, cell_cols)
    sig_df.to_csv(HERE / "roi_signatures.csv", index=False)

    print("[wilc] section comparisons")
    sec_df = section_stats(sig_df)
    sec_df.to_csv(HERE / "section_stats.csv", index=False)
    print(sec_df.to_string(index=False))

    print("[wilc] per-cell-type a vs b")
    per_ct = per_celltype_wilcoxon(sig_df, cell_cols, "a", "b")
    per_ct.to_csv(HERE / "per_celltype_wilcoxon.csv", index=False)
    sig_at_p05 = int((per_ct["p_bh"] < 0.05).sum())
    print(f"       p_bh<0.05: {sig_at_p05}/80")
    print(per_ct.head(10)[["cell_type","mean_a","mean_b","delta","p","p_bh"]]
          .to_string(index=False))

    print("[chek] pre-registered proteomics marker hypotheses")
    pm = proteomics_marker_check(per_ct)
    pm.to_csv(HERE / "proteomics_marker_matches.csv", index=False)
    print(pm.to_string(index=False))

    print(f"[moran] within-ROI subgraph (47 tube nodes, k={MORAN_KNN})")
    R = moran_r_subgraph(sig_df, cell_cols, tube_centers, MORAN_KNN)
    rows = []
    for i in range(len(cell_cols)):
        for j in range(i, len(cell_cols)):
            rows.append({"A": cell_cols[i], "B": cell_cols[j],
                         "R": float(R[i, j])})
    moran_df = pd.DataFrame(rows)
    moran_df.to_csv(HERE / "moran_within_roi.csv", index=False)
    off = moran_df[moran_df.A != moran_df.B]
    print(f"       diag mean={moran_df[moran_df.A==moran_df.B]['R'].mean():.3f}")
    print("       top 5 positive off-diag:")
    print(off.nlargest(5, "R").to_string(index=False))
    print("       top 5 negative off-diag:")
    print(off.nsmallest(5, "R").to_string(index=False))

    print("[plot] overlay / boxplots / subgraph / top10 / groups / "
          "immune-vs-epithelial / Moran clustermap")
    plot_overlay(XY, roi, sig_df, HERE / "spatial_roi_overlay.png")
    plot_boxplots(sig_df, HERE / "section_boxplots.png")
    plot_subgraph(sig_df, tube_centers, HERE / "roi_subgraph.png")
    plot_top10_roi(sig_df, cell_cols, tube_centers,
                   HERE / "spatial_top10_celltypes.png")
    plot_groups_roi(sig_df, groups, cell_cols, tube_centers,
                    HERE / "spatial_group_heatmaps.png")
    plot_immune_vs_epithelial_roi(sig_df, tube_centers,
                                  HERE / "spatial_immune_vs_epithelial.png")
    plot_moran_clustermap_roi(R, cell_cols,
                              HERE / "moran_r_clustermap.png")

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
