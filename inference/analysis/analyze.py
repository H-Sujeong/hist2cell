"""Per-slide downstream analysis on Hist2Cell predictions.

Inputs (already produced by `prep/prepare_wsi_for_inference_v2.py` +
`inference/infer.py`):
    - predictions.csv : spot_id, X, Y, <80 cell-type abundance columns>
    - <slide>_coords.h5 : level-0 tile coordinates and metadata
    - cell_type_groups.csv : lineage grouping of the 80 lung cell types
                             + `is_strict_proxy` / `is_broad_proxy` flags
                             (see EPITHELIAL_PROXY_METHODOLOGY.md)

Outputs (written to --output):
    - abundance_by_celltype.csv : per-cell-type mean / median / max / fraction-nonzero
    - abundance_by_group.csv    : per-group total mean / spot-fraction; also two
                                  pseudo-groups: "Strict epithelial-proliferative
                                  proxy" (3 types) and "Broad epithelial-activity
                                  proxy" (5 types)
    - spatial_top10_celltypes.png : 10-panel scatter of top-mean cell types
    - spatial_group_heatmaps.png  : one panel per group, abundance summed over the group
    - spatial_immune_vs_epithelial.png : 1×3 panel — immune total / strict proxy / broad proxy
    - moran_r_pairs.csv         : bivariate Moran's R for every cell-pair (3160 rows)
    - moran_r_clustermap.png    : hierarchical clustermap of the 80×80 R matrix

Usage:
    python inference/analysis/analyze.py \\
        --predictions inference/slide1_085_12_v2/predictions.csv \\
        --coords      inference/slide1_085_12_v2/slide1_085_12_coords.h5 \\
        --groups      inference/analysis/cell_type_groups.csv \\
        --output      inference/analysis/slide1_085_12_v2

Caveat: the model was trained on healthy human lung; cell-type column names
are lung-specific. The two epithelial-activity proxies are NOT tumor
detectors — they are lung-derived spatial proxies whose breast-context
validity is itself a hypothesis to be tested by an external breast-trained
model (CUCA her2st). See EPITHELIAL_PROXY_METHODOLOGY.md and README.md.
"""

import argparse
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
from scipy.stats import norm


# ---------------- I/O ----------------

def load_predictions(csv_path: Path):
    df = pd.read_csv(csv_path)
    meta_cols = ["spot_id", "X", "Y"]
    cell_cols = [c for c in df.columns if c not in meta_cols]
    return df[meta_cols].copy(), df[cell_cols].values.astype(np.float32), cell_cols


def load_groups(csv_path: Path, expected_celltypes):
    g = pd.read_csv(csv_path)
    miss = set(expected_celltypes) - set(g["cell_type"])
    extra = set(g["cell_type"]) - set(expected_celltypes)
    if miss:
        raise SystemExit(f"groups CSV is missing {len(miss)} cell types: {sorted(miss)[:5]}...")
    if extra:
        raise SystemExit(f"groups CSV has {len(extra)} unknown cell types: {sorted(extra)[:5]}...")
    g = g.set_index("cell_type").loc[expected_celltypes].reset_index()
    return g


def load_coords_h5(h5_path: Path):
    """Returns dict of metadata + coords array (top-left tile xy)."""
    with h5py.File(h5_path, "r") as f:
        coords = f["coords"][:]
        meta = {k: (v.decode() if isinstance(v, bytes) else v)
                for k, v in f["metadata"].attrs.items()}
    return coords, meta


# ---------------- per-celltype / per-group stats ----------------

def per_celltype_stats(preds, cell_cols):
    rows = []
    for j, name in enumerate(cell_cols):
        col = preds[:, j]
        rows.append(dict(
            cell_type=name,
            mean=float(col.mean()),
            median=float(np.median(col)),
            max=float(col.max()),
            fraction_nonzero=float((col > 0).mean()),
        ))
    out = pd.DataFrame(rows).sort_values("mean", ascending=False).reset_index(drop=True)
    return out


def per_group_stats(preds, groups_df, cell_cols):
    """Sum predictions within each lineage group; also produce two pseudo-groups
    for the strict (3-type) and broad (5-type) epithelial-activity proxies.
    See EPITHELIAL_PROXY_METHODOLOGY.md for the rationale of the two scores."""
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    rows = []
    grouping = list(groups_df.groupby("group"))
    grouping.append(("Strict epithelial-proliferative proxy",
                     groups_df[groups_df["is_strict_proxy"] == 1]))
    grouping.append(("Broad epithelial-activity proxy",
                     groups_df[groups_df["is_broad_proxy"] == 1]))
    for gname, g in grouping:
        idx = [name_to_idx[c] for c in g["cell_type"].tolist()]
        sub = preds[:, idx]
        spot_sum = sub.sum(axis=1)
        rows.append(dict(
            group=gname,
            n_celltypes=len(idx),
            mean_per_spot=float(spot_sum.mean()),
            max_per_spot=float(spot_sum.max()),
            fraction_spots_nonzero=float((spot_sum > 0).mean()),
            sum_total=float(spot_sum.sum()),
        ))
    out = pd.DataFrame(rows).sort_values("mean_per_spot", ascending=False).reset_index(drop=True)
    return out


# ---------------- spatial plots ----------------

def _scatter(ax, X, Y, c, title, vmin=None, vmax=None):
    sc = ax.scatter(X, Y, c=c, s=1, cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(sc, ax=ax, fraction=0.04)


def plot_top10(meta_df, preds, cell_cols, ct_stats_df, out_path):
    top = ct_stats_df.head(10)
    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    for ax, (_, row) in zip(axes.flat, top.iterrows()):
        j = name_to_idx[row["cell_type"]]
        _scatter(ax, meta_df["X"], meta_df["Y"], preds[:, j],
                 f"{row['cell_type']}  μ={row['mean']:.2f}")
    plt.tight_layout()
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_groups(meta_df, preds, groups_df, cell_cols, out_path):
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    glist = sorted(groups_df["group"].unique())
    n = len(glist)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows))
    axes = axes.flat if rows > 1 else axes
    for ax, gname in zip(axes, glist):
        members = groups_df[groups_df["group"] == gname]["cell_type"].tolist()
        idx = [name_to_idx[c] for c in members]
        gsum = preds[:, idx].sum(axis=1)
        _scatter(ax, meta_df["X"], meta_df["Y"], gsum,
                 f"{gname} (n={len(members)})  μ={gsum.mean():.2f}")
    for ax in list(axes)[n:]:
        ax.axis("off")
    plt.tight_layout()
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_immune_vs_epithelial(meta_df, preds, groups_df, cell_cols, out_path):
    """3-panel: immune total / strict epithelial-proliferative / broad epithelial-activity.

    The two epithelial scores are NOT tumor detectors — they are lung-derived
    proxies (strict = Basal/Dividing_AT2/Dividing_Basal; broad = +AT2/Suprabasal).
    See EPITHELIAL_PROXY_METHODOLOGY.md."""
    name_to_idx = {n: i for i, n in enumerate(cell_cols)}
    immune_members = groups_df[groups_df["group"].isin(["Immune-lymphoid", "Immune-myeloid"])]["cell_type"].tolist()
    strict_members = groups_df[groups_df["is_strict_proxy"] == 1]["cell_type"].tolist()
    broad_members  = groups_df[groups_df["is_broad_proxy"]  == 1]["cell_type"].tolist()
    immune_idx = [name_to_idx[c] for c in immune_members]
    strict_idx = [name_to_idx[c] for c in strict_members]
    broad_idx  = [name_to_idx[c] for c in broad_members]
    immune_sum = preds[:, immune_idx].sum(axis=1)
    strict_sum = preds[:, strict_idx].sum(axis=1)
    broad_sum  = preds[:, broad_idx].sum(axis=1)
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))
    _scatter(axes[0], meta_df["X"], meta_df["Y"], immune_sum,
             f"Immune total (n={len(immune_idx)})  μ={immune_sum.mean():.2f}")
    _scatter(axes[1], meta_df["X"], meta_df["Y"], strict_sum,
             f"Strict epithelial-proliferative proxy (n={len(strict_idx)})  μ={strict_sum.mean():.2f}")
    _scatter(axes[2], meta_df["X"], meta_df["Y"], broad_sum,
             f"Broad epithelial-activity proxy (n={len(broad_idx)})  μ={broad_sum.mean():.2f}")
    fig.suptitle("Immune vs strict / broad epithelial-activity proxy\n"
                 "(lung-trained model → see EPITHELIAL_PROXY_METHODOLOGY.md)", fontsize=11)
    plt.tight_layout()
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


# ---------------- Moran's R ----------------

def build_knn_weight_matrix(coords_xy, k=20):
    """Sparse row-normalized kNN weight matrix W (csr). Self excluded."""
    n = len(coords_xy)
    tree = cKDTree(coords_xy)
    _, nn = tree.query(coords_xy, k=min(k + 1, n))   # includes self at column 0
    nn = nn[:, 1:]
    rows = np.repeat(np.arange(n), nn.shape[1])
    cols = nn.ravel()
    data = np.ones(rows.size, dtype=np.float32)
    W = csr_matrix((data, (rows, cols)), shape=(n, n))
    # symmetrize (union)
    W = (W + W.T)
    W.data[:] = 1.0
    # row-normalize
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    W = W.multiply(1.0 / row_sum[:, None]).tocsr()
    return W


def moran_r_pairs(preds, W):
    """Bivariate Moran's R for every cell-type pair (j<k).

    For row-standardized W and z-scored variables x, y:
        R(x,y) = z_x^T W z_y / n
    Analytic z under randomization is approximated using the diagonal
    of W^T W (Cliff & Ord) — for inference scale we treat each pair
    as independent.
    """
    n, m = preds.shape
    # z-score
    means = preds.mean(axis=0, keepdims=True)
    stds = preds.std(axis=0, keepdims=True)
    stds[stds == 0] = 1.0
    Z = (preds - means) / stds                    # (n, m)
    WZ = W @ Z                                     # (n, m), row-normalized lag

    # R = (Z^T (W Z)) / n  — bivariate Moran's I/R
    R = (Z.T @ WZ) / n                              # (m, m)

    # variance approximation under randomization: V0 = (n-1)/(n+1) ... too involved;
    # fall back to permutation-free heuristic SE = sqrt(2 * (W**2).sum()) / n
    W_sq = W.multiply(W).sum()
    se = float(np.sqrt(2.0 * W_sq) / n) if W_sq > 0 else 1.0
    Z_score = R / se
    P_val = 2.0 * norm.sf(np.abs(Z_score))

    return R, Z_score, P_val


def write_moran_pairs_csv(R, Z, P, cell_cols, out_path):
    m = R.shape[0]
    rows = []
    for i in range(m):
        for j in range(i, m):
            rows.append(dict(
                A=cell_cols[i], B=cell_cols[j],
                R=float(R[i, j]),
                z=float(Z[i, j]),
                p=float(P[i, j]),
            ))
    pd.DataFrame(rows).to_csv(out_path, index=False)


def plot_moran_clustermap(R, cell_cols, out_path):
    df = pd.DataFrame(R, index=cell_cols, columns=cell_cols)
    cg = sns.clustermap(
        df, cmap="vlag", center=0, vmin=-0.3, vmax=0.3,
        figsize=(14, 14), xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.10, 0.10), cbar_pos=(0.02, 0.92, 0.05, 0.06),
    )
    cg.ax_heatmap.tick_params(axis="x", labelsize=6, rotation=90)
    cg.ax_heatmap.tick_params(axis="y", labelsize=6, rotation=0)
    cg.fig.suptitle("Bivariate Moran's R (kNN-weighted) on Hist2Cell predictions",
                    y=1.01, fontsize=12)
    cg.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(cg.fig)


# ---------------- main ----------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", required=True, type=Path)
    ap.add_argument("--coords", required=True, type=Path)
    ap.add_argument("--groups", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--knn", type=int, default=20,
                    help="k for the weight matrix used by Moran's R (default 20)")
    return ap.parse_args()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    print(f"[load] predictions  = {args.predictions}")
    meta_df, preds, cell_cols = load_predictions(args.predictions)
    print(f"       n_spots = {len(meta_df)}, n_celltypes = {len(cell_cols)}")
    print(f"[load] coords h5    = {args.coords}")
    coords_h5, h5_meta = load_coords_h5(args.coords)
    if len(coords_h5) != len(meta_df):
        print(f"  WARN: coords h5 has {len(coords_h5)} entries, predictions have {len(meta_df)}")
    print(f"[load] groups       = {args.groups}")
    groups_df = load_groups(args.groups, cell_cols)
    print(f"       groups: {sorted(groups_df['group'].unique())}")
    print(f"       strict proxy (3): {groups_df[groups_df['is_strict_proxy']==1]['cell_type'].tolist()}")
    print(f"       broad  proxy (5): {groups_df[groups_df['is_broad_proxy']==1]['cell_type'].tolist()}")

    # --- per-celltype stats ---
    print("[stats] per-cell-type")
    ct_stats = per_celltype_stats(preds, cell_cols)
    ct_stats.to_csv(args.output / "abundance_by_celltype.csv", index=False)

    # --- per-group stats ---
    print("[stats] per-group")
    g_stats = per_group_stats(preds, groups_df, cell_cols)
    g_stats.to_csv(args.output / "abundance_by_group.csv", index=False)

    # --- spatial plots ---
    print("[plot] top-10 cell types")
    plot_top10(meta_df, preds, cell_cols, ct_stats, args.output / "spatial_top10_celltypes.png")
    print("[plot] groups")
    plot_groups(meta_df, preds, groups_df, cell_cols, args.output / "spatial_group_heatmaps.png")
    print("[plot] immune vs epithelial (strict + broad)")
    plot_immune_vs_epithelial(meta_df, preds, groups_df, cell_cols,
                              args.output / "spatial_immune_vs_epithelial.png")

    # --- Moran's R ---
    print(f"[moran] kNN(k={args.knn}) weight matrix")
    coords_xy = meta_df[["X", "Y"]].values.astype(np.float64)
    W = build_knn_weight_matrix(coords_xy, k=args.knn)
    print(f"        weight nnz = {W.nnz}")
    print("[moran] bivariate Moran's R for 80×80 cell-type pairs")
    R, Z, P = moran_r_pairs(preds, W)
    print(f"        R range [{R.min():.3f}, {R.max():.3f}], diag mean = {np.diag(R).mean():.3f}")
    write_moran_pairs_csv(R, Z, P, cell_cols, args.output / "moran_r_pairs.csv")
    print("[plot] Moran's R clustermap")
    plot_moran_clustermap(R, cell_cols, args.output / "moran_r_clustermap.png")

    print(f"\nDone. Outputs in {args.output}")
    for p in sorted(args.output.iterdir()):
        if p.is_file():
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
