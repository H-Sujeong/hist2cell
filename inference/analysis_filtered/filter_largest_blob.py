"""Crop a v2 slide's spots to the X-range of its largest spatial blob.

Algorithm:
  1. Load predictions.csv + <slide>_coords.h5 from inference/<slide>_v2/.
  2. Build a symmetrized kNN graph (k=6) over the spot (X, Y) positions.
  3. Find connected components; pick the one with the most spots.
  4. Take that component's [Xmin, Xmax].
  5. Filter ALL spots — keep iff their X is in [Xmin, Xmax]. (Y is not
     constrained: per request, only the X-range of the largest blob.)
  6. Write filtered predictions.csv + <slide>_coords.h5 + spots.csv into
     inference/analysis_filtered/<slide>_v2/.

Usage:
    python inference/analysis_filtered/filter_largest_blob.py \\
        --in-dir   inference/slide1_085_12_v2 \\
        --out-dir  inference/analysis_filtered/slide1_085_12_v2 \\
        --slide    slide1_085_12
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--slide", required=True, help="slide name (used for *_coords.h5 filename)")
    ap.add_argument("--knn", type=int, default=6, help="k for connected-component graph")
    return ap.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    preds_csv = args.in_dir / "predictions.csv"
    coords_h5 = args.in_dir / f"{args.slide}_coords.h5"
    spots_csv = args.in_dir / "spots.csv"

    print(f"[load] {preds_csv}")
    preds = pd.read_csv(preds_csv)
    n0 = len(preds)
    XY = preds[["X", "Y"]].values.astype(np.float64)
    print(f"       n_spots = {n0}")

    print(f"[graph] kNN(k={args.knn}) symmetrized for connected components")
    k = min(args.knn + 1, n0)
    _, nn = cKDTree(XY).query(XY, k=k)
    nn = nn[:, 1:]  # drop self
    rows = np.repeat(np.arange(n0), nn.shape[1])
    cols = nn.ravel()
    data = np.ones(rows.size, dtype=np.uint8)
    G = csr_matrix((data, (rows, cols)), shape=(n0, n0))
    G = G + G.T  # symmetrize (union)

    n_comp, labels = connected_components(G, directed=False)
    sizes = np.bincount(labels)
    largest = int(np.argmax(sizes))
    print(f"        n_components = {n_comp}")
    top5 = np.argsort(sizes)[::-1][:5]
    for r, c in enumerate(top5):
        print(f"        comp[{int(c)}] size = {int(sizes[c])} ({100.0*sizes[c]/n0:.2f}%)"
              + ("  <-- largest" if c == largest else ""))

    blob_X = XY[labels == largest, 0]
    xmin, xmax = float(blob_X.min()), float(blob_X.max())
    print(f"[xrange] largest blob X: [{xmin:.0f}, {xmax:.0f}]  (level-0 px)")

    keep = (XY[:, 0] >= xmin) & (XY[:, 0] <= xmax)
    n1 = int(keep.sum())
    print(f"[keep ] {n1}/{n0} spots ({100.0*n1/n0:.2f}%)")

    # filtered predictions
    filt_idx = np.where(keep)[0]
    filt_preds = preds.iloc[filt_idx].reset_index(drop=True)
    filt_preds.to_csv(args.out_dir / "predictions.csv", index=False)

    # filtered coords.h5 (preserve metadata, add filter info)
    with h5py.File(coords_h5, "r") as fin:
        coords_arr = fin["coords"][:]
        meta = {k_: v_ for k_, v_ in fin["metadata"].attrs.items()}
    if len(coords_arr) != n0:
        # coords.h5 should have the same row count as predictions
        raise SystemExit(f"coords.h5 has {len(coords_arr)} rows, predictions {n0}")
    filt_coords = coords_arr[keep]
    out_h5 = args.out_dir / f"{args.slide}_coords.h5"
    with h5py.File(out_h5, "w") as fout:
        fout.create_dataset("coords", data=filt_coords, compression="gzip")
        m = fout.create_group("metadata")
        for k_, v_ in meta.items():
            m.attrs[k_] = v_
        m.attrs["total_tiles"] = int(n1)
        m.attrs["filter_kind"] = "largest_blob_xrange"
        m.attrs["filter_xmin"] = xmin
        m.attrs["filter_xmax"] = xmax
        m.attrs["filter_kept"] = int(n1)
        m.attrs["filter_orig_total"] = int(n0)

    # filtered spots.csv (if present)
    if spots_csv.exists():
        spots = pd.read_csv(spots_csv)
        if len(spots) == n0:
            spots.iloc[filt_idx].reset_index(drop=True).to_csv(
                args.out_dir / "spots.csv", index=False)

    # also save a tiny summary CSV
    summary = pd.DataFrame([{
        "slide": args.slide,
        "n_orig": n0,
        "n_filtered": n1,
        "kept_pct": round(100.0 * n1 / n0, 3),
        "n_components": int(n_comp),
        "largest_comp_size": int(sizes[largest]),
        "largest_comp_pct": round(100.0 * sizes[largest] / n0, 3),
        "xmin": xmin,
        "xmax": xmax,
    }])
    summary.to_csv(args.out_dir / "filter_summary.csv", index=False)

    print(f"\nDone:")
    print(f"  {args.out_dir/'predictions.csv'}")
    print(f"  {out_h5}")
    if spots_csv.exists():
        print(f"  {args.out_dir/'spots.csv'}")
    print(f"  {args.out_dir/'filter_summary.csv'}")


if __name__ == "__main__":
    main()
