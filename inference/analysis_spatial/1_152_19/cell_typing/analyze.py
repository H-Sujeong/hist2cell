"""Minimal Hist2Cell ROI-level analysis for slide2 (1_152_19) —
only outputs needed by ../proofs/core_proofs.py (Claim 1 + Claim 2).

Section labels (slide2):
  e → High-risk Tumor
  f → Low-risk Tumor
  g → High-risk T-cell
  h → Low-risk T-cell
  v → Middle-risk Tumor (control)
  w → Middle-risk T-cell (control)  (absent in this pkl)

Inputs
  ../1_152_19_ROI_groups.pkl              48 tubes, 182 patches @ 1024×1024
  ../meteo_1_152_19_coords.npy            8,668 candidate 512-tile top-lefts
  /home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv
                                          40,502 spots × 80 cell types
  ../../../analysis/cell_type_groups.csv  strict / broad / lineage flags

Outputs (this folder)
  roi_signatures.csv               per-tube 48 × (80 cell types + 3 scores)
  roi_spot_counts.csv              per-tube n_patches / n_tiles / n_spots
  section_stats.csv                e vs f and g vs h Wilcoxon (3 scores)
  per_celltype_wilcoxon.csv        80-row e vs f Wilcoxon + BH-FDR
  marker_hypotheses.csv            pre-registered marker-celltype direction check
"""
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


HERE      = Path(__file__).resolve().parent
PARENT    = HERE.parent
PRED_CSV  = Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv")
ROI_PKL   = PARENT / "1_152_19_ROI_groups.pkl"
NPY_FULL  = PARENT / "meteo_1_152_19_coords.npy"
GROUPS_CSV = Path("/home/sjhong/hist2cell/inference/analysis/cell_type_groups.csv")

NPY_TILE     = 512
ROI_PATCH    = 1024
ROI_OFFSETS  = [(0, 0), (512, 0), (0, 512), (512, 512)]

SECTION_LABEL = {
    "e": "High-risk Tumor",
    "f": "Low-risk Tumor",
    "g": "High-risk T-cell",
    "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)",
    "w": "Middle-risk T-cell (ctrl)",
}

# Pre-registered marker-celltype hypotheses (re-use slide1 set; both
# slides are TNBC-context KBSMC samples and the lung-derived proxy
# convention applies the same way).
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


def per_tile_signature(npy_tiles, XY_spots, P, tile_size=NPY_TILE):
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


def per_tube_signature(roi, npy_tiles, tile_sig, n_spots_tile, cell_cols):
    nset = {tuple(v): i for i, v in enumerate(npy_tiles.tolist())}
    rows_sig, rows_cnt = [], []
    for tid in sort_tubes(roi.keys()):
        idx = []
        spots_total = 0
        for px, py in roi[tid]:
            for dx, dy in ROI_OFFSETS:
                k = (float(px + dx), float(py + dy))
                if k in nset:
                    j = nset[k]
                    idx.append(j)
                    spots_total += int(n_spots_tile[j])
        sig = (tile_sig[idx].mean(axis=0) if idx else
               np.zeros(tile_sig.shape[1], dtype=np.float32))
        rec_common = {"tube_id": tid, "section": tid[0],
                      "section_label": SECTION_LABEL.get(tid[0], "?"),
                      "n_patches": len(roi[tid]),
                      "n_tiles": len(idx),
                      "n_spots": int(spots_total)}
        rows_sig.append({**rec_common,
                         **{c: float(sig[i]) for i, c in enumerate(cell_cols)}})
        rows_cnt.append(rec_common)
    return pd.DataFrame(rows_sig), pd.DataFrame(rows_cnt)


def add_scores(sig_df, groups, cell_cols):
    strict = groups[groups.is_strict_proxy == 1]["cell_type"].tolist()
    broad  = groups[groups.is_broad_proxy  == 1]["cell_type"].tolist()
    immune = groups[groups.group.isin(["Immune-lymphoid", "Immune-myeloid"])]["cell_type"].tolist()
    sig_df["score_strict_proxy"] = sig_df[strict].sum(axis=1)
    sig_df["score_broad_proxy"]  = sig_df[broad].sum(axis=1)
    sig_df["score_immune_total"] = sig_df[immune].sum(axis=1)
    return sig_df


def mw(values, sections, sa, sb):
    a = np.asarray([v for v, s in zip(values, sections) if s == sa])
    b = np.asarray([v for v, s in zip(values, sections) if s == sb])
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
            rows.append({"comparison": label, "score": score,
                         **mw(sig_df[score].tolist(), sec, sa, sb)})
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
                     "mean_a": float(a.mean()),
                     "mean_b": float(b.mean()),
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


def main():
    print("[load] preds + ROI pkl + npy + groups")
    preds = pd.read_csv(PRED_CSV)
    cell_cols = [c for c in preds.columns if c not in ("spot_id", "X", "Y")]
    P = preds[cell_cols].values.astype(np.float32)
    XY = preds[["X", "Y"]].values.astype(np.float64)
    with open(ROI_PKL, "rb") as f:
        roi = pickle.load(f)
    npy = np.load(NPY_FULL).astype(np.float64)
    groups = pd.read_csv(GROUPS_CSV)
    print(f"  spots={len(preds)} cell_types={len(cell_cols)} "
          f"tubes={len(roi)} patches={sum(len(v) for v in roi.values())} "
          f"npy_tiles={len(npy)}")

    print("[A] per-npy-tile signatures (8,668 tiles)")
    tile_sig, n_spots_tile = per_tile_signature(npy, XY, P)
    print(f"  tiles with ≥1 spot inside: {int((n_spots_tile>0).sum())}/{len(n_spots_tile)}")

    print("[B] per-tube signatures (48 tubes)")
    sig_df, cnt_df = per_tube_signature(roi, npy, tile_sig, n_spots_tile, cell_cols)
    sig_df = add_scores(sig_df, groups, cell_cols)
    sig_df.to_csv(HERE / "roi_signatures.csv", index=False)
    cnt_df.to_csv(HERE / "roi_spot_counts.csv", index=False)
    print(cnt_df.groupby("section_label")[["n_patches","n_tiles","n_spots"]].agg(
        ["count","sum","mean"]).round(1).to_string())

    print("[C] section_stats — Wilcoxon e vs f and g vs h")
    sec_df = section_stats(sig_df)
    sec_df.to_csv(HERE / "section_stats.csv", index=False)
    print(sec_df.to_string(index=False))

    print("[D] per-cell-type e vs f (80 type) + BH-FDR")
    per_ct = per_celltype_wilcoxon(sig_df, cell_cols, "e", "f")
    per_ct.to_csv(HERE / "per_celltype_wilcoxon.csv", index=False)
    print(f"  p_bh<0.05: {int((per_ct['p_bh']<0.05).sum())}/80")
    print(per_ct.head(10)[["cell_type","mean_a","mean_b","delta","p","p_bh"]]
          .to_string(index=False))

    print("[E] pre-registered marker hypothesis check (Hist2Cell side)")
    mk = marker_check(per_ct)
    mk.to_csv(HERE / "marker_hypotheses.csv", index=False)
    print(mk[["protein_marker","hist2cell_type","predicted_direction",
              "observed_direction","matches_hypothesis","delta","p","p_bh"]]
          .to_string(index=False))

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".py"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
