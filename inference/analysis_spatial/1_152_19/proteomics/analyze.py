"""Minimal proteomics differential analysis for slide2 (1_152_19) —
only outputs needed by ../proofs/core_proofs.py (Claim 1).

Slide2 columns in gg_matrix: e1-e13 (13) + f1-f15 (15) + g1-g7 (7) +
h1-h8 (8) + v1-v5 (5) = 48 samples.  (The duplicate 'e2' has already
been removed at the gg_matrix level — see commit a03c5d0.)

Pipeline
  1. Filter to slide2 columns, log2-transform (NaN preserved).
  2. Quality filter — drop genes that aren't detected in ≥30% of either
     comparison group (e for the Tumor pair).
  3. Per-gene Mann-Whitney U for e vs f (Tumor) and g vs h (T-cell),
     each comparison BH-FDR corrected.
  4. Marker hypothesis check on the e vs f Tumor table — same
     pre-registered set as slide1 so the two slides are directly
     comparable.

Outputs (this folder)
  tumor_e_vs_f_genes.csv          BH-sorted per-gene table
  tcell_g_vs_h_genes.csv          same for T-cell pair
  marker_hypothesis_check.csv     MYH11 / TAGLN / KIF20A / KIF22 / INCENP
                                  direction check (Tumor e vs f)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


HERE = Path(__file__).resolve().parent
GG   = Path("/home/sjhong/hist2cell/inference/analysis_spatial/report.gg_matrix (1).tsv")

SECTION_LABEL = {
    "e": "High-risk Tumor", "f": "Low-risk Tumor",
    "g": "High-risk T-cell", "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)", "w": "Middle-risk T-cell (ctrl)",
}

# Mirror slide1 pre-registered set so the cross-slide comparison is
# apples-to-apples; predicted direction for the high vs low Tumor pair
# is e>f (slide2's e is High-risk Tumor).
MARKER_HYPOTHESES = [
    ("KIF20A",  "e>f"),
    ("KIF22",   "e>f"),
    ("INCENP",  "e>f"),
    ("MYH11",   "e>f"),
    ("TAGLN",   "e>f"),
    ("NCAM1",   "e>f"),
    ("APOBEC3C","e>f"),
]


def quality_filter(log2_df, cols, sec_a, sec_b, min_detect=0.30):
    a_cols = [c for c in cols if c[0] == sec_a]
    b_cols = [c for c in cols if c[0] == sec_b]
    fa = log2_df[a_cols].notna().mean(axis=1)
    fb = log2_df[b_cols].notna().mean(axis=1)
    return (fa >= min_detect) & (fb >= min_detect)


def per_gene_mw(log2_df, slide2_cols, sec_a, sec_b):
    a_cols = [c for c in slide2_cols if c[0] == sec_a]
    b_cols = [c for c in slide2_cols if c[0] == sec_b]
    rows = []
    for g in log2_df.index:
        ra = log2_df.loc[g, a_cols].dropna().values
        rb = log2_df.loc[g, b_cols].dropna().values
        if len(ra) < 3 or len(rb) < 3:
            rows.append({"gene": g, "n_a": len(ra), "n_b": len(rb),
                         "log2_mean_a": float(ra.mean()) if len(ra) else np.nan,
                         "log2_mean_b": float(rb.mean()) if len(rb) else np.nan,
                         "log2_fc": np.nan, "U": np.nan, "p": np.nan})
            continue
        try:
            U, p = mannwhitneyu(ra, rb, alternative="two-sided")
        except ValueError:
            U, p = np.nan, 1.0
        rows.append({"gene": g, "n_a": int(len(ra)), "n_b": int(len(rb)),
                     "log2_mean_a": float(ra.mean()),
                     "log2_mean_b": float(rb.mean()),
                     "log2_fc": float(ra.mean() - rb.mean()),
                     "U": float(U), "p": float(p)})
    df = pd.DataFrame(rows)
    valid = df["p"].notna()
    df.loc[valid, "p_bh"] = multipletests(df.loc[valid, "p"],
                                          method="fdr_bh")[1]
    return df.sort_values("p").reset_index(drop=True)


def marker_check(per_gene_df, hypotheses):
    rows = []
    for gene, predicted in hypotheses:
        row = per_gene_df[per_gene_df["gene"] == gene]
        if len(row) == 0:
            rows.append({"gene": gene, "predicted_direction": predicted,
                         "observed_direction": "(not measured / filtered)",
                         "matches_hypothesis": False,
                         "log2_fc": np.nan, "p": np.nan, "p_bh": np.nan})
            continue
        r = row.iloc[0]
        observed = "e>f" if r["log2_fc"] > 0 else "e<f"
        rows.append({"gene": gene, "predicted_direction": predicted,
                     "observed_direction": observed,
                     "matches_hypothesis": predicted == observed,
                     "log2_fc": float(r["log2_fc"]),
                     "p": float(r["p"]),
                     "p_bh": float(r.get("p_bh", np.nan))})
    return pd.DataFrame(rows)


def main():
    print("[load] gg_matrix slide2 columns")
    m = pd.read_csv(GG, sep="\t")
    slide2_cols = [c for c in m.columns[3:] if c[0] in "efghv"]
    print(f"  total genes: {len(m)}, slide2 samples: {len(slide2_cols)}")
    sec_counts = {s: sum(1 for c in slide2_cols if c[0] == s)
                  for s in "efghv"}
    print(f"  per-section counts: {sec_counts}")

    raw = m[slide2_cols].replace(0, np.nan)
    log2 = np.log2(raw); log2.index = m["Genes"].astype(str)
    print(f"  log2 range: [{np.nanmin(log2.values):.2f}, {np.nanmax(log2.values):.2f}]")

    # ── Tumor e vs f ──
    print("\n[1] Tumor e vs f (High-risk vs Low-risk)")
    keep = quality_filter(log2, slide2_cols, "e", "f")
    print(f"  pass detect-≥30%: {int(keep.sum())}/{len(log2)} "
          f"(n_e=13, n_f=15 samples)")
    per_ef = per_gene_mw(log2[keep], slide2_cols, "e", "f")
    per_ef.to_csv(HERE / "tumor_e_vs_f_genes.csv", index=False)
    print(f"  BH-FDR<0.05: {int((per_ef['p_bh']<0.05).sum())}")
    print(per_ef.head(10)[["gene", "n_a", "n_b", "log2_fc", "p", "p_bh"]]
          .to_string(index=False))

    # ── T-cell g vs h ──
    print("\n[2] T-cell g vs h")
    keep_gh = quality_filter(log2, slide2_cols, "g", "h")
    per_gh = per_gene_mw(log2[keep_gh], slide2_cols, "g", "h")
    per_gh.to_csv(HERE / "tcell_g_vs_h_genes.csv", index=False)
    print(f"  pass: {int(keep_gh.sum())} / BH<0.05: {int((per_gh['p_bh']<0.05).sum())}")

    # ── Marker hypothesis check ──
    print("\n[3] pre-registered marker hypothesis check (Tumor e vs f)")
    mk = marker_check(per_ef, MARKER_HYPOTHESES)
    mk.to_csv(HERE / "marker_hypothesis_check.csv", index=False)
    print(mk.to_string(index=False))

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".py"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
