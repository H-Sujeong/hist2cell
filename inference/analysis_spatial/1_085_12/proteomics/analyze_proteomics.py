"""Proteomics differential analysis on slide1 (1_085_12) ROI tubes
using the gg_matrix from `analysis_spatial/report.gg_matrix (1).tsv`.

Section labels (per user spec):
  a → High-risk Tumor
  b → Low-risk Tumor
  c → High-risk T-cell
  d → Low-risk T-cell
  t → Middle-risk Tumor (control)

Slide1 columns in the gg_matrix: a2-a10 (9 expected, 8 present — a5 missing)
+ b1-b21 (21) + c1-c5 (5) + d1-d9 (9) + t1-t3 (3) = 46 samples.

Pipeline
  1. Filter to slide1 columns, log2-transform (NaN preserved).
  2. Quality filter — drop genes that are not measurable in at least
     30% of either comparison group ('a' or 'b' for tumor pair, 'c' or
     'd' for T-cell pair).
  3. Per-gene Mann-Whitney U test for two pre-registered comparisons:
       - a (High-risk Tumor)   vs b (Low-risk Tumor)
       - c (High-risk T-cell)  vs d (Low-risk T-cell)
     Each comparison BH-FDR corrected independently.
  4. Pre-registered marker hypothesis check (KIF20A / KIF22 / INCENP
     for mitosis; MYH11 / TAGLN for smooth muscle — same set used by
     the cell_typing analysis so the two modalities can be cross-
     compared at protein and cell-type levels.)
  5. PCA of samples coloured by section + descriptive label.
  6. Volcano plots, sample×top-marker heatmap.

Inputs
  ../../report.gg_matrix (1).tsv         (TSV, ~7800 genes × 95 samples)
  ../cell_typing/roi_signatures.csv      (47 ROI tubes × 80 cell types
                                           + 3 scores — cross-modality
                                           overlay)

Outputs (this folder)
  slide1_columns_summary.csv      sample → section + n_detected + median
  log2_intensity_matrix.csv       46 sample × n_gene log2 matrix (only
                                  genes passing the quality filter)
  tumor_a_vs_b_genes.csv          ranked per-gene Wilcoxon table (BH)
  tcell_c_vs_d_genes.csv          same for c vs d
  marker_hypothesis_check.csv     pre-registered markers' direction + p
  pca_samples.csv                 sample PCs + section label
  volcano_tumor_a_vs_b.png        log2FC vs -log10 p
  volcano_tcell_c_vs_d.png
  pca_samples.png                 PC1 vs PC2 coloured by section
  top_markers_heatmap.png         sample × top-20 a-vs-b markers
  section_protein_summary.png     per-section mean # detected genes +
                                  median intensity
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


HERE       = Path(__file__).resolve().parent
GG_MATRIX  = Path("/home/sjhong/hist2cell/inference/analysis_spatial/report.gg_matrix (1).tsv")
ROI_SIG    = HERE.parent / "cell_typing" / "roi_signatures.csv"

SECTION_LABEL = {
    "a": "High-risk Tumor",
    "b": "Low-risk Tumor",
    "c": "High-risk T-cell",
    "d": "Low-risk T-cell",
    "t": "Middle-risk Tumor (ctrl)",
}
SECTION_COLOR = {"a": "#d62728", "b": "#1f77b4", "c": "#2ca02c",
                 "d": "#9467bd", "t": "#7f7f7f"}
SECTION_ORDER = ["a", "b", "c", "d", "t"]

# Pre-registered markers (from proteomics_분석.pdf + existing findings)
MARKER_HYPOTHESES = {
    "Tumor a vs b (high vs low)": [
        ("KIF20A",  "a>b"),
        ("KIF22",   "a>b"),
        ("INCENP",  "a>b"),
        ("MYH11",   "a>b"),
        ("TAGLN",   "a>b"),
        ("NCAM1",   "a>b"),
        ("APOBEC3C","a>b"),
    ],
}


# ---- I/O ----

def load_gg_matrix():
    m = pd.read_csv(GG_MATRIX, sep="\t")
    slide1_cols = [c for c in m.columns[3:] if c[0] in "abcdt"]
    # de-dupe stable order, but here all unique
    return m, slide1_cols


def sample_section(col):
    return col[0]


# ---- main pipeline ----

def log2_intensity(m, slide1_cols):
    # raw intensities, NaN preserved
    raw = m[slide1_cols].copy()
    # log2 with NaN preserved; treat 0 as NaN (zero detection ≠ low)
    raw[raw == 0] = np.nan
    return np.log2(raw)


def quality_filter(log2_df, slide1_cols, sec_a, sec_b, min_detect=0.30):
    """Keep genes with ≥ min_detect non-NaN fraction in EACH of sec_a and sec_b."""
    a_cols = [c for c in slide1_cols if c.startswith(sec_a)]
    b_cols = [c for c in slide1_cols if c.startswith(sec_b)]
    fa = log2_df[a_cols].notna().mean(axis=1)
    fb = log2_df[b_cols].notna().mean(axis=1)
    keep = (fa >= min_detect) & (fb >= min_detect)
    return keep, len(a_cols), len(b_cols)


def per_gene_mw(log2_df, genes, slide1_cols, sec_a, sec_b):
    a_cols = [c for c in slide1_cols if c.startswith(sec_a)]
    b_cols = [c for c in slide1_cols if c.startswith(sec_b)]
    rows = []
    for i, g in enumerate(genes):
        ra = log2_df.iloc[i][a_cols].dropna().values
        rb = log2_df.iloc[i][b_cols].dropna().values
        if len(ra) < 3 or len(rb) < 3:
            rows.append({"gene": g, "n_a": len(ra), "n_b": len(rb),
                         "log2_mean_a": float(ra.mean()) if len(ra) else np.nan,
                         "log2_mean_b": float(rb.mean()) if len(rb) else np.nan,
                         "log2_fc": np.nan,
                         "U": np.nan, "p": np.nan})
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
        observed = "a>b" if r["log2_fc"] > 0 else "a<b"
        rows.append({"gene": gene, "predicted_direction": predicted,
                     "observed_direction": observed,
                     "matches_hypothesis": predicted == observed,
                     "log2_fc": float(r["log2_fc"]),
                     "p": float(r["p"]),
                     "p_bh": float(r.get("p_bh", np.nan))})
    return pd.DataFrame(rows)


# ---- PCA / plots ----

def pca_samples(log2_df, slide1_cols, n_components=3):
    """PCA on samples — impute NaN with row mean (gene mean across samples)."""
    X = log2_df[slide1_cols].copy()
    row_means = X.mean(axis=1)
    X = X.apply(lambda col: col.fillna(row_means))     # NaN → gene mean
    Xt = X.T  # sample × gene
    Xt = Xt.fillna(Xt.mean(axis=0))                    # any residual NaN
    pca = PCA(n_components=n_components)
    PCs = pca.fit_transform(Xt.values)
    pcs_df = pd.DataFrame(PCs,
                          columns=[f"PC{i+1}" for i in range(n_components)],
                          index=Xt.index)
    pcs_df["section"]       = [c[0] for c in pcs_df.index]
    pcs_df["section_label"] = [SECTION_LABEL[s] for s in pcs_df["section"]]
    pcs_df["explained_var"] = ", ".join(f"{v*100:.1f}%" for v in pca.explained_variance_ratio_)
    return pcs_df, pca.explained_variance_ratio_


def plot_pca(pcs_df, evr, out_path):
    fig, ax = plt.subplots(figsize=(9, 7))
    for s in SECTION_ORDER:
        sub = pcs_df[pcs_df.section == s]
        if len(sub) == 0:
            continue
        ax.scatter(sub.PC1, sub.PC2, s=130, c=SECTION_COLOR[s],
                   edgecolor="black", linewidth=0.4,
                   label=f"{SECTION_LABEL[s]} (n={len(sub)})", alpha=0.85)
        for idx, row in sub.iterrows():
            ax.annotate(idx, (row.PC1, row.PC2), fontsize=7,
                        ha="center", va="center")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% var)")
    ax.set_title("PCA of slide1 proteomics samples — coloured by section",
                 fontsize=12)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(per_gene_df, title, out_path, label_top=12):
    df = per_gene_df.dropna(subset=["log2_fc", "p"]).copy()
    df["neglog10_p"] = -np.log10(df["p"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(10, 7))
    sig = df["p_bh"] < 0.05
    ax.scatter(df.loc[~sig, "log2_fc"], df.loc[~sig, "neglog10_p"],
               s=8, c="#bbbbbb", alpha=0.6)
    ax.scatter(df.loc[sig, "log2_fc"], df.loc[sig, "neglog10_p"],
               s=14, c="#d62728", alpha=0.85, label=f"BH-FDR < 0.05 (n={sig.sum()})")
    # annotate top
    top_up = df[sig & (df["log2_fc"] > 0)].nsmallest(label_top, "p")
    top_dn = df[sig & (df["log2_fc"] < 0)].nsmallest(label_top, "p")
    for _, r in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(r.gene, (r.log2_fc, r.neglog10_p), fontsize=7,
                    ha="center", va="center")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="grey", linestyle="-", linewidth=0.8)
    ax.set_xlabel("log2 fold change (a − b)")
    ax.set_ylabel("-log10(raw p)")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_top_markers_heatmap(per_gene_df, log2_df, slide1_cols, out_path,
                             top_n=20):
    top = per_gene_df.head(top_n)
    if len(top) == 0:
        return
    rows = top["gene"].tolist()
    matrix = log2_df.loc[log2_df.index[log2_df.index.isin(top.index)], slide1_cols]
    # gene-wise z-score for visual normalization
    m = matrix.values.copy()
    means = np.nanmean(m, axis=1, keepdims=True)
    stds  = np.nanstd(m, axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    Z = (m - means) / stds
    # sort columns by section then index
    col_order = sorted(slide1_cols, key=lambda c: (SECTION_ORDER.index(c[0]),
                                                    int(c[1:]) if c[1:].isdigit() else 0))
    col_idx = [slide1_cols.index(c) for c in col_order]
    Z = Z[:, col_idx]

    fig, ax = plt.subplots(figsize=(15, max(6, 0.35 * len(rows))))
    im = ax.imshow(Z, aspect="auto", cmap="vlag", vmin=-2.5, vmax=2.5)
    ax.set_xticks(np.arange(len(col_order)))
    ax.set_xticklabels(col_order, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows, fontsize=8)
    # section bar
    for i, c in enumerate(col_order):
        ax.add_patch(plt.Rectangle((i-0.5, -1.5), 1, 1,
                                    color=SECTION_COLOR[c[0]], clip_on=False))
    ax.set_title(f"Top-{len(rows)} markers (BH-sorted) — sample×gene z-score (a vs b)",
                 fontsize=11, pad=20)
    plt.colorbar(im, ax=ax, fraction=0.02, label="z-score")
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                          markerfacecolor=SECTION_COLOR[s], markersize=10,
                          label=SECTION_LABEL[s]) for s in SECTION_ORDER]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.18, 1.0),
              fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_section_protein_summary(log2_df, slide1_cols, out_path):
    n_detect = (log2_df[slide1_cols].notna()).sum(axis=0)
    med_int  = log2_df[slide1_cols].median(axis=0)
    sec = [c[0] for c in slide1_cols]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    for ax, vals, ylabel in zip(axes,
                                [n_detect.values, med_int.values],
                                ["# detected genes (non-NaN)",
                                 "median log2 intensity"]):
        data = [vals[[i for i, s in enumerate(sec) if s == ss]]
                for ss in SECTION_ORDER]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        labels=[SECTION_LABEL[s].replace(" ", "\n")
                                for s in SECTION_ORDER])
        for patch, s in zip(bp["boxes"], SECTION_ORDER):
            patch.set_facecolor(SECTION_COLOR[s]); patch.set_alpha(0.55)
        for i, s in enumerate(SECTION_ORDER):
            ys = [v for v, sec_ in zip(vals, sec) if sec_ == s]
            xs = np.random.normal(loc=i+1, scale=0.04, size=len(ys))
            ax.scatter(xs, ys, s=14, c="black", alpha=0.55, zorder=3)
        ax.set_ylabel(ylabel)
    fig.suptitle("Per-section sample quality — slide1 proteomics", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---- main ----

def main():
    print(f"[load] gg_matrix")
    m, slide1_cols = load_gg_matrix()
    print(f"       total genes: {len(m)}, slide1 samples: {len(slide1_cols)}")
    sec_counts = {s: sum(1 for c in slide1_cols if c[0] == s) for s in SECTION_ORDER}
    print(f"       per-section sample counts: {sec_counts}")

    # log2 transform
    log2_df = log2_intensity(m, slide1_cols)
    log2_df.index = m["Genes"].astype(str)        # gene as row label
    # also append the gene metadata
    print(f"       log2 matrix range: {np.nanmin(log2_df.values):.2f} .. "
          f"{np.nanmax(log2_df.values):.2f}")

    # per-sample summary
    summary = pd.DataFrame({
        "sample": slide1_cols,
        "section": [c[0] for c in slide1_cols],
        "section_label": [SECTION_LABEL[c[0]] for c in slide1_cols],
        "n_detected": [int((~log2_df[c].isna()).sum()) for c in slide1_cols],
        "median_log2": [float(log2_df[c].median(skipna=True))
                         for c in slide1_cols],
    })
    summary.to_csv(HERE / "slide1_columns_summary.csv", index=False)

    # ── Comparison 1: Tumor a vs b ──
    print(f"\n[stat] Tumor a vs b (High-risk vs Low-risk)")
    keep_ab, n_a, n_b = quality_filter(log2_df, slide1_cols, "a", "b")
    print(f"      genes passing detect-≥30% in both groups: {int(keep_ab.sum())}/{len(log2_df)}  "
          f"(n_a={n_a} samples, n_b={n_b} samples)")
    log2_ab = log2_df[keep_ab]
    per_g_ab = per_gene_mw(log2_ab, log2_ab.index.tolist(),
                           slide1_cols, "a", "b")
    per_g_ab.to_csv(HERE / "tumor_a_vs_b_genes.csv", index=False)
    sig_ab = int((per_g_ab["p_bh"] < 0.05).sum())
    print(f"      genes with BH-FDR<0.05: {sig_ab}")
    print(f"      top 10 by raw p:")
    print(per_g_ab.head(10)[["gene","n_a","n_b","log2_fc","p","p_bh"]]
          .to_string(index=False))

    # ── Comparison 2: T-cell c vs d ──
    print(f"\n[stat] T-cell c vs d (High-risk vs Low-risk)")
    keep_cd, n_c, n_d = quality_filter(log2_df, slide1_cols, "c", "d")
    print(f"      genes passing: {int(keep_cd.sum())}/{len(log2_df)}  "
          f"(n_c={n_c}, n_d={n_d})")
    log2_cd = log2_df[keep_cd]
    per_g_cd = per_gene_mw(log2_cd, log2_cd.index.tolist(),
                           slide1_cols, "c", "d")
    per_g_cd.to_csv(HERE / "tcell_c_vs_d_genes.csv", index=False)
    sig_cd = int((per_g_cd["p_bh"] < 0.05).sum())
    print(f"      genes with BH-FDR<0.05: {sig_cd}")
    print(f"      top 10:")
    print(per_g_cd.head(10)[["gene","n_c"if False else "n_a", "n_b", "log2_fc","p","p_bh"]]
          .to_string(index=False))

    # ── Marker hypothesis check ──
    print(f"\n[chek] pre-registered marker hypothesis check (Tumor a vs b)")
    mk = marker_check(per_g_ab, MARKER_HYPOTHESES["Tumor a vs b (high vs low)"])
    mk.to_csv(HERE / "marker_hypothesis_check.csv", index=False)
    print(mk.to_string(index=False))

    # ── PCA ──
    print(f"\n[pca]  sample PCA")
    pcs_df, evr = pca_samples(log2_df.loc[keep_ab], slide1_cols)
    pcs_df.to_csv(HERE / "pca_samples.csv")
    print(f"      explained variance (PC1-PC3): "
          f"{', '.join(f'{v*100:.1f}%' for v in evr)}")

    # ── Save the cleaned log2 matrix ──
    log2_df.loc[keep_ab | keep_cd].to_csv(HERE / "log2_intensity_matrix.csv")

    # ── Plots ──
    print(f"\n[plot] PCA / volcano / heatmap / sample-quality")
    plot_pca(pcs_df, evr, HERE / "pca_samples.png")
    plot_volcano(per_g_ab, "Tumor: High-risk (a) vs Low-risk (b)",
                 HERE / "volcano_tumor_a_vs_b.png")
    plot_volcano(per_g_cd, "T-cell: High-risk (c) vs Low-risk (d)",
                 HERE / "volcano_tcell_c_vs_d.png")
    plot_top_markers_heatmap(per_g_ab, log2_ab, slide1_cols,
                             HERE / "top_markers_heatmap.png")
    plot_section_protein_summary(log2_df, slide1_cols,
                                 HERE / "section_protein_summary.png")

    print(f"\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
