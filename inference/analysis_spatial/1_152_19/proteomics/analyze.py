"""Proteomics differential analysis for slide2 (1_152_19) — full PNG set
matching slide1's heavy version (PCA / volcanos / top markers heatmap /
section protein summary).

Slide2 columns in gg_matrix: e1-e13 + f1-f15 + g1-g7 + h1-h8 + v1-v5 = 48.
(Duplicate 'e2' already removed at the source — commit a03c5d0.)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


HERE = Path(__file__).resolve().parent
GG   = Path("/home/sjhong/hist2cell/inference/analysis_spatial/report.gg_matrix (1).tsv")

SECTION_LABEL = {
    "e": "High-risk Tumor", "f": "Low-risk Tumor",
    "g": "High-risk T-cell", "h": "Low-risk T-cell",
    "v": "Middle-risk Tumor (ctrl)", "w": "Middle-risk T-cell (ctrl)",
}
SECTION_COLOR = {"e":"#d62728","f":"#1f77b4","g":"#2ca02c","h":"#9467bd","v":"#7f7f7f","w":"#bcbd22"}
SECTION_ORDER = ["e", "f", "g", "h", "v"]

MARKER_HYPOTHESES = [
    ("KIF20A",  "e>f"), ("KIF22", "e>f"), ("INCENP", "e>f"),
    ("MYH11",   "e>f"), ("TAGLN", "e>f"),
    ("NCAM1",   "e>f"), ("APOBEC3C", "e>f"),
]


def quality_filter(log2_df, cols, sa, sb, min_detect=0.30):
    a_cols = [c for c in cols if c[0] == sa]
    b_cols = [c for c in cols if c[0] == sb]
    fa = log2_df[a_cols].notna().mean(axis=1)
    fb = log2_df[b_cols].notna().mean(axis=1)
    return (fa >= min_detect) & (fb >= min_detect), len(a_cols), len(b_cols)


def per_gene_mw(log2_df, slide2_cols, sa, sb):
    a_cols = [c for c in slide2_cols if c[0] == sa]
    b_cols = [c for c in slide2_cols if c[0] == sb]
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
    df.loc[valid, "p_bh"] = multipletests(df.loc[valid, "p"], method="fdr_bh")[1]
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


def pca_samples(log2_df, slide2_cols, n_components=3):
    X = log2_df[slide2_cols].copy()
    row_means = X.mean(axis=1)
    X = X.apply(lambda col: col.fillna(row_means))
    Xt = X.T
    Xt = Xt.fillna(Xt.mean(axis=0))
    pca = PCA(n_components=n_components)
    PCs = pca.fit_transform(Xt.values)
    pcs = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(n_components)],
                       index=Xt.index)
    pcs["section"] = [c[0] for c in pcs.index]
    pcs["section_label"] = [SECTION_LABEL[s] for s in pcs["section"]]
    return pcs, pca.explained_variance_ratio_


def plot_pca(pcs, evr, out):
    fig, ax = plt.subplots(figsize=(9, 7))
    for s in SECTION_ORDER:
        sub = pcs[pcs.section == s]
        if len(sub) == 0: continue
        ax.scatter(sub.PC1, sub.PC2, s=130, c=SECTION_COLOR[s],
                   edgecolor="black", linewidth=0.4,
                   label=f"{SECTION_LABEL[s]} (n={len(sub)})", alpha=0.85)
        for idx, row in sub.iterrows():
            ax.annotate(idx, (row.PC1, row.PC2), fontsize=7,
                        ha="center", va="center")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% var)")
    ax.set_title("PCA of slide2 proteomics samples — coloured by section", fontsize=12)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_volcano(per_gene_df, title, out, label_top=12):
    df = per_gene_df.dropna(subset=["log2_fc", "p"]).copy()
    df["neglog10_p"] = -np.log10(df["p"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(10, 7))
    sig = df["p_bh"] < 0.05
    ax.scatter(df.loc[~sig, "log2_fc"], df.loc[~sig, "neglog10_p"],
               s=8, c="#bbbbbb", alpha=0.6)
    ax.scatter(df.loc[sig, "log2_fc"], df.loc[sig, "neglog10_p"],
               s=14, c="#d62728", alpha=0.85,
               label=f"BH-FDR < 0.05 (n={sig.sum()})")
    top_up = df[(df["log2_fc"] > 0)].nsmallest(label_top, "p")
    top_dn = df[(df["log2_fc"] < 0)].nsmallest(label_top, "p")
    for _, r in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(r.gene, (r.log2_fc, r.neglog10_p), fontsize=7,
                    ha="center", va="center")
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="grey", linestyle="-", linewidth=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(raw p)")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="best")
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_top_markers_heatmap(per_gene_df, log2_df, slide2_cols, out, top_n=20):
    top = per_gene_df.head(top_n)
    rows = top["gene"].tolist()
    log2_unique = log2_df[~log2_df.index.duplicated(keep="first")]
    present = [g for g in rows if g in log2_unique.index]
    matrix = log2_unique.loc[present, slide2_cols]
    rows = present
    m = matrix.values.astype(float).copy()
    means = np.nanmean(m, axis=1, keepdims=True)
    stds = np.nanstd(m, axis=1, keepdims=True); stds[stds == 0] = 1.0
    Z = (m - means) / stds
    Z = np.ma.masked_invalid(Z)
    col_order = sorted(slide2_cols,
                       key=lambda c: (SECTION_ORDER.index(c[0]) if c[0] in SECTION_ORDER else 99,
                                      int(c[1:]) if c[1:].isdigit() else 0))
    col_idx = [slide2_cols.index(c) for c in col_order]
    Z = Z[:, col_idx]

    fig, ax = plt.subplots(figsize=(15, max(6, 0.35 * len(rows))))
    cmap = plt.get_cmap("vlag").copy()
    cmap.set_bad(color="#dddddd")
    im = ax.imshow(Z, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5)
    ax.set_xticks(np.arange(len(col_order)))
    ax.set_xticklabels(col_order, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows, fontsize=8)
    for i, c in enumerate(col_order):
        ax.add_patch(plt.Rectangle((i-0.5, -1.5), 1, 1,
                                    color=SECTION_COLOR[c[0]], clip_on=False))
    ax.set_title(f"Top-{len(rows)} markers (BH-sorted) — slide2 sample×gene z-score (e vs f)",
                 fontsize=11, pad=20)
    plt.colorbar(im, ax=ax, fraction=0.02, label="z-score")
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                          markerfacecolor=SECTION_COLOR[s], markersize=10,
                          label=SECTION_LABEL[s]) for s in SECTION_ORDER]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.18, 1.0),
              fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_section_protein_summary(log2_df, slide2_cols, out):
    n_detect = (log2_df[slide2_cols].notna()).sum(axis=0)
    med_int = log2_df[slide2_cols].median(axis=0)
    sec = [c[0] for c in slide2_cols]
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
    fig.suptitle("Per-section sample quality — slide2 proteomics", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    print("[load] gg_matrix slide2 columns")
    m = pd.read_csv(GG, sep="\t")
    slide2_cols = [c for c in m.columns[3:] if c[0] in "efghv"]
    print(f"  total genes: {len(m)}, slide2 samples: {len(slide2_cols)}")
    raw = m[slide2_cols].replace(0, np.nan)
    log2 = np.log2(raw); log2.index = m["Genes"].astype(str)

    summary = pd.DataFrame({
        "sample": slide2_cols,
        "section": [c[0] for c in slide2_cols],
        "section_label": [SECTION_LABEL[c[0]] for c in slide2_cols],
        "n_detected": [int((~log2[c].isna()).sum()) for c in slide2_cols],
        "median_log2": [float(log2[c].median(skipna=True)) for c in slide2_cols],
    })
    summary.to_csv(HERE/"slide2_columns_summary.csv", index=False)

    print("\n[1] Tumor e vs f")
    keep_ef, n_e, n_f = quality_filter(log2, slide2_cols, "e", "f")
    print(f"  pass detect-≥30%: {int(keep_ef.sum())}/{len(log2)} (n_e={n_e}, n_f={n_f})")
    per_ef = per_gene_mw(log2[keep_ef], slide2_cols, "e", "f")
    per_ef.to_csv(HERE/"tumor_e_vs_f_genes.csv", index=False)
    print(f"  BH<0.05: {int((per_ef.p_bh<.05).sum())}")

    print("\n[2] T-cell g vs h")
    keep_gh, n_g, n_h = quality_filter(log2, slide2_cols, "g", "h")
    per_gh = per_gene_mw(log2[keep_gh], slide2_cols, "g", "h")
    per_gh.to_csv(HERE/"tcell_g_vs_h_genes.csv", index=False)
    print(f"  pass: {int(keep_gh.sum())} / BH<0.05: {int((per_gh.p_bh<.05).sum())}")

    print("\n[3] marker hypothesis check")
    mk = marker_check(per_ef, MARKER_HYPOTHESES)
    mk.to_csv(HERE/"marker_hypothesis_check.csv", index=False)
    print(mk.to_string(index=False))

    print("\n[4] PCA")
    pcs, evr = pca_samples(log2.loc[keep_ef], slide2_cols)
    pcs.to_csv(HERE/"pca_samples.csv")
    print(f"  evr PC1-3: {', '.join(f'{v*100:.1f}%' for v in evr)}")

    # Save log2 matrix for downstream use
    log2.loc[keep_ef | keep_gh].to_csv(HERE/"log2_intensity_matrix.csv")

    print("\n[5] plots")
    plot_pca(pcs, evr, HERE/"pca_samples.png")
    plot_volcano(per_ef, "slide2 Tumor: High-risk (e) vs Low-risk (f)",
                 HERE/"volcano_tumor_e_vs_f.png")
    plot_volcano(per_gh, "slide2 T-cell: High-risk (g) vs Low-risk (h)",
                 HERE/"volcano_tcell_g_vs_h.png")
    plot_top_markers_heatmap(per_ef, log2[keep_ef], slide2_cols,
                             HERE/"top_markers_heatmap.png")
    plot_section_protein_summary(log2, slide2_cols,
                                 HERE/"section_protein_summary.png")

    print("\nDone. Outputs:")
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.suffix in {".csv", ".png", ".py"}:
            print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
