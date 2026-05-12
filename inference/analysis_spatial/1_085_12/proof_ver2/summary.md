# slide1 (1_085_12) — proof_ver2 summary

Data-driven cross-modality validation, **ignoring** the collaborator's
pre-selected marker panel. Every claim below is derived from the slide's
own 46 ROIs (Hist2Cell × proteomics).

## Inputs
- Hist2Cell ROI signatures: 46 ROIs × 80 cell types
- Proteomics: 46 ROIs × 4216 genes (detect ≥ 50%)
- Sections: a/b/c/d/t (Tumor h/l + T-cell h/l + Tumor ctrl)

---

## Claim 1 — cross-modality positive correlation (data-driven)

### CCA (PCA→CCA, 10 PCs → 3 canonical pairs)
| axis | canonical r |
|------|-------------|
| 1    | **+0.936**  |
| 2    | +0.875      |
| 3    | +0.710      |

Hist2Cell PCA first 3: 55.1%, 20.5%, 15.3%
Proteomics PCA first 3: 22.6%, 19.9%, 7.3%

### Permutation null (1000 reps)
- observed top r = **+0.936**
- null mean = +0.778, 95% range = [+0.683, +0.863]
- empirical p (two-sided) = **0.0000**
  → observed value is **outside the 95% range** of the null.

### All-pair Pearson + BH-FDR
- 400 positive pairs with BH-FDR < 0.05
- 400 negative pairs with BH-FDR < 0.05
- top 10 discovered positive marker–celltype pairs:

| cell_type            | gene     | r      | p_bh    |
|----------------------|----------|--------|---------|
| B_plasma_IgA         | HSPA1L   | +0.730 | 6.2e-7  |
| Fibro_adventitial    | PRDX6    | +0.730 | 6.2e-7  |
| B_plasma_IgA         | SLC25A13 | +0.719 | 1.1e-6  |
| Fibro_adventitial    | NME2     | +0.714 | 1.2e-6  |
| B_plasma_IgA         | DDX3X    | +0.713 | 1.2e-6  |
| Chondrocyte          | DBN1     | +0.712 | 1.2e-6  |
| AT2                  | COLGALT1 | +0.710 | 1.2e-6  |
| Chondrocyte          | CDH1     | +0.710 | 1.2e-6  |
| Macro_AW_CX3CR1      | DDX3X    | +0.705 | 1.5e-6  |
| Macro_AW_CX3CR1      | CC2D1A   | +0.702 | 1.6e-6  |

### Per-ROI cross-modality cosine similarity (top-3 discovered markers / cell type)
- mean = **+0.555**, range = [+0.491, +0.600]
- All 46 ROIs > 0 → consistent direction agreement at the ROI level.

**Verdict:** Claim 1 is positively supported by independent CCA, permutation,
and per-ROI cosine evidence. Three different reductions of the same
46×46 pairing all show positive cross-modality coupling.

---

## Claim 2 — per-ROI top cell types
Already produced under the prior pipeline at
[`../proofs/roi_top_celltypes.csv`](../proofs/roi_top_celltypes.csv)
and [`../proofs/roi_top_celltypes_heatmap.png`](../proofs/roi_top_celltypes_heatmap.png).
No re-computation needed for Claim 2; the data-driven proof_ver2 work
is added evidence for Claim 1.

---

## Outputs in this folder
- `cca_summary.csv` — train r's + permutation 95% range + p
- `cca_scatter.png` — 3-axis canonical scatter, points coloured by section
- `permutation_null.png` — null histogram + observed line
- `cca_loadings_axis1.png` — top ± Hist2Cell celltype / proteomics gene loaders
- `discovered_marker_pairs.csv` — 800 BH<0.05 marker–celltype pairs
- `top_discovered_pairs.png` — top-20 positive bar
- `per_roi_cosine_similarity.csv` — 46-ROI cosine table
- `per_roi_cosine.png` — per-ROI bar coloured by section

## Honest caveats
- N = 46 ROIs is small; CCA on PCA-reduced features inflates train r
  (null mean = +0.78 with no signal). The observed +0.94 is real but
  the *magnitude* should be compared against the permutation null, not
  against 1.0.
- 400 BH<0.05 pairs is a post-hoc count after BH correction — any
  individual pair is reproducible only on this slide; cross-slide
  validation would require slide2's discovered pairs (see
  `../../1_152_19/proof_ver2/`).
- Hist2Cell predicts *lung* cell-type abundances on a breast slide —
  positive cross-modality coupling here means the lung lineage labels
  carry tissue-organisational signal usable on the breast section, not
  that the labels are biologically correct on breast tissue.
