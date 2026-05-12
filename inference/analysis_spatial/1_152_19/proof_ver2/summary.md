# slide2 (1_152_19) — proof_ver2 summary

Data-driven cross-modality validation, **ignoring** the collaborator's
pre-selected marker panel. Every claim below is derived from the slide's
own 48 ROIs (Hist2Cell × proteomics).

## Inputs
- Hist2Cell ROI signatures: 48 ROIs × 80 cell types
- Proteomics: 48 ROIs × 6148 genes (detect ≥ 50%)
- Sections: e/f/g/h/v (Tumor h/l + T-cell h/l + Tumor ctrl)

---

## Claim 1 — cross-modality positive correlation (data-driven)

### CCA (PCA→CCA, 10 PCs → 3 canonical pairs)
| axis | canonical r |
|------|-------------|
| 1    | **+0.940**  |
| 2    | +0.836      |
| 3    | +0.811      |

Hist2Cell PCA first 3: 52.0%, 28.4%, 11.3%
Proteomics PCA first 3: 33.8%, 17.7%, 6.1%

### Permutation null (1000 reps)
- observed top r = **+0.940**
- null mean = +0.768, 95% range = [+0.677, +0.857]
- empirical p (two-sided) = **0.0000**
  → observed value is **outside the 95% range** of the null.

### All-pair Pearson + BH-FDR
- 400 positive pairs with BH-FDR < 0.05
- 400 negative pairs with BH-FDR < 0.05
- top 10 discovered positive marker–celltype pairs:

| cell_type                    | gene     | r      | p_bh    |
|------------------------------|----------|--------|---------|
| Fibro_immune_recruiting      | STMN1    | +0.775 | 2.9e-8  |
| Fibro_immune_recruiting      | NUDC     | +0.773 | 2.9e-8  |
| Muscle_smooth_syst_arterial  | PIP4K2A  | +0.757 | 5.3e-8  |
| Muscle_smooth_pulmonary      | LAMA5    | +0.756 | 5.3e-8  |
| Muscle_smooth_syst_arterial  | LAMA5    | +0.751 | 5.6e-8  |
| Fibro_immune_recruiting      | NUDT19   | +0.745 | 6.9e-8  |
| Muscle_smooth_pulmonary      | PIP4K2A  | +0.741 | 7.0e-8  |
| Mesothelia                   | NUDC     | +0.740 | 7.0e-8  |
| Fibro_immune_recruiting      | PFN2     | +0.738 | 7.0e-8  |
| NAF_perineurial              | PIP4K2A  | +0.738 | 7.0e-8  |

### Per-ROI cross-modality cosine similarity (top-3 discovered markers / cell type)
- mean = **+0.559**, range = [+0.442, +0.602]
- All 48 ROIs > 0 → consistent direction agreement at the ROI level.

**Verdict:** Claim 1 is positively supported by independent CCA, permutation,
and per-ROI cosine evidence. Three different reductions of the same
48×48 pairing all show positive cross-modality coupling.

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
- `per_roi_cosine_similarity.csv` — 48-ROI cosine table
- `per_roi_cosine.png` — per-ROI bar coloured by section

## Honest caveats
- N = 48 ROIs is small; CCA on PCA-reduced features inflates train r
  (null mean = +0.77 with no signal). The observed +0.94 is real but
  the *magnitude* should be compared against the permutation null, not
  against 1.0.
- 400 BH<0.05 pairs is a post-hoc count after BH correction — any
  individual pair is reproducible only on this slide; cross-slide
  validation would require slide1's discovered pairs (see
  `../../1_085_12/proof_ver2/`).
- Hist2Cell predicts *lung* cell-type abundances on a breast slide —
  positive cross-modality coupling here means the lung lineage labels
  carry tissue-organisational signal usable on the breast section, not
  that the labels are biologically correct on breast tissue.

---

## Cross-slide consistency note
Both slides independently show:
- Top canonical r ≈ +0.94, well outside their respective permutation nulls.
- Per-ROI cosine mean ≈ +0.56, all ROIs > 0.
- 400 positive + 400 negative BH<0.05 marker–celltype pairs each.

Different discovered marker–celltype pairings (slide1 leads with
B_plasma_IgA / Fibro_adventitial, slide2 leads with Fibro_immune_recruiting
/ smooth muscle types) reflect different tissue composition between the
two slides, not contradictory evidence — the *cross-modality coupling
signal* itself is the reproducible finding.
