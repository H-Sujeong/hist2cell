# slide2 (1_152_19) — proof_ver2 요약

## 이 분석이 무엇인지

협업자가 사전 선정한 마커 패널(MYH11, KRT8, COL1A1 …)을 **무시하고**,
이 슬라이드 자체의 Hist2Cell × proteomics 행렬에서 두 modality 사이의
공통 신호를 처음부터 다시 끌어내본 분석. 이전(`../proofs/`) 결과는
사전 마커에 묶여 약한 상관계수만 보여줬는데, 데이터 자체에 더 강한
교차-modality 결합이 있는지 검증하는 것이 목적이다.

분석 방식·해석 틀은 slide1 (`../../1_085_12/proof_ver2/summary.md`) 와
동일하다. 본 문서는 slide2 의 수치와 슬라이드별 차이만 기술한다.

## 분석 대상

| 항목 | 수치 |
|---|---|
| ROI 개수 (양쪽 modality 모두 있는 것) | 48 |
| Hist2Cell cell type 개수 | 80 |
| Proteomics gene 개수 (detect ≥ 50% 필터) | 6148 |
| Section 라벨 | e/f/g/h/v (Tumor h, Tumor l, T-cell h, T-cell l, Tumor ctrl) |

slide1 보다 proteomics 가 1.5배 가량 더 많은 유전자를 검출했고 (4216 →
6148), ROI 도 2개 더 많다 (46 → 48). 즉 slide2 는 신호 대비 자유도가
조금 더 여유로운 케이스.

---

## Claim 1 — "교차-modality 양의 상관관계가 존재한다" 의 데이터-기반 검증

### 1. CCA — 두 modality 의 공통 축

PCA 로 각 modality 를 10차원까지 축약 후 CCA. 3개 canonical pair:

| canonical axis | r (Pearson) |
|---|---|
| 1 | **+0.940** |
| 2 | +0.836 |
| 3 | +0.811 |

PCA 분산 (참고):
- Hist2Cell PC1‒3: 52.0% / 28.4% / 11.3% (= 91.7% 까지 설명)
- Proteomics PC1‒3: 33.8% / 17.7% / 6.1% — slide1 (22.6/19.9/7.3) 대비
  PC1 비중이 훨씬 크다. **slide2 proteomics 가 단일 방향으로 더 강하게
  변동한다**는 의미.

→ 두 modality 모두 PC1 이 분산을 압도하는 구조이고, 그 PC1 두 개가
거의 같이 움직인다 (canonical 1 r = +0.94). 그래프는 `cca_scatter.png`.

### 2. Permutation null — 우연 수준 비교

proteomics ROI 축 1000회 shuffle 후 CCA 재실행.

| 지표 | 값 |
|---|---|
| 관측 top canonical r | +0.940 |
| 영가설 평균 | +0.768 |
| 영가설 95% 구간 | [+0.677, +0.857] |
| empirical p-value (양측) | 0.0000 |

→ slide1 (관측 +0.936, null 평균 +0.778) 과 거의 동일한 패턴. **두
슬라이드가 서로 독립적으로 같은 결론**(관측치가 null 95% 구간 바깥) 을
주고 있다는 점이 이번 분석의 가장 강한 메시지 중 하나.

### 3. All-pair Pearson + BH-FDR — 어떤 cell type / gene 이 결합을 만드는가

80 × 6148 모든 (cell type × gene) 조합에서 ROI 48개의 Pearson r 계산
→ cell type 별 top-5 양/음 후보 800개 → BH-FDR 보정.

| 항목 | 값 |
|---|---|
| 후보 pair | 800 |
| BH-FDR < 0.05 통과 양의 pair | 400 / 400 |
| BH-FDR < 0.05 통과 음의 pair | 400 / 400 |

상위 10개 양의 pair:

| cell type | gene | r | p_bh |
|---|---|---|---|
| Fibro_immune_recruiting | STMN1 | +0.775 | 2.9e-8 |
| Fibro_immune_recruiting | NUDC | +0.773 | 2.9e-8 |
| Muscle_smooth_syst_arterial | PIP4K2A | +0.757 | 5.3e-8 |
| Muscle_smooth_pulmonary | LAMA5 | +0.756 | 5.3e-8 |
| Muscle_smooth_syst_arterial | LAMA5 | +0.751 | 5.6e-8 |
| Fibro_immune_recruiting | NUDT19 | +0.745 | 6.9e-8 |
| Muscle_smooth_pulmonary | PIP4K2A | +0.741 | 7.0e-8 |
| Mesothelia | NUDC | +0.740 | 7.0e-8 |
| Fibro_immune_recruiting | PFN2 | +0.738 | 7.0e-8 |
| NAF_perineurial | PIP4K2A | +0.738 | 7.0e-8 |

**해석.**
- 상위권이 **fibroblast / smooth muscle / mesothelia** 의 stromal 축으로
  쏠려 있다. slide1 이 B_plasma 같은 면역 형질세포 축이 가장 강했던
  것과 대비된다. 두 슬라이드의 *드라이브하는 조직 구성* 자체가 다르다는
  것을 보여줌 (slide2 의 tumor microenvironment 가 stromal-heavy).
- **Fibro_immune_recruiting** head 가 STMN1·NUDC·NUDT19·PFN2 와 강하게
  결합 — STMN1 / PFN2 모두 cytoskeleton-remodeling 단백질로, 활성
  fibroblast 의 *이동성·증식* 단백 시그니처와 일관됨.
- **Muscle_smooth_***(pulmonary, syst_arterial) head 가 LAMA5 (laminin
  α5, basement membrane) 와 PIP4K2A (phosphatidylinositol kinase) 에
  강하게 결합 — 평활근 - 기저막 단백 페어링은 생물학적으로 깨끗한 신호.
  lung 모델의 smooth muscle head 들이 breast 조직의 평활근/혈관주위
  영역을 잘 잡고 있다고 볼 수 있다.

### 4. Per-ROI cosine similarity

cell type 별 top-3 양의 마커로 합성한 proteomics 점수 vs Hist2Cell 벡터의
ROI 별 cosine.

| 지표 | 값 |
|---|---|
| 평균 cosine | **+0.559** |
| 범위 | [+0.442, +0.602] |
| 음수 ROI 개수 | **0 / 48** (전부 양수) |

Section 별 평균:

| section | label | n | 평균 cosine |
|---|---|---|---|
| e | High-risk Tumor | 13 | +0.563 |
| f | Low-risk Tumor | 15 | +0.575 |
| g | High-risk T-cell | 7 | +0.504 |
| h | Low-risk T-cell | 8 | +0.560 |
| v | Middle-risk Tumor (ctrl) | 5 | +0.577 |

→ 모든 ROI 가 양의 일치. section g (High-risk T-cell) 의 cosine 이 약간
낮은 것은 stromal 축에서 골라낸 마커가 T-cell-rich 영역에서는 자연히
신호가 약해지기 때문 — 그래도 여전히 모두 양수임.

### 종합 — Claim 1 결론

slide1 과 동일한 결론. 세 가지 reduction (CCA, all-pair Pearson, ROI
cosine) 이 모두 같은 방향으로 양의 cross-modality 결합을 가리킨다.
slide1·slide2 가 **각자 다른 ROI 들로 독립적으로** 같은 결론에 도달했
다는 사실이 단일 슬라이드의 우연으로 결론을 설명할 수 없게 만든다.

---

## Claim 2 — ROI 별 top cell type 목록

이전 파이프라인 산출물 그대로 사용:

- [`../proofs/roi_top_celltypes.csv`](../proofs/roi_top_celltypes.csv) —
  48 ROI × top1‒top5 cell type + lineage group
- [`../proofs/roi_top_celltypes_heatmap.png`](../proofs/roi_top_celltypes_heatmap.png) —
  ROI × union-of-top z-score 히트맵

---

## 산출 파일

| 파일 | 내용 |
|---|---|
| `cca_summary.csv` | 3개 canonical axis 의 train r + null 95% 구간 + p |
| `cca_scatter.png` | canonical pair 1/2/3 산점도 (section 색) |
| `permutation_null.png` | null 히스토그램 + 관측 r 수직선 |
| `cca_loadings_axis1.png` | canonical axis 1 의 ± top loader |
| `discovered_marker_pairs.csv` | 800 후보 pair (top-5 양/음 × 80 cell type, BH-FDR 포함) |
| `top_discovered_pairs.png` | 상위 20개 양의 pair 막대 |
| `per_roi_cosine_similarity.csv` | 48 ROI × cosine 점수 |
| `per_roi_cosine.png` | per-ROI cosine 막대 (section 색) |

---

## 두 슬라이드 비교 (slide1 vs slide2)

| 지표 | slide1 (1_085_12) | slide2 (1_152_19) |
|---|---|---|
| ROI 개수 | 46 | 48 |
| Proteomics gene 수 | 4216 | 6148 |
| CCA top r | +0.936 | +0.940 |
| Null 평균 | +0.778 | +0.768 |
| Null 95% 상한 | +0.863 | +0.857 |
| Per-ROI cosine 평균 | +0.555 | +0.559 |
| 음수 cosine ROI | 0 / 46 | 0 / 48 |
| 발견 1순위 cell type | B_plasma_IgA / Fibro_adventitial | Fibro_immune_recruiting / Muscle_smooth_* |

**핵심 관찰.**
- 전체 강도(CCA r, per-ROI cosine)는 슬라이드 사이에 매우 유사 —
  두 슬라이드가 같은 *cross-modality coupling 강도* 를 보여준다.
- 어떤 cell type / gene 이 결합을 주도하는지는 **다르다** — slide1 은
  면역 형질세포가, slide2 는 stromal (섬유아세포·평활근) 축이 가장 강하다.
  이는 두 종양의 미세환경 자체가 다르기 때문이지 모순되는 결과가 아니다.

---

## 솔직한 caveat

1. **N=48 의 한계.** CCA 자체가 작은 N 에서 train r 을 부풀린다.
   관측치는 *절대 크기* 가 아니라 **null 분포 대비 위치** 로 해석해야
   한다 (95% 구간 바깥, p < 1/1000).

2. **BH-FDR 범위.** `discovered_marker_pairs.csv` 의 BH 는 *cell type 별
   top-5* 로 추려진 800개 후보 위의 보정이다. 80×6148 ≈ 491k 전체 pair
   에 대한 전역 보정은 아니며, 그 경우 생존 pair 수는 적어질 수 있다.

3. **post-hoc 선택 편향.** 마커가 같은 행렬에서 학습·검증된다. 진정한
   외부 검증은 slide1 (`../../1_085_12/proof_ver2/`) 의 발견 pair 와
   비교하는 것 — 두 슬라이드 공통 pair 가 가장 신뢰할 만한 신호다.
   현 결과에서 발견된 상위 pair 는 슬라이드별로 갈리지만, *coupling 자체
   가 있다는 사실* 은 양쪽 슬라이드에서 독립적으로 재현되었다.

4. **Hist2Cell 모델의 도메인 갭.** lung 단일세포 라벨로 학습된 모델을
   breast 슬라이드에 적용한 상황. 양의 cross-modality 결합은 "lung
   라벨의 분류기 출력이 breast 조직 구성 차이를 의미있게 따라간다"는
   뜻이지, lung cell type 이 breast 에서 같은 정체성을 갖는다는 뜻이
   아니다. 라벨 해석은 "기능 모듈" 수준으로 받아들이는 게 안전.
