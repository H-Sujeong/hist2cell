# slide1 (1_085_12) — proof_ver2 요약

## 이 분석이 무엇인지

협업자가 사전 선정한 마커 패널(MYH11, KRT8, COL1A1 …)을 **무시하고**,
이 슬라이드 자체의 Hist2Cell × proteomics 행렬에서 두 modality 사이의
공통 신호를 처음부터 다시 끌어내본 분석. 이전(`../proofs/`) 결과는
사전 마커에 묶여 Pearson r 이 최대 +0.38 수준에서 멈췄던 것이
"correlation 이 약하다"는 인상을 줬는데, 데이터 자체에는 더 강한
교차-modality 결합이 있는지 검증하는 것이 목적이다.

## 분석 대상

| 항목 | 수치 |
|---|---|
| ROI 개수 (양쪽 modality 모두 있는 것) | 46 |
| Hist2Cell cell type 개수 | 80 |
| Proteomics gene 개수 (detect ≥ 50% 필터) | 4216 |
| Section 라벨 | a/b/c/d/t (Tumor h, Tumor l, T-cell h, T-cell l, Tumor ctrl) |

---

## Claim 1 — "교차-modality 양의 상관관계가 존재한다" 의 데이터-기반 검증

### 1. CCA (Canonical Correlation Analysis) — 두 modality 의 *공통 축* 찾기

**무엇을 했나.**
46×80 Hist2Cell 행렬, 46×4216 proteomics 행렬을 각각 PCA 로 10차원까지
축약한 뒤(고차원 두 modality 를 직접 CCA 에 넣으면 자유도가 모자라
trivial 한 r=1 이 나오기 때문에 PCA 가 필수), 두 modality 사이에서
**상관관계가 가장 큰 선형결합 쌍**을 찾는다. 정확히 이런 쌍 3개를
canonical pair 1/2/3 으로 보고한다.

**결과.**

| canonical axis | r (Pearson) |
|---|---|
| 1 | **+0.936** |
| 2 | +0.875 |
| 3 | +0.710 |

PCA 가 잡은 분산 (참고용):
- Hist2Cell PC1‒3: 55.1% / 20.5% / 15.3% (= 90.9% 까지 설명)
- Proteomics PC1‒3: 22.6% / 19.9% / 7.3% (proteomics 가 더 분산되어 있음)

→ **두 modality 가 단순히 같은 변수에 평행하게 의존하는 게 아니라,
공통의 잠재 축이 실제로 존재한다는 첫 번째 증거**.
그래프는 `cca_scatter.png`. 점을 section 별 색으로 찍었기 때문에,
canonical axis 1 이 Tumor vs T-cell 같은 큰 조직 구성 차이를 따라가는지도
시각적으로 확인 가능.

### 2. Permutation null — "이 r 이 정말 우연 이상인가"

**무엇을 했나.**
proteomics 의 ROI 축을 1000번 무작위로 섞은 뒤 매번 PCA→CCA 를 다시
돌려서, 신호가 없는 상태에서 자연히 나오는 top r 의 분포를 구함.
관측치와 비교하면 "관측된 +0.94 가 우연 분포의 어디에 있는가" 가 나옴.

**결과.**

| 지표 | 값 |
|---|---|
| 관측 top canonical r | +0.936 |
| 영가설 평균 | +0.778 |
| 영가설 95% 구간 | [+0.683, +0.863] |
| empirical p-value (양측) | 0.0000 (1000회 중 0회) |

→ CCA 자체가 작은 N (=46) 에서 *신호 없이도* 평균 +0.78 까지 부풀리는
경향이 있다는 점은 솔직히 인정해야 한다. 핵심은 +0.94 가 **null 의 95%
범위 바깥**에 있다는 것 — null 의 0.94 분위는 약 0.86 이므로 관측치는
그보다 멀리 떨어져 있다. 그래프는 `permutation_null.png` (회색 히스토
그램 + 빨간 수직선).

### 3. All-pair Pearson + BH-FDR — *어떤 유전자/세포타입* 이 결합을 만드나

**무엇을 했나.**
80 × 4216 의 모든 cell type × gene 조합에 대해 ROI 46개 위에서 Pearson r
을 계산. 그 중 각 cell type 별로 상관이 가장 큰 양/음 5개씩만 골라낸 뒤
이 800개 후보에 BH-FDR 다중검정 보정을 적용.

**결과.**

| 항목 | 값 |
|---|---|
| 후보 pair (top-5/cell type × 양음) | 800 |
| BH-FDR < 0.05 통과 양의 pair | 400 / 400 |
| BH-FDR < 0.05 통과 음의 pair | 400 / 400 |

상위 10개 양의 pair:

| cell type | gene | r | p_bh |
|---|---|---|---|
| B_plasma_IgA | HSPA1L | +0.730 | 6.2e-7 |
| Fibro_adventitial | PRDX6 | +0.730 | 6.2e-7 |
| B_plasma_IgA | SLC25A13 | +0.719 | 1.1e-6 |
| Fibro_adventitial | NME2 | +0.714 | 1.2e-6 |
| B_plasma_IgA | DDX3X | +0.713 | 1.2e-6 |
| Chondrocyte | DBN1 | +0.712 | 1.2e-6 |
| AT2 | COLGALT1 | +0.710 | 1.2e-6 |
| Chondrocyte | CDH1 | +0.710 | 1.2e-6 |
| Macro_AW_CX3CR1 | DDX3X | +0.705 | 1.5e-6 |
| Macro_AW_CX3CR1 | CC2D1A | +0.702 | 1.6e-6 |

**해석.**
- **B_plasma_IgA** ↔ HSPA1L / SLC25A13 / DDX3X — 점막성 IgA 형질세포가
  풍부한 ROI 가 일관되게 미토콘드리아·RNA 스트레스 단백질이 풍부한 ROI 와
  맞아떨어진다. lung 모델의 B-cell head 가 breast 조직에서도 plasma-cell
  rich 영역을 잡고 있다는 신호.
- **Fibro_adventitial** ↔ PRDX6 / NME2 — adventitial fibroblast head 가
  활성 fibroblast 영역(에너지대사·산화환원 단백 풍부)과 함께 움직임.
- **AT2** ↔ COLGALT1 — collagen 후번역 변형 효소가 AT2 head 와 같이
  움직이는 것은 breast 조직 맥락에서는 *상피 분비 활성* 으로 해석하는 게
  자연스럽다 (lung AT2 라벨이 breast 의 분비상피를 대신 잡는 케이스).
- 80 cell type 모두 BH<0.05 후보가 5개씩 살아남았다는 것은 **cell type
  마다 적어도 다수의 유전자와 강하게 상관된다**는 의미. 단, 이 BH 보정은
  *이미 top-5 로 추려진* 800개 위에서 돌린 것이지 80×4216 ≈ 337k 전체
  pair 에 대한 보정이 아니다 (caveat 참조).

### 4. Per-ROI cosine similarity — ROI 단위로도 일관되는가

**무엇을 했나.**
3번에서 cell type 별로 가장 강하게 양의 상관을 보인 마커 3개씩을 골라
proteomics 점수 벡터를 합성 (= "각 cell type 의 데이터-기반 단백질
점수"). 그 다음 ROI 별로 Hist2Cell 의 80-cell-type 벡터와 이 proteomics
점수 벡터의 cosine similarity 를 잰다. ROI 마다 두 modality 가 같은
방향을 가리키는지 보는 ROI-local check.

**결과.**

| 지표 | 값 |
|---|---|
| 평균 cosine | **+0.555** |
| 범위 | [+0.491, +0.600] |
| 음수인 ROI 개수 | **0 / 46** (전부 양수) |

Section 별 평균:

| section | label | n | 평균 cosine |
|---|---|---|---|
| a | High-risk Tumor | 8 | +0.524 |
| b | Low-risk Tumor | 21 | +0.563 |
| c | High-risk T-cell | 5 | +0.558 |
| d | Low-risk T-cell | 9 | +0.563 |
| t | Middle-risk Tumor (ctrl) | 3 | +0.549 |

→ Section 라벨과 무관하게 모든 ROI 가 양의 일치를 보임. CCA 의 큰 r
이 전역 통계의 artefact 가 아니라 **ROI 하나하나 단위에서도 두 modality
가 같은 방향을 가리킨다**는 것을 확인하는 sanity check 역할.

그래프는 `per_roi_cosine.png` (ROI bar, section 색).

### 종합 — Claim 1 결론

세 가지 서로 다른 reduction —
**(1)** 전체 행렬을 잠재공간으로 펼친 CCA,
**(2)** cell-type × gene 단위의 미세 상관관계 검정,
**(3)** ROI 단위의 cosine —
모두 양의 cross-modality 결합을 가리킨다.

이전 `../proofs/cross_modality_correlations.csv` 에서 본 r=+0.38 은
*사전에 사람이 고른 6개 마커 그룹 + 사전에 고른 cell type 그룹* 으로
사이즈를 줄여놓은 상태의 상관계수였다. 데이터 기반으로 다시 보면
**더 큰 결합 신호가 분명히 있다**는 게 이번 분석의 핵심 메시지.

---

## Claim 2 — ROI 별 top cell type 목록

이 부분은 이미 이전 파이프라인 `../proofs/` 에 산출되어 있어서
proof_ver2 에서 재계산하지 않았다. 참조 파일:

- [`../proofs/roi_top_celltypes.csv`](../proofs/roi_top_celltypes.csv) —
  46 ROI × top1‒top5 cell type + lineage group
- [`../proofs/roi_top_celltypes_heatmap.png`](../proofs/roi_top_celltypes_heatmap.png) —
  ROI × union-of-top z-score 히트맵

proof_ver2 는 Claim 1 (교차-modality 양의 상관) 에 대해 더 강한 증거를
얹는 역할이다.

---

## 산출 파일

| 파일 | 내용 |
|---|---|
| `cca_summary.csv` | 3개 canonical axis 의 train r + null 95% 구간 + p |
| `cca_scatter.png` | canonical pair 1/2/3 산점도 (section 색) |
| `permutation_null.png` | null 히스토그램 + 관측 r 수직선 |
| `cca_loadings_axis1.png` | canonical axis 1 의 ± top loader (Hist2Cell 셀타입 / proteomics 유전자) |
| `discovered_marker_pairs.csv` | 800 후보 pair (각 cell type 별 top-5 양/음, BH-FDR 포함) |
| `top_discovered_pairs.png` | 상위 20개 양의 pair 막대 |
| `per_roi_cosine_similarity.csv` | 46 ROI × cosine 점수 |
| `per_roi_cosine.png` | per-ROI cosine 막대 (section 색) |

---

## 솔직한 caveat

1. **N=46 의 한계.** CCA 는 작은 N 에서 신호가 없어도 자연히 train r 이
   부풀려진다. null 평균 +0.78 이 그것을 보여준다. 따라서 관측 +0.94 의
   *절대치* 자체에 의미부여하지 말고, **null 분포 대비 얼마나 멀리 있느냐**
   로 읽어야 한다 (95% 구간 바깥, p < 1/1000).

2. **BH-FDR 범위.** `discovered_marker_pairs.csv` 의 BH 는 *각 cell type
   별 top-5* 로 추려진 800개 후보 위에 적용된 것이다. 즉 "이미 가장 강한
   상관관계만 들어간 pool 안에서의 보정". 더 엄밀하게는 80×4216 ≈ 337k
   전체 pair 에 BH 를 적용해야 하며, 그 경우 살아남는 pair 수는
   현저히 줄어들 수 있다. 본 분석은 "각 cell type 마다 어느 정도 결합되는
   유전자가 적어도 일부 있는지 빠르게 훑는" 목적의 검정이다.

3. **post-hoc 선택 편향.** 발견된 마커 pair 는 같은 ROI 행렬에서 학습되고
   같은 행렬에서 검증된다. 진정한 외부 검증은 slide2 (`../../1_152_19/
   proof_ver2/`) 의 결과와 비교하는 것 — 두 슬라이드에서 공통으로 나오는
   pair 가 있다면 그것이 가장 신뢰할 수 있는 신호.

4. **Hist2Cell 모델의 도메인 갭.** Hist2Cell 은 *폐* 조직 단일세포 라벨로
   학습된 모델인데, 본 슬라이드는 *유방* 조직이다. 따라서 양의 cross-
   modality 결합은 "lung 라벨의 분류기 출력이 breast 조직의 구성적 차이
   (상피/면역/기질/혈관) 를 의미있게 따라간다"는 의미이지, **lung cell
   type 이 breast 에서도 같은 정체성을 갖는다는 뜻은 아니다**. 마커 해석
   시에는 라벨명을 "기능 모듈" 정도로 받아들이는 게 안전하다.
