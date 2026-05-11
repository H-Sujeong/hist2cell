# 원본 analysis vs largest-blob x-range filtered analysis

## 필터링 방식

각 v2 슬라이드 spot (X, Y) 좌표에 kNN(k=6) 그래프 → connected components →
가장 큰 blob 의 [Xmin, Xmax] 만 유지하도록 모든 spot 의 X 를 범위 내로 잘랐다.
(Y 는 미제약 — 사용자 요구사항대로 X-range 만.)

| slide | n_orig | n_filtered | kept % | n_components | largest_comp % | xmin | xmax |
|---|---:|---:|---:|---:|---:|---:|---:|
| slide1_085_12 | 35,821 | 21,659 | 60.5% | 12 | 59.9% | 12,600 | 137,400 |
| slide2_152_19 | 40,502 | 27,339 | 67.5% | 26 | 65.6% | 44,600 | 176,600 |

두 슬라이드 모두 **여러 개의 분리된 조직 덩어리** 가 있었음. 가장 큰 덩어리만
60–66% 차지. 나머지 (각각 26%/13% (slide1), 30% (slide2) 등) 가 큰 덩어리의
x-range 바깥에 있어서 필터에서 빠짐. (큰 덩어리 안 spot 외에도 같은 x-range 안에
들어간 일부 작은 component spot 은 살아 남음 — 그래서 kept % 가 largest comp %
보다 약간 큼.)

## 주요 차이점

### 1. cell-type / group mean 의 일률적 상승 (denominator effect)

대부분 cell type 의 mean_per_spot 이 +20~65% 증가. 이건 진짜 신호 증가가
아니라, 필터로 빠진 ~33–40% spot 이 대부분 가장자리/저신호 spot 이라서
**전체 평균을 희석시키던 분모가 작아진 결과**. group level 에서도 거의 모든
group 이 +20–60% 상승.

### 2. slide1 — 패턴 보존, 절대값만 상승

| group | orig mean/spot | filt mean/spot | Δ% |
|---|---:|---:|---:|
| Stromal-muscle | 2.227 | 3.576 | +60.6% |
| Stromal-fibroblast | 1.814 | 2.745 | +51.3% |
| Epithelial-alveolar | 1.458 | 1.897 | +30.1% |
| Immune-lymphoid | 1.245 | 1.717 | +37.9% |
| Broad epithelial-activity proxy | 1.013 | 1.293 | +27.6% |

**그룹 순위 변화 거의 없음** — Stromal-muscle 1위, Stromal-fibroblast 2위 패턴
유지. broad epithelial-activity proxy 비중도 비슷. 즉 slide1 의 "stromal-rich, broad-proxy 중간"
결론은 robust.

### 3. slide2 — Secretory_Goblet 이 -54%, 압축적 Epithelial-airway 우세 약화

| cell type | orig | filt | Δ% |
|---|---:|---:|---:|
| **Secretory_Goblet** | 0.382 | 0.177 | **-53.8%** |
| Ciliated | 1.215 | 1.098 | -9.7% |
| AT2 | 1.102 | 1.415 | +28.4% |
| Fibro_alveolar | 0.812 | 1.120 | +37.9% |

| group | orig | filt | Δ% |
|---|---:|---:|---:|
| Epithelial-airway | 2.706 | 2.591 | **-4.3%** |
| Epithelial-alveolar | 1.802 | 2.379 | +32.0% |
| Broad epithelial-activity proxy | 1.427 | 1.730 | +21.3% |

**slide2 는 큰 변화** — Secretory_Goblet 이 -54%, Ciliated 가 -10% 떨어졌고
Epithelial-airway 그룹 전체가 -4% (다른 그룹은 +30~40% 상승하는 와중). 즉
**가장 큰 덩어리 바깥의 작은 섹션이 mucinous/ciliated airway 에 강한 신호를
주고 있었음**. 원본 분석에서 "Epithelial-airway 1위, mucinous compartment 의
Goblet ↔ immune mutual exclusion" 결론은 그 작은 섹션의 영향이 컸다는 정황.

### 4. Moran R 공간 자기상관: 두 슬라이드 모두 약화

| slide | orig diag mean | filt diag mean | range orig | range filt |
|---|---:|---:|---:|---:|
| slide1 | 0.683 | 0.626 | [0.193, 0.830] | [0.149, 0.770] |
| slide2 | 0.665 | 0.535 | [0.401, 0.843] | [0.284, 0.746] |

좁은 영역으로 잘린 만큼 spatial scale 이 짧아져 R 의 분포가 약간 좁아짐.
slide2 의 감소폭 (0.665→0.535) 이 더 큰 것은 위 §3 의 Goblet/Ciliated 섹션
제거가 큰 spatial 신호를 잃게 했기 때문.

### 5. top cell-pair Moran R 의 community 재배열

slide1:
- orig top: Monocyte_CD16/NKT, Macrophage_intermediate/Monocyte_CD16, B_memory/Monocyte_CD16 (myeloid+lymphoid 혼합)
- filt top: NK_CD16hi/NK_CD11d, Monocyte_CD16/NKT, B_naive/NK_CD16hi (NK 중심 클러스터)
- 같은 immune-co-clustering 테마 유지, 구체 페어는 재편

slide2:
- orig top: B_memory/DC_1, B_memory/Monocyte_CD14, B_memory/Monocyte_CD16 (B-myeloid axis)
- filt top: DC_1/Macro_int, DC_1/Macro_interstitial, B_memory/DC_1 (myeloid 우세)
- B_memory 가 떨어지고 DC/Macro 가 올라옴 → B-cell-rich 영역이 잘려나간 정황

## 결론

| 항목 | slide1 | slide2 |
|---|---|---|
| 그룹 순위 | 변화 거의 없음 | Epithelial-airway 약화, alveolar 강화 |
| broad epithelial-activity proxy 패턴 | 유지 (+27.6%) | 유지 (+21.3%) — 단 broad-dominant 비율은 17.7% → 3.6% 로 큰 감소 |
| Goblet/Ciliated | 변화 없음 | **명확히 감소** (Goblet -54%) |
| Moran R | 약간 감소 (-8%) | 명확히 감소 (-20%) |
| 가장 큰 정성적 변화 | 거의 없음 | mucinous compartment 신호 손실 |

원본 분석의 정성적 결론은 두 슬라이드 모두 대체로 robust (epithelial-activity proxy /
immune-cluster 패턴 유지). 다만 **slide2 의 "Goblet–immune mutual exclusion"
은 큰 덩어리 바깥의 ciliated/mucinous 섹션의 기여가 크다** — 가장 큰 덩어리
하나만 보면 그 신호가 약해진다는 점이 의미 있는 발견.

mean_per_spot 절대값 비교는 denominator 효과 때문에 직접 비교 불가.
**원본 vs 필터 비교는 그룹 ratio / Moran R / top-pair identity** 같은
denominator-invariant 지표로만 해석할 것.
