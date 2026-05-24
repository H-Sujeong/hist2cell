# 80×80 cell-cell Moran R 분석 — 방법론과 의의

> **이 문서의 위치**
>
> `inference/analysis_spatial/analyze.py` 의 `build_knn_weight_matrix()` + `moran_r_pairs()` 가 슬라이드별 산출물 `moran_r_pairs.csv` (대각 포함 3,240 row) + `moran_r_clustermap.png` (80×80 hierarchical clustermap) 을 만드는 절차의 수학적 / 생물학적 근거. `findings.md` §3.4 의 해석이 본 문서를 기반으로 한다.

---

## 1. 무엇을 측정하나

본 분석은 슬라이드 한 장의 spot-level Hist2Cell 예측 (`predictions.csv`, N × 80) 에 대해 80 cell type *쌍 (pair)* 의 **공간 공국 (spatial co-localization)** 정도를 정량화한다. 결과는 80×80 matrix R 로:

- **대각 원소** R(x, x) — 단일 cell type 의 공간 자기상관 (univariate Moran's I): 그 type 의 abundance 가 *얼마나 응집되어 있는가*.
- **비대각 원소** R(x, y) — 두 cell type 의 공간 공국 (bivariate Moran's R/I): type x 의 hot-spot 과 type y 의 hot-spot 이 *공간적으로 일치하는가 / 분리되는가*.

R 값의 범위는 row-standardized weight 하에서 대략 [-1, 1]:
- R > 0: positive spatial association — 같은 영역에 공존
- R ≈ 0: spatial 독립 — 무관
- R < 0: negative spatial association — mutual exclusion (anatomical 분리)

이는 단순 abundance correlation (Pearson ρ on spot-by-spot abundance) 과 다르다 — Pearson 은 spot 단위 동시 발현을, Moran's R 은 spot + **그 이웃** 까지 묶어 spatial neighborhood 안의 co-occurrence 를 잡는다.

---

## 2. 수학적 정의

### 2.1 공간 가중 matrix W

n 개 spot 의 (X, Y) 좌표를 입력으로 `cKDTree` 로 각 spot 의 **k=20 nearest neighbors** 을 찾는다 (default `--knn 20`, spot 수에 따라 조정 — `findings.md` 의 §4.4 가이드 참고).

```
W_raw[i, j] = 1   if j ∈ kNN(i)  (i 제외)
            = 0   otherwise
```

→ **symmetrize (union)**: `W = W_raw + W_raw.T` → 이웃 관계가 한쪽만 들어 있는 경우도 양방향으로 잡음.

→ **row-normalize**: 각 row 의 합을 1 로 정규화. 이로 인해 모든 spot 이 동등한 spatial influence 를 갖게 됨 (Cliff & Ord 표준).

```
W[i, j] = W[i, j] / Σⱼ W[i, j]
```

`analyze.py:build_knn_weight_matrix()` 가 이를 sparse `csr_matrix` 로 보존 — N=20-40k spot 에서도 메모리 효율적.

### 2.2 z-score normalization

각 cell type column 을 평균 0, 표준편차 1 로 표준화:

```
Z[i, m] = (preds[i, m] - mean(preds[:, m])) / std(preds[:, m])
```

표준편차 0 인 cell type (전부 0 예측 등) 은 분모를 1 로 두어 NaN 방지.

### 2.3 Bivariate Moran's R (matrix form)

표준 univariate Moran's I 는:

```
I = (n / S₀) · Σᵢⱼ wᵢⱼ (xᵢ - x̄)(xⱼ - x̄) / Σᵢ (xᵢ - x̄)²
```

여기서 S₀ = Σᵢⱼ wᵢⱼ. row-standardized W 의 경우 S₀ = n 이고, z-scored x 의 분산은 1 이라 다음으로 단순화:

```
I = z^T W z / n
```

본 분석은 이를 *2 변수* 로 확장 — **bivariate Moran's R**:

```
R(x, y) = z_x^T W z_y / n
```

80 cell type 에 대해 한번에 계산:

```python
WZ = W @ Z          # (n, 80)  — 각 type 의 spatial lag
R  = Z.T @ WZ / n   # (80, 80)
```

`analyze.py:moran_r_pairs()` 의 핵심 두 줄. sparse W 와 dense Z 의 곱이라 O(n × k × m + n × m²) ≈ N·20·80 + N·80² 으로 35k spot × 80 type 에서 약 1 초.

### 2.4 통계적 유의성 — SE 근사 + p-value

Cliff & Ord 의 **무작위화 (randomization) 가정 하의 분산 근사** 를 사용:

```
Var(R) ≈ 2 · ||W||²_F / n²
SE     = √(Var(R)) = √(2 · Σᵢⱼ wᵢⱼ²) / n
```

`||W||_F` = Frobenius norm of W (row-normalized 후의 weight 들의 제곱합).

```
z-score = R / SE
p-value = 2 · (1 - Φ(|z-score|))   # two-sided normal approx
```

⚠️ 본 SE 는 **단순 근사**. 엄밀한 통계 검증을 원하면 **permutation-based p-value** (각 cell type column 을 random 으로 spatial permutation 한 후 R 의 null distribution 구함) 가 더 보수적. 본 분석은 *탐색용* 으로 sufficient — 정량 검증 단계에서는 permutation 권장.

`moran_r_pairs.csv` 의 컬럼: `A, B, R, z, p` — 80 type × (80+1)/2 = 3,240 row.

---

## 3. 80×80 matrix 의 의미 — 산출물 해석

### 3.1 대각 (diagonal) — univariate Moran's I

각 cell type 의 공간 자기상관:

| R(x, x) 범위 | 해석 |
|---|---|
| 0.7 이상 | 강한 spatial blob — 응집된 hot-spot 형성. 후속 ROI 검증 1순위 |
| 0.4 – 0.7 | 적당한 blob — 부분 응집 |
| 0.2 – 0.4 | 약한 응집 — dispersed |
| ~0 | 무작위 분포 |
| 음수 | checkerboard 패턴 (드물지만 일부 noise label 에서 발생) |

`findings.md` §3.4 의 "cancer-proxy 5 종의 자기상관" 표가 이 대각 값을 인용한다.

### 3.2 비대각 (off-diagonal) — bivariate Moran's R

cell type 쌍의 공간 공국:

| R(x, y) 범위 | 해석 |
|---|---|
| 0.6 이상 | 강한 co-localization — 같은 영역에 동시 발현 (예: B + T + DC 의 TLS 신호) |
| 0.3 – 0.6 | 적당한 공국 — lineage 가족 또는 community |
| ~0 | 공간 독립 |
| -0.2 이하 | mutual exclusion — anatomical compartment 분리 (예: 상피 ↔ 간질) |

R 의 절대값 크기는 W 의 k 와 row-standardization 가정에 의존 — 절대값 자체보다 **상대 순위 + 부호** 가 해석의 핵심.

### 3.3 hierarchical clustermap

80×80 R matrix 를 seaborn `clustermap` (Ward linkage on row distance) 로 그리면:
- **빨간 block** = 같은 community (co-localized cell types) — 보통 lineage group (immune, stromal, epithelial) 과 일치
- **파란 block / row** = 그 cell type 이 거의 모든 다른 type 과 분리됨 — anatomical outlier
- **블록 구조 자체** = 슬라이드의 spatial organization 의 신호

이로써 **clustering-based community detection** 을 별도 알고리즘 없이 시각적으로 수행 가능.

### 3.4 clustermap 읽는 법 — 시각적 가이드

처음 보면 헷갈리는 게 정상. 본 산출물 (`moran_r_clustermap.png`) 의 *모든 요소* 가 어디서 왔고 무엇을 의미하는지 step-by-step.

#### 3.4.1 산출 과정 한눈에

```
Step 1: predictions.csv
        spot 1   AT1=0.5  AT2=0.3  Basal=0.0  ...  (80 cell type)
        spot 2   AT1=0.1  AT2=0.0  Basal=0.0  ...
        ...
        spot N   ...
        → matrix P (N × 80)

Step 2: 각 column z-score
        Z[:, m] = (P[:, m] - mean) / std
        → matrix Z (N × 80, 평균 0 / 표준편차 1)

Step 3: 공간 가중치 W 구성
        cKDTree 로 각 spot 의 k=20 NN → sparse symmetric row-normalized matrix
        → matrix W (N × N, sparse)

Step 4: bivariate Moran's R 계산
        WZ = W @ Z          (각 cell type 의 spatial lag)
        R  = Z.T @ WZ / N   → matrix R (80 × 80)

Step 5: hierarchical clustering 으로 row/col 순서 재배열
        seaborn clustermap 이 Ward linkage 로 자동 처리

Step 6: heatmap 그림 (vlag cmap, vmin=-0.3, vmax=+0.3)
        → 양옆에 dendrogram, 가운데 80×80 색 매트릭스
```

R(i, j) 는 **"cell type i 의 hot-spot 과 cell type j 의 hot-spot 이 공간적으로 얼마나 같이 있느냐"** 의 정량 수치.

#### 3.4.2 가로 / 세로 읽는 법

```
                    column = cell type B
                    ┌─────────────────────────────────┐
              row = │                                 │
              cell  │      (i, j) 의 색 = R(A, B)     │
              type  │                                 │
              A     │                                 │
                    └─────────────────────────────────┘
```

- **row 라벨 (왼쪽 / 오른쪽)** = cell type A
- **column 라벨 (위 / 아래)** = cell type B
- **셀 색깔** = 두 타입의 Moran R 값
  - 🔴 빨강 (R > 0) → 공간적으로 **공국 (co-localize)**, 같은 영역에 둘 다 있음
  - 🔵 파랑 (R < 0) → 공간적으로 **분리 (mutually exclusive)**, 한쪽이 있으면 다른쪽은 없음
  - ⚪ 흰색 (R ≈ 0) → 공간 무관
- **대각선** (i = j): R(A, A) = univariate Moran's I — 단일 type 의 응집 정도. 항상 양수, 보통 빨강 줄
- **대칭성**: R(A, B) = R(B, A) 이므로 대각선 기준 위 / 아래가 거울상
- **80 type 의 순서**: 원본 알파벳 순이 아니라 **hierarchical clustering (Ward) 으로 재정렬**. 비슷한 spatial 패턴 가진 type 들이 인접하게 배치. 양옆의 tree (dendrogram) 가 이 clustering 결과의 hierarchy

#### 3.4.3 시각 패턴 → 의미 표

| 헤어맵에서 보이는 것 | 의미 |
|---|---|
| 대각선의 빨간 점들 (왼위 → 오른아래 줄) | 각 cell type 의 자기 spatial blob (응집) |
| **대각 근처의 빨간 사각형 block** | 같은 community (co-localized cell types) — 보통 lineage group (immune / stromal / epithelial) |
| **빨간 사각형 안에 더 진한 작은 block** | sub-community (예: immune 안의 B-Mono-DC TLS-like) |
| **파란 row 또는 column 줄** | 그 cell type 이 거의 모든 다른 type 과 분리 — anatomical outlier |
| 파란 사각형 block | 두 community 가 서로 분리 — compartment 경계 (예: epithelial ↔ stromal) |
| 색 거의 없음 (흰색 영역) | spatial 독립 — 무관한 type 들 |

#### 3.4.4 구체 예 — slide2 의 clustermap 읽기

`inference/analysis_spatial/slide2_152_19_v2/moran_r_clustermap.png` 에서:

1. **`Secretory_Goblet` 라벨 찾기** (왼쪽 row 또는 위 column)
2. **그 row 를 따라 가로로 스캔** — Goblet 의 모든 column 과의 R 값 보기
3. **B_memory / CD4_naive_CM / NKT / CD8_EM_EMRA 등 immune column 들의 자리가 모두 파랑** 으로 보임 → 이 부분이 우리가 `findings.md` 에서 말한 "Secretory_Goblet ↔ immune mutual exclusion top 5" 의 시각적 표현
4. **immune type 들 자체끼리는** clustermap 가운데 어딘가의 큰 **빨간 block** — B / T / Mono / DC 가 서로 co-localized = TLS-like 응집

같은 방식으로 slide1 에서는:
- `Deuterosomal` row 를 따라가면 Muscle / Fibro / Vascular column 자리가 파랑 → "상피 (Deuterosomal) compartment 와 stromal 분리"
- immune type 들 끼리 가운데 큰 빨간 block

#### 3.4.5 색 강도의 한계

`vmin=-0.3, vmax=+0.3` 으로 색 범위가 **클램프** 되어 있음 — 즉 R > 0.3 인 강한 양의 값은 모두 같은 진한 빨강으로 표시. 절대값 차이는 색만으론 못 봄. **정확한 값은 `moran_r_pairs.csv` 의 R 컬럼 확인** (또는 top 5 표를 findings.md 가 인용).

이 클램프는 의도된 것 — 극단값이 색 스케일을 잡아먹어 중간 영역이 안 보이는 현상 방지. 대신 절대 정량은 CSV 로 가능하도록 분리.

---

## 4. 생물학적 / 분석적 의의

### 4.1 Hist2Cell × proteomics 매칭에 주는 가치

본 분석은 cell type abundance 의 **mean / max / fraction-nonzero** (`abundance_by_celltype.csv`) 만으로는 보이지 않는 **공간 조직 (spatial organization)** 신호를 잡는다. 같은 mean 을 가진 두 cell type 이라도 한쪽은 hot-spot 응집, 다른 쪽은 dispersed 일 수 있는데, 본 R 의 diagonal 이 이를 구분.

proteomics 매칭 관점에서:
- proteomics 의 high-risk Tumor 마커 (e.g., KIF20A/22/INCENP for slide1) 의 *spatial 위치* 는 Hist2Cell 의 **proliferation-like signal blob (Dividing_AT2 의 high I)** 의 위치와 후속 검증 가능.
- proteomics 의 immune-mixed 영역 (e.g., slide2 의 GZMH/LCK) 은 Hist2Cell 의 **myeloid community** (DC/Macro top R block) 와 spatial overlap 검증 가능.
- 그 외 *공간 안 분리* (Goblet ↔ immune 같은 mutual exclusion) 는 단일-marker proteomics 로는 잡기 어렵고 본 분석이 추가 정보 제공.

### 4.2 cell community detection

본 분석으로 추출한 *community*:

| 슬라이드 | top community | 해석 |
|---|---|---|
| slide1 (원본) | Monocyte_CD16 ↔ NKT ↔ Macrophage_intermediate ↔ B_memory (R ≈ 0.80) | TLS-like immune cluster |
| slide1 (필터) | NK_CD16hi ↔ NK_CD11d ↔ B_naive ↔ Monocyte_CD16 (R ≈ 0.77) | NK-편향 immune 응집 |
| slide2 (원본) | B_memory ↔ DC_1 ↔ Monocyte_CD14/CD16 ↔ CD8_EM_EMRA (R ≈ 0.78) | B-cell 중심 TLS |
| slide2 (필터) | DC_1 ↔ Macro_int / Macro_interstitial / Macro_CCL (R ≈ 0.62) | myeloid 중심 (B-cell 측부 의존) |

→ filter 적용 전후로 community 가 재구성되는 점이 "측부 덩어리에 특정 immune subset 이 응집되어 있었음" 을 정량적으로 드러내는 evidence. 단일 abundance 비교만으론 잡기 어려움.

### 4.3 mutual exclusion 으로 본 compartment 구조

- slide1 원본: Deuterosomal ↔ Muscle/Fibro/Vascular 음수 R 줄 → **상피 (Deuterosomal) compartment 와 간질/혈관 compartment 가 anatomical 으로 분리**
- slide2 원본: Secretory_Goblet ↔ 다수 immune cell (top 5 negative pair 전부 Goblet 관련) → **mucinous airway/ductal-glandular 영역과 immune 영역의 강한 분리**
- slide2 필터: Goblet 관련 음수 R 이 top 5 에서 사라짐 → Goblet 신호가 측부 덩어리 의존이라는 점 직접 입증

→ 단순 abundance 로는 잡을 수 없는 **공간적 격리** 가설을 R < 0 로 정량 가설화.

---

## 5. 한계 (caveats)

### 5.1 방법론적 한계

1. **k=20 kNN 선택의 영향**: k 가 작으면 (≤8) 짧은 spatial scale 만 잡아 *local* noise 증가. k 가 너무 크면 (≥40) distant spot 까지 평균하여 R 이 평탄. `--knn` 가이드 (cookbook §4.4) 의 spot 수 기반 선정 사용.
2. **row-standardization 가정**: 각 spot 의 spatial influence 가 동등하다는 가정. 격자가 불균일하면 (예: WSI 가장자리 spot 의 이웃 수 부족) bias 발생 가능. 본 분석은 `prep_v2` 의 균일 격자라 영향 미미하지만 다른 prep 에선 확인 필요.
3. **SE 근사의 보수성**: 본 분석의 Cliff-Ord 근사는 *대략적* — 엄밀 검증에는 permutation p-value 권장. 본 분석의 z/p 값은 *탐색 단계용*.
4. **다중비교 보정 부재**: 3,240 pair 의 p-value 가 raw — Bonferroni / BH 보정 시 살아남는 유의 페어 수 급감. 본 분석은 *순위 / 부호* 기반 해석에 집중.

### 5.2 데이터 측면 한계

5. **lung-trained 모델 출력**: Hist2Cell 예측 자체가 lung 라벨이므로 R 값은 *모델 출력의 spatial 구조* 이지 실제 cell type spatial 구조의 직접 측정이 아님. lung-derived label 합 (epithelial-activity proxy) 의 cross-tissue 해석은 `EPITHELIAL_PROXY_METHODOLOGY.md` 참조.
6. **error correlation 인플레이션**: 모델이 비슷한 morphology 에서 비슷한 prediction 을 내면 그 type 들 간 R 이 *모델 일반화 능력의 한계* 때문에 인플레이션 될 수 있음 — 즉 진짜 생물학적 co-localization 외에 모델 confusion 도 R 에 기여.
7. **slide 내 heterogeneity 고려 안 함**: R 은 슬라이드 전체 평균. compartment 별 (예: 큰 덩어리 vs 측부) sub-analysis 는 별도 (필터 분석 비교 시 노출됨).
8. **n=1 슬라이드 단위 통계**: 본 R 은 슬라이드 한 장 안의 spatial 구조. cross-slide / cross-patient generalization 은 별도 검증 (현재 n=2 환자로 불가).

---

## 6. 본 분석의 슬라이드 1 / 2 결과 요약

### slide1

- 대각 (응집): AT2 0.745, Dividing_AT2 0.749, Dividing_Basal 0.691 — 강한 blob 형성 (원본).
- top positive: Monocyte_CD16 / NKT / Macrophage_intermediate / B_memory — TLS-like 응집.
- top negative: Deuterosomal ↔ stromal/vascular — 상피 ↔ 간질 anatomical 분리.
- 필터 후: 패턴 보존, NK 편향으로 community 재구성.

### slide2

- 대각: AT2 0.682, Dividing_AT2 0.629 — slide1 보다 약간 약하나 여전히 blob.
- top positive: B_memory ↔ DC_1 / Monocyte / CD8_EM_EMRA — B-cell 중심 TLS.
- top negative: Secretory_Goblet ↔ immune — mucinous compartment 의 명확한 분리.
- 필터 후: B_memory 중심 ↔ myeloid (DC/Macro) 중심으로 community 재구성, **Goblet ↔ immune 음수 페어 모두 top 5 에서 사라짐** (Goblet 자체가 측부 덩어리 의존).

→ 두 슬라이드 모두 *공간 조직* 신호가 정량화되어 후속 ROI proteomics 매칭의 검증 가설을 제공.

---

## 7. 다음 단계 — 정량 검증으로의 확장

ROI 좌표 (`.tmpprotocol`) 도착 시:

1. **ROI 안 spot 들의 평균 R 비교**: 각 ROI tube 내 Hist2Cell spot 들의 cell-type-pair R 평균이 proteomics 의 marker-pair 상관과 매칭되는지.
2. **permutation p-value**: high-rank R pair (e.g., B_memory ↔ DC_1 in slide2) 를 spatial permutation test 로 정밀 검증.
3. **두 modality joint factor**: CCA / MOFA 로 R matrix 와 ROI proteomics matrix 간의 shared latent axis 추출.
4. **CUCA her2st 도착 후 cross-validation**: lung-trained R 의 top community vs breast-trained R 의 top community 의 spatial overlap — proxy 해석의 사후 검증.

---

## 8. 참고문헌

[1] **Moran PAP**. *Notes on continuous stochastic phenomena*. **Biometrika** 1950;37(1/2):17–23. doi:10.2307/2332142 — 원본 Moran's I 정의.

[2] **Cliff AD**, Ord JK. *Spatial Autocorrelation*. London: Pion, 1973. — randomization assumption 하 분산 근사 + 통계 추론 표준.

[3] **Wartenberg D**. *Multivariate spatial correlation: a method for exploratory geographical analysis*. **Geogr Anal** 1985;17(4):263–283. doi:10.1111/j.1538-4632.1985.tb00849.x — bivariate Moran's I 의 multivariate 확장 원본.

[4] **Anselin L**. *Local indicators of spatial association — LISA*. **Geogr Anal** 1995;27(2):93–115. doi:10.1111/j.1538-4632.1995.tb00338.x — local Moran's I (본 분석은 global 만 사용, 후속 local 확장 가능).

[5] **Pebesma E**. *Multivariable geostatistics in S: the gstat package*. **Computers & Geosciences** 2004;30:683–691. — sparse W 기반 large-scale Moran 구현의 reference.

[6] 본 분석 코드: `inference/analysis_spatial/analyze.py` 의 `build_knn_weight_matrix()` + `moran_r_pairs()`.

---

## 9. 본 문서와 연결된 파일

- 분석 코드 / W 구성: `inference/analysis_spatial/analyze.py:build_knn_weight_matrix()`
- R 계산: `inference/analysis_spatial/analyze.py:moran_r_pairs()`
- 산출 CSV / PNG: `inference/analysis_spatial/1_{085_12,152_19}/cell_typing/moran_within_roi.csv` + `moran_r_clustermap.png`
- 슬라이드별 해석: `1_{085_12,152_19}/findings.md` §3.4
- 본 분석의 lung-→breast proxy 한계 (R 해석의 caveat): `EPITHELIAL_PROXY_METHODOLOGY.md`
- 격자 + kNN guideline: `report/04_WSI에서_분석까지_쿡북.md` §4.4
