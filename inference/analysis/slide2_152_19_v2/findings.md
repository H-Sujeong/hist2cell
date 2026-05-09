# slide2_152_19 — 분석 소견 (sanity-check, lung-trained × breast)

> **⚠️ 본 소견은 "모델이 출력한 신호" 의 패턴 기술이지 breast 조직에 대한 임상 진단이 아니다.**
>
> Hist2Cell 가중치는 **healthy human lung** 학습본이고 입력 슬라이드는 KBSMC **breast** SVS 다. 80개 cell type 라벨은 모두 lung 분류이므로 절대값/세부 sub-type 해석은 무의미하고, 다음 두 가지만 의미가 있다.
> 1. **그룹 단위 (immune / epithelial / stromal / vascular …) 의 상대 분포와 spatial pattern**
> 2. **다른 modality (proteomics 등) 와의 spatial co-registration 시 reference signal**
>
> 따라서 아래 소견은 모두 "모델이 이렇게 보고했다" 라는 함의로 읽고, biological 결론은 proteomics/IHC 같은 ground-truth modality 로 검증 필요.

---

## 1. 한 줄 요약

**epithelial-rich, 활발한 immune+proliferative signal 슬라이드.** Epithelial-airway 가 모든 그룹 중 압도적 (μ=2.71/spot), Epithelial-alveolar 가 그 다음. Immune signal 도 **slide1 보다 명확히 강함** (lymphoid 1.64 + myeloid 1.09 = 2.73), proliferative epithelial proxy 도 1.43 으로 slide1 (1.01) 의 1.4배. 공간적으로는 **mucus-secreting epithelial (Goblet) 영역과 immune cell 영역의 강한 배타성**이 두드러진다.

스폿: **40,502** (tile 400×400 px @ mpp=0.2615μm/px → 약 105 μm 격자)

---

## 2. Group composition (top→bottom)

| group | n_celltypes | mean per spot | sum_total | 비고 |
|---|---:|---:|---:|---|
| **Epithelial-airway** | 14 | **2.71** | 109,600 | 압도적 1위 |
| Epithelial-alveolar | 3 | 1.80 | 73,000 | AT1+AT2+Dividing_AT2 |
| Immune-lymphoid | 20 | 1.64 | 66,288 | slide1 보다 31% 강함 |
| Stromal-fibroblast | 6 | 1.49 | 60,265 | |
| **Cancer-proxy** | 5 | **1.43** | 57,796 | slide1 (1.01) 대비 41% 높음 |
| Vascular | 7 | 1.31 | 53,142 | |
| Stromal-muscle | 6 | 1.19 | 48,359 | slide1 (2.23) 대비 절반 |
| Immune-myeloid | 16 | 1.09 | 44,012 | slide1 (0.62) 대비 76% 강함 |
| Stromal-other | 4 | 0.14 | 5,857 | |
| Other-blood | 2 | 0.09 | 3,743 | |
| Neural | 2 | 0.08 | 3,285 | |

→ **epithelial (airway + alveolar) 합산 4.51** 로 모든 그룹/메타-그룹 중 최강. immune 합 2.73 도 매우 높음. stromal 합산은 2.83 로 slide1 (4.22) 보다 낮음.

상위 개별 cell type (top 10): Ciliated (μ=1.22), AT2, Fibro_alveolar, AT1, Endothelia_vascular_Cap_a, Muscle_smooth_syst_arterial, Fibro_adventitial, Endothelia_vascular_Cap_g, Secretory_Goblet, Muscle_airway.

→ Ciliated + Goblet 같은 **mucosal/airway epithelial** 신호가 강하다는 점이 lung-trained 모델로서는 자연스럽지만, breast 조직에서는 **luminal epithelium** 류로 read 되는 것일 가능성.

---

## 3. Immune ↔ cancer-proxy 관계

- **spot-level Pearson (immune_total vs cancer_proxy_total): ρ = 0.816** ← 강한 양의 상관 (slide1 의 0.936 보다는 낮음)
- 82.3% spots: immune > cancer-proxy
- 17.7% spots: cancer-proxy > immune  ← slide1 (10.9%) 의 1.6배

**해석**: slide2 는 cancer-proxy 우세 영역이 slide1 보다 더 큰 비율을 차지한다. 즉 "proliferative epithelial 이 immune 을 압도하는" spot 이 더 많음. 이는 **활발한 epithelial proliferation 영역이 inflammatory infiltrate 와 부분적으로 분리** 됨을 시사. 이 17.7% spot 들의 spatial 위치 (`spatial_immune_vs_cancer.png` 우측 panel 의 hot-spot) 가 후속 검증의 1차 ROI.

slide1 (cancer 10.9%, ρ=0.94) → slide2 (cancer 17.7%, ρ=0.82) 방향성:
- slide2 가 epithelial-proliferative compartment 와 immune compartment 가 더 분리되어 있다.
- 두 슬라이드가 같은 환자라면 이 차이는 sampling 부위/병변 stage 의 차이로 해석 가능.

---

## 4. 공간 구조 (Moran's R 하이라이트)

### 가장 강한 co-localization (R = 양수, 같이 분포)

slide1 과 마찬가지로 **immune cell 들끼리의 cluster** 가 압도적:
- B_memory ↔ DC_1 (R=0.780)
- B_memory ↔ Monocyte_CD14 (0.779)
- B_memory ↔ Monocyte_CD16 (0.778)
- B_memory ↔ CD8_EM_EMRA (0.775)
- DC_1 ↔ Macro_int (0.774)
- B_memory ↔ NKT (0.772)
- DC_1 ↔ Macro_interstitial (0.770)
- CD8_EM_EMRA ↔ NKT (0.769)

→ **B 세포 중심의 immune cluster** (B-cell + DC + Monocyte + CD8 + NK) 가 slide2 의 가장 두드러진 spatial co-occurrence pattern. tertiary lymphoid structure 와 가까운 배치.

### 가장 강한 mutual exclusion (R = 음수)

```
Secretory_Goblet ↔ CD4_naive_CM     R=-0.357
Secretory_Goblet ↔ NKT              R=-0.342
Secretory_Goblet ↔ B_memory         R=-0.340
Secretory_Goblet ↔ CD8_EM_EMRA      R=-0.340
Secretory_Goblet ↔ Monocyte_CD14    R=-0.338
Secretory_Goblet ↔ Fibro_alveolar   R=-0.334
Secretory_Goblet ↔ Erythrocyte      R=-0.329
Secretory_Goblet ↔ Macro_int        R=-0.328
Secretory_Goblet ↔ Monocyte_CD16    R=-0.328
B_plasma_IgG ↔ Schwann_nonmyelinating R=-0.325
```

→ **Secretory_Goblet (mucus-producing epithelial) 이 거의 모든 immune cell 과 강하게 공간 분리**. lung-context 에서는 mucosal goblet epithelium 영역이 immune infiltrate 영역과 다른 anatomical compartment 에 있음을 의미. breast 슬라이드에서 이 패턴이 나오는 건 모델이 **mucin-rich / luminal-secretory 영역**을 Goblet 으로 read 했고 그 영역이 실제로 immune 침윤이 적은 곳일 가능성. **proteomics 의 mucin marker (MUC1/MUC5AC) 와의 매칭 검증 1순위**.

### Cancer-proxy 자기상관 (single-var Moran I = diag R)

| cell type | I |
|---|---:|
| AT2 | 0.682 |
| Dividing_AT2 | 0.629 |
| Dividing_Basal | 0.579 |
| Suprabasal | 0.523 |
| Basal | 0.475 |

→ slide1 보다 전반적으로 단일-cell 자기상관이 약간 낮음 (대각 평균 0.578 vs slide1 0.560). 그러나 모든 5개 type 이 0.5 이상 → 여전히 **공간적으로 응집된 blob 형태**. AT2/Dividing_AT2 가 가장 큰 blob 을 만든다.

---

## 5. 첨부 그림 어떻게 읽나

- `spatial_top10_celltypes.png` — Ciliated 와 AT2 가 중앙 조직 영역 거의 전체에 분포. Secretory_Goblet 은 더 sparse 하고 부분 영역 (mucus-rich pocket) 에 응집.
- `spatial_group_heatmaps.png` — **Epithelial-airway** panel 이 가장 강하고 광범위. **Immune-lymphoid** panel 은 epithelial 영역과 부분적으로 겹치되 더 sparse 한 hot-spot 형태.
- `spatial_immune_vs_cancer.png` — 좌(immune, μ=1.86) 와 우(cancer-proxy, μ=1.01) 모두 중앙 조직에서 활성. immune 쪽이 더 밝고 spread 큼. 우측 panel 의 hot-spot 들이 cancer-proxy ROI.
- `moran_r_clustermap.png` — 좌상단의 큰 immune block (B/T/NK/Macro/Mono/DC) 이 명확. 우하단의 epithelial block 안에서 Goblet 만 separately 어두운 row/column 으로 보일 것 (위 mutual exclusion 패턴 반영).

---

## 6. slide1 과의 비교 (한 환자/연구 단위로 묶어서 보면)

| 지표 | slide1 (085-12) | slide2 (152-19) | 방향 |
|---|---:|---:|---|
| spots | 35,821 | 40,502 | slide2 가 더 큼 |
| Top group | Stromal-muscle (2.23) | **Epithelial-airway (2.71)** | 정반대 |
| Stromal 합 | 4.22 | 2.83 | slide1 ↑ |
| Epithelial 합 | 2.67 | **4.51** | slide2 ↑ |
| Immune 합 | 1.86 | **2.73** | slide2 ↑ |
| Cancer-proxy μ | 1.01 | **1.43** | slide2 ↑ |
| immune↔cancer ρ | 0.94 | 0.82 | slide2 가 분리도 더 높음 |
| cancer-우세 spot % | 10.9% | 17.7% | slide2 ↑ |

→ **두 슬라이드는 같은 modality/조직이지만 패턴이 매우 다르다.** slide1 은 "stromal/quiescent" 표현, slide2 는 "epithelial-rich + active immune + active proliferation". 같은 환자라면 sampling 부위 차이 (lesion vs surrounding) 가능성, 다른 환자라면 baseline 차이.

---

## 7. 한계 및 후속 제안

1. **slide label false positive (~5%)** — 좌측 라벨 sticker 가 일부 조직으로 잡힘 (`spot_view.jpg` 좌측 빨간 영역). proteomics 매칭 전 X<8000 정도 컷 권장.
2. **lung 분포 학습** — Ciliated/Goblet 같은 airway epithelial 라벨은 breast 의 luminal/secretory epithelium 으로 읽힌 결과일 가능성 높음. 동료가 proteomics 의 luminal marker (CK8/18, EpCAM, MUC1) 와 spatial 일치 검증해 주면 본 신호의 의미 확정 가능.
3. **cancer-proxy 17.7% 영역** — 이 영역이 slide2 분석의 핵심 ROI. spatial map 상의 위치를 좌표 (`predictions.csv` 에서 `cancer_proxy > immune_total` 필터) 로 추출해 proteomics 의 proliferation/CK 양성 영역과 매칭 권장.
4. **Goblet vs immune 강한 mutual exclusion** — proteomics 의 mucin (MUC1/MUC5AC) marker 분포를 봤을 때 Goblet-rich 영역이 실제로 immune-poor 인지 검증. 일치하면 본 모델의 "mucinous compartment" detection 이 통계적으로 유의함.
5. **slide1 과 함께 보기** — `inference/analysis/slide1_085_12_v2/findings.md` 와 비교. 차이점이 환자/sample 변이의 spatial 시그니처일 가능성.

---

## 8. proteomics 매칭 첫걸음 (체크리스트, slide1 과 동일 패턴)

- [ ] proteomics 데이터의 좌표계 확인 (같은 SVS / consecutive section / bulk)
- [ ] 같은 SVS 면 `predictions.csv` 의 `(X, Y)` 직접 사용 (μm 변환 시 × 0.2615)
- [ ] **immune marker** (CD3/CD8/CD20/CD68) ↔ **immune_total** spatial 일치 검증
- [ ] **proliferation marker** (Ki67/MCM2/PCNA) ↔ **cancer_proxy_total** 일치 검증
- [ ] **mucin marker** (MUC1/MUC5AC) ↔ **Secretory_Goblet** 일치 검증 (slide2 특이적 추가 검증)
- [ ] 위 검증 통과 후 cell-pair Moran's R 패턴까지 신뢰 확장 가능

---

*분석 입력: `predictions.csv` / `predictions.npy` / `slide2_152_19_coords.h5` (한 디렉토리 위), `cell_type_groups.csv` (analysis 디렉토리)*
*분석 코드: `inference/analysis/analyze.py`*
*전체 caveat 및 사용법: `inference/analysis/README.md`*
