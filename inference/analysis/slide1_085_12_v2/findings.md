# slide1_085_12 — 분석 소견 (sanity-check, lung-trained × breast)

> **⚠️ 본 소견은 "모델이 출력한 신호" 의 패턴 기술이지 breast 조직에 대한 임상 진단이 아니다.**
>
> Hist2Cell 가중치는 **healthy human lung** 학습본이고 입력 슬라이드는 KBSMC **breast** SVS 다. 80개 cell type 라벨은 모두 lung 분류이므로 절대값/세부 sub-type 해석은 무의미하고, 다음 두 가지만 의미가 있다.
> 1. **그룹 단위 (immune / epithelial / stromal / vascular …) 의 상대 분포와 spatial pattern**
> 2. **다른 modality (proteomics 등) 와의 spatial co-registration 시 reference signal**
>
> 따라서 아래 소견은 모두 "모델이 이렇게 보고했다" 라는 함의로 읽고, biological 결론은 proteomics/IHC 같은 ground-truth modality 로 검증 필요.

---

## 1. 한 줄 요약

**stromal-rich, 비교적 quiescent 한 슬라이드.** Stromal-muscle 계열이 단일 그룹 중 가장 강한 신호 (μ=2.23/spot), epithelial-airway/alveolar 가 그 다음. immune signal 은 중간 수준 (lymphoid 1.25 + myeloid 0.62), proliferative epithelial proxy 도 중간 (μ=1.01). 공간적으로는 **epithelial 영역과 stromal 영역의 분리**가 가장 두드러진다.

스폿: **35,821** (tile 400×400 px @ mpp=0.2615μm/px → 약 105 μm 격자)

---

## 2. Group composition (top→bottom)

| group | n_celltypes | mean per spot | sum_total | 비고 |
|---|---:|---:|---:|---|
| **Stromal-muscle** | 6 | **2.23** | 79,774 | 가장 강한 그룹 |
| Stromal-fibroblast | 6 | 1.81 | 64,971 | |
| Epithelial-alveolar | 3 | 1.46 | 52,245 | AT1+AT2+Dividing_AT2 |
| Immune-lymphoid | 20 | 1.25 | 44,604 | T/B/NK 류 |
| Epithelial-airway | 14 | 1.22 | 43,525 | Ciliated/Goblet/SMG/… |
| Vascular | 7 | 1.20 | 42,997 | 내피 (capillary 우세) |
| **Cancer-proxy** | 5 | **1.01** | 36,293 | AT2/Basal/Suprabasal/Dividing_* |
| Immune-myeloid | 16 | 0.62 | 22,158 | Macro/DC/Mono/Mast |
| Stromal-other | 4 | 0.18 | 6,382 | Chondrocyte/Mesothelia/NAF |
| Neural | 2 | 0.12 | 4,376 | Schwann |
| Other-blood | 2 | 0.07 | 2,368 | Erythrocyte/Megakaryocyte |

→ **stromal (muscle + fibroblast + other) 가 그룹 합산으로 가장 큼** (μ=4.22). 이어 epithelial (airway+alveolar, 2.67), immune (lymphoid+myeloid, 1.86), vascular (1.20) 순.

상위 개별 cell type (top 10): Muscle_smooth_syst_arterial (μ=0.96), AT2, Fibro_adventitial, Fibro_alveolar, AT1, Muscle_airway, Muscle_smooth_pulmonary, Endothelia_vascular_Cap_a, Fibro_myofibroblast, Ciliated.

---

## 3. Immune ↔ cancer-proxy 관계

- **spot-level Pearson (immune_total vs cancer_proxy_total): ρ = 0.936** ← 매우 강한 양의 상관
- 89.1% spots: immune > cancer-proxy
- 10.9% spots: cancer-proxy > immune

**해석**: 모델이 immune-rich 라고 본 영역과 proliferative-epithelial 이 강한 영역이 거의 겹친다. 이는 breast cancer 의 일반적 관찰 (tumor + tumor-infiltrating lymphocyte/myeloid co-occurrence) 과 정성적으로 부합하지만, 본 모델은 lung 분포에서도 같은 패턴을 학습했을 가능성 높음 (alveolar epithelium 주변 macrophage clustering 등). **proteomics 매칭 시 immune-marker 와 proliferation-marker (Ki67 류) 의 공동 분포 검증 권장.**

cancer-proxy 우세 영역 (10.9%) 은 spatial map (`spatial_immune_vs_cancer.png` 우측) 에서 따로 추출 가능 — 후속 검증의 우선 ROI.

---

## 4. 공간 구조 (Moran's R 하이라이트)

### 가장 강한 co-localization (R = 양수, 같이 분포)

immune cell 들끼리의 cluster 가 압도적:
- Monocyte_CD16 ↔ NKT (R=0.802)
- Macrophage_intermediate ↔ Monocyte_CD16 (0.801)
- B_memory ↔ Monocyte_CD16 (0.798)
- Macrophage_intermediate ↔ NKT (0.794)
- DC_1, MAIT, gdT, CD8_EM_EMRA 등도 같은 그룹

→ **classical "immune cell cluster"** 가 슬라이드 내 일정한 영역에 모여 있음. tertiary lymphoid structure 또는 inflammatory infiltrate 의 spatial proxy 가능.

### 가장 강한 mutual exclusion (R = 음수)

```
Deuterosomal ↔ Muscle_airway              R=-0.288
Deuterosomal ↔ Muscle_smooth_pulmonary    R=-0.287
Deuterosomal ↔ Fibro_myofibroblast        R=-0.281
Deuterosomal ↔ Muscle_smooth_syst_arterial R=-0.278
Deuterosomal ↔ Endothelia_vascular_venous_systemic R=-0.276
Deuterosomal ↔ Fibro_peribronchial        R=-0.270
```

→ **상피 (Deuterosomal: 분화 중인 ciliated 전구) 와 stromal/vascular 가 공간적으로 배타적**. 즉 epithelial compartment 와 stromal compartment 가 슬라이드에서 분리된 영역으로 존재. 조직학적 ductal/lobular vs interstitium 구분과 부합 가능.

### Cancer-proxy 자기상관 (single-var Moran I = diag R)

| cell type | I |
|---|---:|
| Dividing_AT2 | 0.749 |
| AT2 | 0.745 |
| Dividing_Basal | 0.691 |
| Suprabasal | 0.333 |
| Basal | 0.280 |

→ AT2/Dividing_AT2 는 **공간적으로 응집된 큰 blob 형태**, 반면 Basal/Suprabasal 은 더 분산되어 있음. proteomics 와 매칭 시 AT2-rich blob 영역이 가장 안정적인 cancer-proxy ROI.

---

## 5. 첨부 그림 어떻게 읽나

- `spatial_top10_celltypes.png` — 평균 abundance 상위 10 type 의 spot scatter. 좌상단 panel (Muscle_smooth_syst_arterial) 의 외곽 strip (좌우 가장자리) 신호는 **slide ink-stain 잔존 false positive** — 무시. 중앙 사각 조직 영역만 해석.
- `spatial_group_heatmaps.png` — 10 그룹 panel. **Stromal-muscle / fibroblast** panel 이 중앙 조직 전반에 강하게 깔린 모습, **Epithelial-airway** 는 더 sparse. **Immune-lymphoid** panel 이 부분 영역에 hot-spot 형성.
- `spatial_immune_vs_cancer.png` — 좌(immune) vs 우(cancer-proxy). 둘이 visually 비슷한 영역에서 강하지만, **immune** 쪽이 절대값과 spread 가 더 큼 (max 14 vs 7).
- `moran_r_clustermap.png` — 좌상단 큰 빨간 block 이 immune-cell cluster (B/T/NK/Macro/Mono/DC). 우하단 stromal+epithelial block 도 분리되어 보임.

---

## 6. 한계 및 후속 제안

1. **slide ink/label false positive (~10%)** — 좌우 가장자리의 vertical strip 은 조직 아님. proteomics 매칭 전 X 좌표 컷 필요 (예: `8000 < X < 200000` 정도, `spot_view.jpg` 보고 결정).
2. **lung 분포 학습** — 위 모든 cell type 라벨은 폐 기준. breast 조직에서:
   - Immune 라벨 (B/T/NK/Macro/DC/Mono) → 보편적, 그룹 합 신뢰 OK
   - Stromal-fibroblast 라벨 → 보편적, 그룹 합 OK
   - **Stromal-muscle** 라벨 → smooth muscle 은 두 조직에 모두 있으나 air-way muscle 등 lung 특이적 sub-type 은 의미 없음. 그룹 합으로만 사용
   - **Epithelial 라벨** → AT1/AT2/Ciliated/Goblet 은 lung 전용. breast 의 ductal/luminal/basal 과 직접 매핑 불가. **AT2/Basal/Suprabasal 합산 cancer-proxy** 만 "분열 가능 epithelial" 의 spatial proxy 로 의미 가짐
3. **AT2/Dividing_AT2-rich blob** — `spatial_immune_vs_cancer.png` 우측 panel 에서 식별되는 hot-spot 들이 후속 validation 의 1차 ROI 후보. proteomics 의 Ki67/proliferation marker, EpCAM/CK 류 신호와 spatial overlap 검증 권장.
4. **slide2 와 비교** — slide2 는 epithelial-rich + 더 강한 immune+cancer-proxy. 두 슬라이드는 같은 patient 의 다른 부위(또는 다른 patient)일 가능성. proteomics 매칭 시 두 슬라이드를 따로 처리 후 비교.

---

## 7. proteomics 매칭 첫걸음 (체크리스트)

- [ ] proteomics 데이터의 좌표계 확인 (같은 SVS / consecutive section / bulk)
- [ ] 같은 SVS 면 `predictions.csv` 의 `(X, Y)` 직접 사용 (μm 변환 시 × 0.2615)
- [ ] consecutive section 면 thumbnail-level affine register 후 좌표 변환
- [ ] proteomics 의 immune marker (CD3/CD8/CD20/CD68 등) 와 본 분석의 **immune_total** spatial 분포 일치 여부 1차 검증
- [ ] proteomics 의 proliferation marker (Ki67, MCM2, PCNA) 와 본 분석의 **cancer_proxy_total** 일치 여부 검증
- [ ] 위 두 검증 통과 시 individual cell-pair Moran's R 패턴까지 신뢰 확장 가능

---

*분석 입력: `predictions.csv` / `predictions.npy` / `slide1_085_12_coords.h5` (한 디렉토리 위), `cell_type_groups.csv` (analysis 디렉토리)*
*분석 코드: `inference/analysis/analyze.py`*
*전체 caveat 및 사용법: `inference/analysis/README.md`*
