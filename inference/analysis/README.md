# `inference/analysis/` — Hist2Cell predictions 후속 분석 (proteomics 매칭 / immune vs cancer 비교용)

> **⚠️ 핵심 caveat (먼저 읽기)**
>
> 본 디렉토리의 모든 분석은 **healthy human lung 으로 학습된 Hist2Cell 가중치** (`humanlung_cell2location_leave_A50_out.pth`) 를 KBSMC **breast** SVS 두 장에 적용한 결과를 가공한 것이다.
>
> - 80개 cell type 라벨은 모두 **폐 전용** (AT1/AT2/Basal/Ciliated/SMG_*/Schwann/Secretory_*/B,T cell 서브타입 등). breast 조직에 대해 절대값을 그대로 해석하면 안 된다.
> - **cancer cell type 자체가 없음.** 이 디렉토리에서 "cancer-proxy" 라고 부르는 그룹은 **proliferative epithelial cell types (AT2, Basal, Suprabasal, Dividing_AT2, Dividing_Basal)** 의 합으로, "분열 가능한 상피세포 신호" 의 spatial proxy 일 뿐 — breast cancer 직접 예측 아님.
> - **immune cell 카테고리** (B/T/NK/Macro/DC/...) 는 조직 보편적이라 spatial 분포의 상대적 의미는 가질 수 있으나, sub-type 단위 (예: CD4_TRM vs CD8_EM) 는 lung-tissue-resident 표현이라 신뢰성 낮음. **그룹 단위 (Immune-lymphoid / Immune-myeloid)** 합산값으로 보는 것이 안전함.
>
> 사용 목적: **proteomics modality 와 spatial registration**, **cell composition pattern 의 first pass exploratory** — 임상 해석 금지.

---

## 1. 디렉토리 구조

```
inference/analysis/
├── README.md                          ← 이 파일
├── cell_type_groups.csv               80 cell type → lineage group + cancer-proxy flag
├── analyze.py                         재실행 스크립트 (한 번에 다 만듦)
├── slide1_085_12_v2/
│   ├── abundance_by_celltype.csv      80 type 별 mean / median / max / fraction-nonzero (mean 내림차순)
│   ├── abundance_by_group.csv         그룹별 합산 통계 (+ cancer-proxy 별도 row)
│   ├── spatial_top10_celltypes.png    상위 10 type 의 spot-level abundance scatter
│   ├── spatial_group_heatmaps.png     10개 group 의 spatial sum heatmap
│   ├── spatial_immune_vs_cancer.png   immune total vs cancer-proxy total 1×2 panel
│   ├── moran_r_pairs.csv              cell-pair (3,240 rows = C(80,2)+80 diag) bivariate Moran's R
│   └── moran_r_clustermap.png         80×80 Moran's R hierarchical clustermap
└── slide2_152_19_v2/                  (동일 구조)
```

상위 디렉토리 `inference/slide{1,2}_*_v2/` 는 prep + inference 산출물:
```
inference/slide1_085_12_v2/                # (and slide2_152_19_v2)
├── predictions.csv      메인 — spot_id, X, Y + 80 cell type 컬럼 (행 N)
├── predictions.npy      float32 [N, 80] (csv 와 같은 행 순서)
├── slide1_085_12_coords.h5  coords[N,2] + metadata (mpp, dims, tile_size 등)
├── spots.csv            spot_id, X, Y, tile_x_topleft, tile_y_topleft
├── spot_view.jpg        thumbnail 위에 tile 박스 overlay (정확 스케일)
├── tissue_mask.png      tissue mask (thumbnail 해상도)
└── top6_cell_heatmaps.png  (legacy v1 시점 — 무시 가능, analysis/ 의 top10 사용)
```

---

## 2. 인덱스 정렬 규칙 (모든 row order 가 동일)

```
predictions.csv  row i  ↔  predictions.npy[i]  ↔  coords.h5["coords"][i]  ↔  spots.csv row i
```
- N (slide1) = 35,821, N (slide2) = 40,502
- spot_id 형식: `<slide>_x<X>y<Y>` (X, Y 는 tile center 의 level-0 픽셀 좌표)

CSV 와 npy 는 같은 데이터를 다른 형식으로 저장한 것이다. 컬럼 80개 순서는 `cell_type_groups.csv` 의 첫 컬럼 순서와 같지 **않을 수 있음** — predictions.csv 헤더를 신뢰. `cell_type_groups.csv` 는 매핑 테이블이지 순서 정의 아님.

---

## 3. 좌표계와 proteomics 매칭

- `X`, `Y` 는 **slide level-0 (full-res) 픽셀 좌표**. tile center.
- 두 슬라이드 모두 Aperio SVS, **mpp = 0.2615 μm/px**, objective = 40×.
- 물리 좌표 변환:
  ```python
  physical_x_um = X * 0.2615        # μm
  physical_y_um = Y * 0.2615
  physical_x_mm = X * 0.2615 / 1000 # mm
  ```

### proteomics 매칭 워크플로 (3 가지 케이스)

#### Case A) proteomics 가 **같은 SVS 의 mass-spec 이미지** 라서 동일 좌표계
- 그냥 `(X, Y)` 의 픽셀 좌표를 그대로 사용. 본 분석의 `predictions.csv` 와 left-join.

#### Case B) proteomics 가 **다른 슬라이드 (consecutive section 등)** — affine register 필요
- 두 슬라이드의 thumbnail 에서 **공통 fiducial** (조직 가장자리, 표지, blood vessel) 식별
- 4-point affine 추정 (`cv2.getAffineTransform` 또는 `skimage.transform.AffineTransform`)
- 본 결과의 `(X, Y)` 를 affine 적용 → proteomics 좌표계로 변환 후 nearest-neighbor 매칭

#### Case C) **bulk proteomics** (whole-slide-level) — slide 평균과 비교
- `abundance_by_group.csv` 의 `mean_per_spot` 또는 `sum_total` 컬럼을 슬라이드 단위 cell composition 으로 사용
- proteomics 의 protein abundance vector 와 correlation / regression

---

## 4. patch 재추출 (model 입력으로 들어간 224×224 이미지를 그대로 다시 보고 싶을 때)

```python
import openslide, pandas as pd
from PIL import Image

slide_path = "/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-085-12,.svs"
df = pd.read_csv("inference/slide1_085_12_v2/predictions.csv")
sl = openslide.OpenSlide(slide_path)

i = 100                                       # 보고 싶은 spot row
X, Y = int(df.loc[i, "X"]), int(df.loc[i, "Y"])
patch = sl.read_region((X - 112, Y - 112), 0, (224, 224)).convert("RGB")
patch.save(f"spot_{df.loc[i, 'spot_id']}.png")
sl.close()
```

원래 prep 단계에서도 같은 좌표 변환 (`x0 = cx - patch_size//2`) 을 사용했으므로 동일 패치가 나온다. ImageNet 정규화는 모델 입력 직전에 적용된 것이라, 시각화용으로는 raw RGB 가 자연스럽다.

---

## 5. cell type grouping (`cell_type_groups.csv`)

10개 lineage group + 1개 cancer-proxy flag:

| group | n | 멤버 예시 |
|---|---:|---|
| Immune-lymphoid | 20 | B_*, CD4_*, CD8_*, NK_*, T_reg, gdT, ILC, MAIT, NKT |
| Immune-myeloid | 16 | DC_*, Macro_* (8), Monocyte_*, Mast_cell, Macrophage_intermediate |
| Epithelial-airway | 14 | Basal, Ciliated, Goblet/Club, SMG_*, Suprabasal, Myoepithelial, … |
| Vascular | 7 | Endothelia_* (lymphatic + 6 vascular sub-types) |
| Stromal-fibroblast | 6 | Fibro_adventitial/alveolar/myofibroblast/peribronchial/perichondrial/immune_recruiting |
| Stromal-muscle | 6 | Muscle_smooth_*, Muscle_airway, pericyte_* |
| Stromal-other | 4 | Chondrocyte, Mesothelia, NAF_endoneurial, NAF_perineurial |
| Epithelial-alveolar | 3 | AT1, AT2, Dividing_AT2 |
| Neural | 2 | Schwann_Myelinating / nonmyelinating |
| Other-blood | 2 | Erythrocyte, Megakaryocyte |

`is_cancer_proxy=1` (총 5개): **AT2, Basal, Suprabasal, Dividing_AT2, Dividing_Basal** — 분열/재생 가능 epithelial 만 모은 그룹. lung 학습 분포에서 'tumor-like proliferative tissue' 와 가장 가까운 신호로 해석 시도. **breast cancer prediction 아님.**

---

## 6. immune vs cancer-proxy 비교 (핵심 결과)

`spatial_immune_vs_cancer.png` 두 panel 비교를 권장. 또는 `predictions.csv` 직접:

```python
import pandas as pd, numpy as np
df = pd.read_csv("inference/slide1_085_12_v2/predictions.csv")
groups = pd.read_csv("inference/analysis/cell_type_groups.csv")

immune = groups.loc[groups["group"].isin(["Immune-lymphoid", "Immune-myeloid"]), "cell_type"]
cancer = groups.loc[groups["is_cancer_proxy"] == 1, "cell_type"]

df["immune_total"] = df[immune.tolist()].sum(axis=1)
df["cancer_proxy_total"] = df[cancer.tolist()].sum(axis=1)

# 영역별 비교 (예: 4분면)
df["quad"] = (df["X"] > df["X"].median()).astype(int) * 2 + (df["Y"] > df["Y"].median()).astype(int)
print(df.groupby("quad")[["immune_total", "cancer_proxy_total"]].mean())

# spot-level scatter
df.plot.scatter("immune_total", "cancer_proxy_total", s=1, alpha=0.3)
```

### 슬라이드별 요약 (현재 분석 결과)

**slide1_085_12** (35,821 spots):
- Top group: **Stromal-muscle** (μ=2.23) > Stromal-fibroblast (1.81) > Epithelial-alveolar (1.46) > Immune-lymphoid (1.25)
- Cancer-proxy μ=1.01, Immune-myeloid μ=0.62

**slide2_152_19** (40,502 spots):
- Top group: **Epithelial-airway** (μ=2.71) > Epithelial-alveolar (1.80) > Immune-lymphoid (1.64) > Stromal-fibroblast (1.49)
- Cancer-proxy μ=1.43, Immune-myeloid μ=1.09

→ slide2 가 epithelial 신호가 강하고 cancer-proxy 도 더 높음. 두 슬라이드의 패턴이 다름은 모델 출력 분포의 차이일 뿐, biological 결론으로 직결시키지 말 것.

---

## 7. spatial autocorrelation (Moran's R)

`moran_r_pairs.csv` — cell-pair 별 bivariate Moran's R.

- **R > 0**: 두 cell type 이 **공간적으로 같이 분포** (co-localization)
- **R < 0**: **상호 배제** (mutual exclusion)
- **R ≈ 0**: 공간 독립

`moran_r_clustermap.png` 에서 hierarchical clustering 으로 묶인 block 이 곧 "함께 다니는 cell types" 묶음.

slide1 R range: [-0.288, 0.830], 대각선 (single-var Moran's I) 평균 0.683 → 대부분의 cell type 이 spatial autocorrelation 강함 (이웃 spot 이 비슷한 abundance 를 가짐).

읽는 법 예시:
```python
import pandas as pd
mr = pd.read_csv("inference/analysis/slide1_085_12_v2/moran_r_pairs.csv")
# 가장 강하게 co-localize 하는 pair (자기 자신 제외)
mr[mr["A"] != mr["B"]].nlargest(20, "R")
# 가장 강하게 mutual exclusion
mr[mr["A"] != mr["B"]].nsmallest(20, "R")
```

---

## 8. analyze.py 재실행

```bash
# 가상환경: /home/sjhong/hist2cell/.venv (uv 관리)
# 의존성: numpy, pandas, h5py, scipy, seaborn, matplotlib

python inference/analysis/analyze.py \
  --predictions inference/slide1_085_12_v2/predictions.csv \
  --coords      inference/slide1_085_12_v2/slide1_085_12_coords.h5 \
  --groups      inference/analysis/cell_type_groups.csv \
  --output      inference/analysis/slide1_085_12_v2

# slide2 동일 패턴 (--predictions / --coords / --output 만 변경)

# Moran's R kNN 변경 (default k=20)
python inference/analysis/analyze.py ... --knn 30
```

소요 시간: 슬라이드당 약 30–60 초 (CPU). 가장 무거운 부분은 80×80 Moran's R 행렬 계산 (sparse 곱).

---

## 9. 자체 분석 코드 작성용 cheat sheet

```python
import pandas as pd, numpy as np, h5py

# 핵심 4개 입력
preds_df = pd.read_csv("inference/slide1_085_12_v2/predictions.csv")
preds = np.load("inference/slide1_085_12_v2/predictions.npy")           # 같은 데이터 npy 형식
with h5py.File("inference/slide1_085_12_v2/slide1_085_12_coords.h5") as f:
    coords = f["coords"][:]                                              # (N, 2) tile top-left
    meta = dict(f["metadata"].attrs.items())                             # mpp, tile_size, etc.
groups = pd.read_csv("inference/analysis/cell_type_groups.csv")

cell_cols = [c for c in preds_df.columns if c not in ("spot_id", "X", "Y")]
assert len(cell_cols) == 80
assert preds.shape == (len(preds_df), 80)

# 그룹 합산
def group_sum(df, g_name):
    members = groups.loc[groups["group"] == g_name, "cell_type"].tolist()
    return df[members].sum(axis=1)

immune_total = group_sum(preds_df, "Immune-lymphoid") + group_sum(preds_df, "Immune-myeloid")
print("immune total per spot — mean:", immune_total.mean())
```

---

## 10. caveats 한 번 더

1. **lung 학습 → breast 적용**: cell type 이름은 lung 기준. immune 그룹 (보편적), epithelial 그룹 (lung 특이적), stromal 그룹 (대체로 보편적), neural/blood (희귀) 순으로 신뢰도 차등.
2. **cancer-proxy ≠ cancer prediction**: 본 모델은 cancer cell 자체를 예측하지 않는다. proliferative epithelial 합산은 어디까지나 spatial reference signal.
3. **mpp mismatch**: 학습 분포는 ~0.5 μm/px (Visium 20×), 본 슬라이드는 0.2615 μm/px (Aperio 40×). 모델이 보는 224×224 시야가 학습 시야보다 절반 수준이라 결과의 절대값보다 **공간 패턴** 으로 해석.
4. **prep 단계 false positive**: tissue mask 가 슬라이드 라벨/inkstain 일부를 조직으로 잡음 (~5–10%). spatial map 의 "외곽" 신호는 무시. proteomics 매칭 시 tissue 영역만 필터하려면 `spot_view.jpg` 를 보고 X 범위 컷.

---

## 11. 관련 파일

- 본 분석 스크립트: `inference/analysis/analyze.py`
- 분류 테이블: `inference/analysis/cell_type_groups.csv`
- prep 코드: `prep/prepare_wsi_for_inference_v2.py`
- inference 코드: `inference/infer.py`
- 전체 보고서: `report/03_breast슬라이드2장_lung가중치_추론결과_v2framework.md`
- 학습 가중치: `model_weights/humanlung_cell2location_leave_A50_out.pth`
- 학습 분포의 cell type 정의: `example_data/humanlung_cell2location/cell_types.pkl`
