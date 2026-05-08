# 원본 Hist2Cell 레퍼런스의 WSI 전처리 파이프라인 — `data_preparation_tutorial.ipynb` 이전 단계 요약

## 1. 이 문서의 목적

`data_preparation_tutorial.ipynb` 는 이미 spot 단위로 잘려진 224×224 patch 와 spot×gene/celltype CSV 가 준비된 상태(`example_data/example_raw_data/<slide>/`) 에서 시작한다. 즉, **WSI 원본 → spot patch + gene expression + cell type ratio** 까지의 단계는 노트북에서 다루지 않고 외부 파이프라인을 가정하고 있다.

이 문서는 README, 노트북의 입력 명세, 그리고 노트북이 인용한 외부 도구를 근거로 **원본 reference 가 WSI 를 어떤 단계로 처리해서 노트북의 입력 형태로 만들었는지** 를 정리한다.

---

## 2. 원본 데이터 출처 및 플랫폼

| 항목 | 내용 |
|---|---|
| 대표 dataset | Healthy Human Lung (Madissoon et al., cell2location 논문 데이터셋) |
| ST 플랫폼 | 10x Genomics **Visium** — hexagonal spot grid, spot 1개당 약 6 이웃 |
| 그 외 dataset | HER2ST (breast cancer), STNet, TCGA, HEST-1k — README §"Datasets and Resources" 참조 |
| 슬라이드 단위 | `WSA_LngSP<...>` 형태의 slide ID, 각 slide 가 하나의 WSI + 하나의 spot grid |

원본 raw 형태는 **Visium SpaceRanger 출력** (WSI `.tiff/.jpg`, `tissue_positions_list.csv`, filtered count matrix `.h5`/`.h5ad`) 에 해당한다고 추정된다. 노트북의 `sp.X_norm5e4_log1p.h5ad` 가 그 통합본이다.

---

## 3. 노트북이 가정하는 입력 (= 원본 파이프라인이 만들어내야 하는 산출물)

`./example_data/example_raw_data/WSA_LngSP9258467/` 안에 있는 파일들이 그대로 **원본 파이프라인의 출력 사양**이다.

### 3.1 이미지 계열
- `WSA_LngSP9258467.jpg` — 원본 WSI (full resolution)
- `WSA_LngSP9258467_low_res.jpg` — 시각화/검수용 다운샘플 이미지
- `spot_view.jpg`, `2x_spot_view.jpg` — spot 위치를 overlay 한 그림 (QC 용)
- `patches/<slide>_<barcode>.jpg` — **spot 좌표 중심 224×224 RGB crop**, spot 수만큼 (예: 422장)
- `2x_patches/...` — super-resolution 용 2× 격자 patch

### 3.2 좌표 계열
- `spots.csv` — `spot_id, X, Y` (WSI **픽셀 좌표**, full-res 기준)
- `2x_spots.csv` — 2× 격자 픽셀 좌표
- `sp.X_norm5e4_log1p.h5ad` — `obs.array_col / obs.array_row` (Visium **격자 좌표**) 포함

### 3.3 발현 계열
- `stdata.csv` — raw count, spot × gene
- `stdata_log1p.csv` — library-normalized + log1p 변환 (full gene set)
- `high_250_stdata.csv`, `high_250_stdata_log1p.csv` — 상위 250개 highly expressed gene 만 추린 버전 (모델 학습 라벨로 사용)

### 3.4 Cell type 계열
- `cell_ratio.csv` — spot × **80 fine-grained cell type abundance ratio**. **cell2location** Bayesian deconvolution 결과를 정규화한 값으로 사용된다. (Visium spot 1개는 보통 1–10 세포를 포함하므로 deconvolution 이 필수)

---

## 4. 원본 WSI → 위 산출물까지의 단계 (재구성)

노트북은 아래 단계를 모두 사전 처리 완료된 것으로 전제한다.

### Step A. WSI 정렬 / 전처리
- SpaceRanger 가 fiducial alignment 후 산출한 `tissue_positions_list.csv` 를 사용해 spot barcode ↔ (array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres) 매핑 확보.
- WSI 자체는 큰 jpg/tiff 한 장으로 보관 (예시에서는 `.jpg`).

### Step B. Spot 단위 patch 추출  ←  **DSMIL-WSI 파이프라인 차용**
README 의 "Troubleshooting" 섹션이 **patch 추출 도구로 [DSMIL-WSI](https://github.com/binli123/dsmil-wsi) 를 명시**한다.
- 각 spot 의 픽셀 좌표 `(X, Y)` 를 중심으로 **224×224 px crop** 을 떠서 `patches/<slide>_<barcode>.jpg` 로 저장.
- in-tissue 필터 (`tissue_positions_list.csv` 의 `in_tissue==1`) 통과한 spot 만 처리.
- ResNet18(ImageNet) 입력에 맞추기 위해 size 224, 최종 학습 시점에서 `transforms.Normalize(mean=[0.485,0.456,0.406], std=[0.229,0.224,0.225])` 적용 (노트북 Step 2).

### Step C. Gene expression 정규화
- raw count → `stdata.csv`
- **library size 5e4 normalization → log1p** (파일명 `norm5e4_log1p` 가 명시) → `stdata_log1p.csv`, AnnData `sp.X_norm5e4_log1p.h5ad`
- 슬라이드 전체에서 발현 상위 250 gene 선정 → `high_250_stdata*.csv`. 이 250개가 모델 라벨의 앞 250 차원이 된다.

### Step D. Cell type deconvolution
- 외부 scRNA-seq reference + spatial count → **cell2location** 으로 spot 별 80개 fine-grained cell type abundance 추정.
- 행 단위로 합 1 이 되도록 비율화하여 `cell_ratio.csv` 생성. 이 80개가 모델 라벨의 뒤 80 차원이 된다.

### Step E. Spot ID 표준화
- 모든 표/이미지의 키를 `<slide>_<barcode>` (예: `WSA_LngSP9258467_AAACCGTTCGTCCAGG-1`) 형태로 통일.
- patch 파일명, `spots.csv` index, `cell_ratio.csv` index, `high_250_stdata*.csv` index, AnnData `obs_names` 가 전부 이 키로 일치해야 노트북 Step 1–4 가 작동.

---

## 5. 노트북이 추가로 수행하는 단계 (이미 처리된 데이터로부터)

원본 파이프라인 산출물이 준비되면, 노트북은 다음만 수행한다 — 이 단계는 외부 도구가 아니라 노트북 자체에서 처리한다.

1. **STDataset / DataLoader** 로 patch 읽고 `transforms.Normalize` 적용 → `[N,3,224,224]` 텐서 (Cell 14–18).
2. `high_250_stdata.csv` + `cell_ratio.csv` merge → `[N, 330]` 라벨 (Cell 14, 17).
3. `sp.X_norm5e4_log1p.h5ad` 의 `array_col/array_row` 를 spot ID 순서대로 정렬 (Cell 22–25).
4. **헥사고날 인접 그래프**: `|Δcol| < 3 and |Δrow| < 2` 를 만족하면 edge, self-loop 포함 → adjacency matrix (Cell 26).
5. `spots.csv` 의 픽셀 좌표를 spot ID 순서대로 정렬 → `pos` (Cell 28).
6. `dense_to_sparse` 로 edge_index 생성 → `torch_geometric.data.Data(x, edge_index, y, pos, spot_id)` 로 묶어 `WSA_LngSP9258467.pt` 로 저장 (Cell 30–31).

따라서 노트북은 **PyG 포맷팅 + 그래프 구성** 단만 책임지고, **WSI 디코딩/등록/패치 추출/발현 정규화/세포 디컨볼루션은 모두 사전 단계**다.

---

## 6. 자체 데이터로 재현하려면 추가로 구현해야 할 것

| 단계 | 도구 / 방법 | 산출 파일 |
|---|---|---|
| WSI + spot 좌표 등록 | SpaceRanger 출력 또는 동등한 Visium pipeline | `tissue_positions_list.csv` 류 |
| Patch 추출 (224×224) | DSMIL-WSI patcher, 또는 OpenSlide + `(X,Y)` 중심 crop | `patches/<spot>.jpg`, `spots.csv` |
| Count → 정규화 | `scanpy.pp.normalize_total(target_sum=5e4)` + `sc.pp.log1p` | `stdata_log1p.csv`, `*.h5ad` |
| 상위 250 gene 선정 | 슬라이드 또는 dataset 단위 mean expression 상위 250 | `high_250_stdata*.csv` |
| Cell deconvolution | **cell2location** + scRNA-seq reference (e.g. Madissoon Human Lung) | `cell_ratio.csv` (spot × 80) |
| Spot ID 일치 검증 | patch 파일명 ↔ csv index ↔ AnnData obs_names | — |

이 표가 채워져야 비로소 `data_preparation_tutorial.ipynb` 가 그대로 돌아간다.

---

## 7. 빠른 체크리스트 (자체 WSI 적용 시)

- [ ] WSI 1장당 in-tissue spot 의 `(barcode, array_row, array_col, pixel_X, pixel_Y)` 가 모두 확보되었는가
- [ ] 모든 spot 에 대해 224×224 patch 가 `<slide>_<barcode>.jpg` 로 저장되어 있는가
- [ ] count matrix 가 `norm5e4 → log1p` 로 변환되어 있고 상위 250 gene 이 추려졌는가
- [ ] cell2location 으로 80 cell type abundance 가 추정되어 row-sum 1 의 ratio 로 저장되어 있는가
- [ ] 위 모든 파일의 spot key 가 `<slide>_<barcode>` 로 통일되어 있는가
- [ ] hexagonal neighbor 임계값(`|Δcol|<3, |Δrow|<2`) 이 사용 중인 ST 플랫폼에 맞는가 (Visium 가 아니면 재정의 필요)

---

## 8. 참고

- 본 repo README §"Project Structure", §"Datasets and Resources", §"Troubleshooting"
- 노트북 `data_preparation_tutorial.ipynb` Cell 12 (raw data structure 표), Cell 22 (AnnData 사용), Cell 24 (hexagonal grid 가정)
- 외부: DSMIL-WSI (https://github.com/binli123/dsmil-wsi), cell2location, scanpy, 10x Visium SpaceRanger
