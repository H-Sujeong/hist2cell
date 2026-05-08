# breast 슬라이드 2장 × lung 가중치 sanity-check 추론 결과 (framework 기반 prep, v2)

> **⚠️ 결과 해석 주의 (caveat)**
>
> 본 추론은 **breast 슬라이드 (KBSMC, 강북삼성)** 2장에 **healthy human lung 으로 학습된 Hist2Cell 가중치** (`humanlung_cell2location_leave_A50_out.pth`) 를 적용한 것이다. 모델이 출력하는 80개 cell type 은 모두 **폐 전용** (AT1, AT2, Basal, Suprabasal, Ciliated, SMG_*, Schwann, Secretory_*, B/T 세포 서브타입 등) 이라 **breast 조직에 대한 생물학적 해석 금지**. 이 보고서의 수치는 **파이프라인 동작 검증 (sanity-check) 용도** 일 뿐이다.

---

## 1. 입력 / 환경

| 항목 | 값 |
|---|---|
| 슬라이드 1 | `Z 2025000042,1-085-12,.svs` (Aperio SVS, 2.5 GB) |
| 슬라이드 2 | `Z 2025000042,1-152-19,,dup1.svs` (Aperio SVS, 2.8 GB) |
| 두 슬라이드 dims | 221,487 × 88,904 px, 4 pyramid level |
| mpp | 0.2615 μm/px (≈ 40× objective) |
| 가중치 | `model_weights/humanlung_cell2location_leave_A50_out.pth` (lung) |
| GPU | 4 × CUDA (분산 inference) |
| openslide | `openslide-bin` 4.0.0.13 + `openslide-python` 1.4.3 (pip wheel, 시스템 lib 불필요) |

---

## 2. 두 단계 prep 비교 (왜 v2 인가)

### v1: hex grid + HSV-saturation Otsu (`prep/prepare_wsi_for_inference.py`)
- 단순 Otsu 가 **슬라이드 라벨 / 잉크 자국** 을 조직으로 잘못 판정 → slide2 에서 실제 조직 대부분 누락 (`spot_view.jpg` 가 라벨 위주만 잡음).
- slide1 결과: 16811 spots / 112533 edges / .pt 10 GB
- slide2 결과: 7545 spots → **상당 부분 라벨 영역, 실제 조직 누락**
- 결정: 사용 안 함.

### v2: framework `ForegroundMasker` + 단순 stride 샘플링 (`prep/prepare_wsi_for_inference_v2.py`)
- **WSI_tile_sampling_framework** 의 `ForegroundMasker` 사용 — YUV (red/green/blue) 다채널 + paraffin/배경 제거 → 라벨/얼룩 회피 정확도 향상.
- framework 의 `TileSampler` 의 contour-based 후처리 (`filter_tiles_by_boundary`) 는 **사용 안 함** — sparse 한 oasis 영역의 작은 contour 를 떨어뜨려 우측 조직을 놓치는 문제 있음.
- 대신 reference (`/home/sjhong/wsitasks/...`) 스타일의 단순 stride 루프 + `min_tissue_frac` 임계 사용.
- thumbnail 은 `openslide.get_thumbnail()` 로 직접 추출 (matplotlib `savefig` 경로의 자동 리스케일 우회).
- 그래프는 hex 그리드 대신 **kNN (k=6) + self-loop + symmetric union** — Visium 가 학습 시 가졌던 6-이웃 위상을 일반 grid 에서도 유지.

### v2 프로세스 단계
1. WSI open (`openslide.OpenSlide`) → level-0 dims 확보
2. thumbnail (long side ≤ `--thumb-max-side` (4000)) → ForegroundMasker 로 0/1 mask
3. `--tile-size` 400 px 단위 stride 로 level-0 좌표 후보 생성, mask 영역 평균 ≥ `--min-tissue-frac` (0.10) 만 keep
4. h5 (`<slide>_coords.h5`) + spots.csv + 정확 스케일 overlay (`spot_view.jpg`) 저장
5. 각 좌표의 중심에서 224×224 패치 추출 + ImageNet 정규화 → `[N, 3, 224, 224]`
6. 좌표로 kNN edge_index 구성 → `Data(x, edge_index, pos, spot_id)` 저장 (`<slide>.pt`)

### framework 수정 (배포본 버그 fix)
- `WSI_tile_sampling_framework/TileSampling.py:159-165` — `TileSampler(tile_size=224, overlap=0, min_tiles=5, ...)` 가 하드코딩되어 사용자 인자 무시 → `self.tile_size / self.overlap / self.min_tiles` 로 변경
- 동일 파일 `__init__` 에 누락된 `self.max_depth = max_depth` 추가 (`run()` 에서 참조됨)
- `ForegroundMasking.py:23-25` — `templates/tcga_brca_template.json` 상대경로 → 모듈 파일 기준 절대경로로 변경 (CWD 와 무관하게 import 가능)

---

## 3. v2 prep 산출물 요약

| 슬라이드 | spots | edges | edges/node | tile_x range | tile_y range | tissue frac (thumb) | .pt 크기 |
|---|---:|---:|---:|---|---|---:|---:|
| slide1_085_12 | **35,821** | 286,755 | 8.0 (incl. self) | 0 – ~210k | 0 – ~85k | — | 21.5 GB |
| slide2_152_19 | **40,502** | 324,604 | 8.0 (incl. self) | 0 – ~205k | 0 – ~80k | 0.298 | 24.4 GB |

QC 이미지:
- slide1: `inference/slide1_085_12_v2/{spot_view.jpg, tissue_mask.png}` — 중앙 분홍 사각 조직 완전 커버, 좌우 inkstain 1 strip 가 false positive 로 일부 포함 (~10%)
- slide2: `inference/slide2_152_19_v2/{spot_view.jpg, tissue_mask.png}` — 중앙 + 우측 조직 모두 커버, 좌측 라벨 sticker false positive (~5–10%)

> v1 (hex+Otsu) 의 slide2 결과와 비교했을 때 **실제 조직 영역 커버리지 대폭 개선됨**. label/inkstain false positive 는 잔여 — 후처리(좌측/우측 margin 컷) 로 정리 가능하나 본 sanity-check 에서는 그대로 진행.

---

## 4. multi-GPU inference (`inference/infer.py`)

- 모델 정의: training tutorial cell 5 의 `Hist2Cell` 클래스 그대로 (cell_dim=80, vit_depth=3). 184 keys, state_dict 완전 매칭.
- 분산 방식: `torch.multiprocessing.spawn` 으로 4 worker, 각자 `cuda:r` 에 binding, `input_nodes = arange(r, n, world_size)` 로 노드 1/4씩 분담.
- 각 GPU 가 자기 shard 만 NeighborLoader 로 (k=2 hop, `[-1,-1]` neighbors, batch_size=16) 처리, fused-head 포함 4 head 평균이 모델 출력.
- 주의 픽스 (작업 중 발견 / 해결):
  - PyG 2.7 의 `NeighborLoader.input_id` 는 `input_nodes` 내 **로컬 인덱스** 를 반환 → 글로벌 인덱스 복구 시 `shard[input_id]` 로 매핑 (코드에 반영)
  - PyG `Data.spot_id` 가 Python 리스트(len==num_nodes)이면 collate_fn 이 거부 → worker 에서 `del data.spot_id` 후 NeighborLoader 에 전달, 메인 프로세스에서 csv 작성 시 다시 사용
  - kNN edge_index 가 numpy `.T` view → non-contiguous → pyg_lib 거부 → `np.ascontiguousarray(edges.T)` 강제

### 실측 시간 (4 × CUDA, 노드 1/4 shard 기준)

| 슬라이드 | 노드 (총) | 노드 (shard) | 전체 wall time | 처리율 (peak/GPU) |
|---|---:|---:|---:|---:|
| slide1_085_12 | 35,821 | 8,956 | **165.6 s** | 53.5 spots/s |
| slide2_152_19 | 40,502 | 10,126 | **251.7 s** | 44.0 spots/s |

---

## 5. 추론 결과 통계

| 슬라이드 | spots | global mean | global max | %zero (ReLU clip) | top-5 cell types (mean) |
|---|---:|---:|---:|---:|---|
| slide1_085_12 | 35,821 | 0.127 | 25.08 | 15.0% | Muscle_smooth_syst_arterial (0.96), AT2 (0.85), Fibro_adventitial (0.71), Fibro_alveolar (0.64), AT1 (0.60) |
| slide2_152_19 | 40,502 | 0.144 | 20.74 | 14.2% | Ciliated (1.22), AT2 (1.10), Fibro_alveolar (0.81), AT1 (0.68), Endothelia_vascular_Cap_a (0.58) |

> top-5 가 모두 lung 세포 — 조직 mismatch 의 직접적 증거. 동일한 lung 가중치를 lung 슬라이드에 적용하면 의미 있는 분포가 나오겠지만, breast 에 적용하면 모델은 결국 학습 분포 안에서 가장 그럴듯한 lung 타입을 출력할 뿐이다.

각 슬라이드 디렉토리에 top-6 cell type spatial heatmap 저장:
- `inference/slide1_085_12_v2/top6_cell_heatmaps.png`
- `inference/slide2_152_19_v2/top6_cell_heatmaps.png`

---

## 6. 산출물 구조

```
inference/
├── slide1_085_12_v2/
│   ├── slide1_085_12.pt              21.5 GB — PyG Data(x, edge_index, pos, spot_id)
│   ├── slide1_085_12_coords.h5       28 KB  — tile coords + metadata (reference style)
│   ├── spots.csv                     1.9 MB
│   ├── tissue_mask.png               48 KB
│   ├── spot_view.jpg                 1.5 MB — 정확 스케일 overlay
│   ├── predictions.csv               31 MB  — 35821 × (spot_id, X, Y, 80 cell types)
│   ├── predictions.npy               11 MB  — float32 [35821, 80]
│   └── top6_cell_heatmaps.png
│
├── slide2_152_19_v2/
│   ├── slide2_152_19.pt              24.4 GB
│   ├── slide2_152_19_coords.h5
│   ├── spots.csv                     2.1 MB
│   ├── tissue_mask.png               89 KB
│   ├── spot_view.jpg                 1.8 MB
│   ├── predictions.csv               35 MB  — 40502 × (spot_id, X, Y, 80 cell types)
│   ├── predictions.npy               13 MB
│   └── top6_cell_heatmaps.png
│
├── slide1_085_12/  (v1 — hex+Otsu, 참조용 잔존)
└── slide2_152_19/  (v1 — 라벨 false positive 가 컸던 케이스)
```

`*.pt` 는 .gitignore 로 차단됨 (개당 20+ GB). predictions.csv/npy 는 작아서 트래킹 가능하나 큰 .pt 와 같은 폴더라 inference/ 전체가 ignore 됨 — 결과 공유는 별도 path 로.

---

## 7. 알려진 한계 및 권장 후속

1. **조직 mismatch (가장 큰 한계)** — breast 데이터에 lung 가중치 적용. 모델 출력은 형식적이며 생물학적 해석 불가. 후속 시 breast 학습 가중치 또는 transfer learning 필요.
2. **Label / inkstain false positive (~5–10% spots)** — 가장 단순한 후처리는 `predictions.csv` 에 `X` 좌표 필터 (slide1: `X > 8000` 등), 또는 prep 단계에서 `--exclude-x-min` 옵션 추가. 본 보고서에서는 정리 안 함.
3. **mpp mismatch (0.26 vs 학습 분포 ~0.5)** — 학습은 ~20× Visium, 입력은 40× Aperio. 모델 입력 224 px 가 절반 크기 (~58 μm) 시야를 보게 됨. 결과 신뢰도에 영향. 학습 분포에 맞추려면 level-1 (4× downsample) 에서 패치 추출하거나 224 입력을 더 큰 crop 후 resize 해야 함.
4. **kNN graph vs hex graph** — k=6 + symmetric + self-loop 로 학습 분포의 ~7 edges/node 에 가깝게 맞췄으나 정확히 같지는 않음. 모델 GAT layer 가 그래프 위상에 어느 정도 견고할 것으로 가정.

---

## 8. 재현 명령

```bash
# 1) prep
python prep/prepare_wsi_for_inference_v2.py \
    --input  "/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-085-12,.svs" \
    --output ./inference/slide1_085_12_v2 \
    --slide-name slide1_085_12 \
    --tile-size 400 \
    --min-tissue-frac 0.10

# 2) inference (4-GPU)
python inference/infer.py \
    --data    inference/slide1_085_12_v2/slide1_085_12.pt \
    --weights model_weights/humanlung_cell2location_leave_A50_out.pth \
    --output  inference/slide1_085_12_v2

# slide2 동일 패턴
```

---

## 9. 관련 파일 정리

| 항목 | 경로 |
|---|---|
| v2 prep 코드 | `prep/prepare_wsi_for_inference_v2.py` |
| v1 prep (deprecated) | `prep/prepare_wsi_for_inference.py` |
| 4-GPU inference | `inference/infer.py` |
| framework (수정본) | `WSI_tile_sampling_framework/TileSampling.py`, `ForegroundMasking.py` |
| 본 보고서 | `report/03_breast슬라이드2장_lung가중치_추론결과_v2framework.md` |
| 선행 보고서 | `report/01_*` (학습 prep 의 외부 파이프라인 정리), `report/02_*` (v1 prep 사용법) |
