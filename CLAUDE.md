# CLAUDE.md — Hist2Cell repo 작업 지도

WSI(.svs) → **tile → graph(.pt) → inference → analysis** 파이프라인.
새 세션이 tiling / graph / inference 작업을 바로 찾을 수 있도록 코드·문서 위치를 정리한다.

## 실행 환경
- Python: **`/home/sjhong/hist2cell/.venv/bin/python`** — openslide·torch·PyG·cv2·h5py 등 설치됨
- conda `base310` 등 다른 env 는 의존성 누락 → 쓰지 말 것
- GPU: 4× CUDA. inference 는 GPU 필요 (`inference/infer.py` 가 CUDA 가정).
  `torch.cuda.is_available()` 가 False 면 `nvidia-smi` 확인 (NVML 죽었으면 머신 restart)

## 1. Tiling — WSI → 224/112 타일 좌표 (HDF5)
- 프레임워크: `WSI_tile_sampling_framework/` — tissue mask + tile 좌표
- 드라이버: `WSI_tile_sampling_framework/run_tiling_tcga_luad.py`
  (인자 `--tile-size`, `--output`. CLI `tile_processing.py` 는 positional-arg 버그 있어 **사용 금지**)
- **how-to 문서**: `report/05_WSI타일링_프레임워크_224타일_사용법및절차.md`
- `Tiling.py` 의 `filter_tiles_by_boundary` centroid 버그는 **수정 완료**
  (`round(tile_size)` → `round(tile_size/downscale_factor)`)

## 2. Graph — 타일/WSI → PyG `Data(.pt)`
- `prep/build_graph_from_tiles.py` — 기존 tiling h5 → `.pt` (h5 의 tile_size 읽어 자동 224 resize)
- `prep/prepare_wsi_for_inference_v2.py` — WSI 한 장 → tissue mask + tiling + graph 한 번에
- `prep/prepare_wsi_for_inference.py` — v1 (deprecated, hex grid)
- **문서**: `report/02`(prep 사용법), `report/04`(WSI→분석 end-to-end 쿡북)
- `.pt` = `Data(x[N,3,224,224] ImageNet-norm, edge_index, pos, spot_id)`.
  로드 시 `weights_only=False` 필수

## 3. Inference — `.pt` → cell-type 예측
- `inference/infer.py` — multi-GPU Hist2Cell 추론 → `predictions.csv` / `.npy`
  - `--data <.pt> --weights <.pth> --output <dir>` (옵션: `--world-size`, `--batch-size`)
- 가중치 `model_weights/`: `humanlung_cell2location_leave_A50_out.pth` (기본),
  `..._A37_out.pth`, `demo_ckpt.pth`
- cell type 이름(80종): `example_data/humanlung_cell2location/cell_types.pkl`

## 4. Analysis
- 공간 분석 산출물: `inference/analysis_spatial/`, `inference/analysis_filtered/`
- 워크플로: `report/04` §4 참조

## report/ 문서 색인
| 파일 | 내용 |
|---|---|
| `report/01_*` | 원본 WSI 전처리 파이프라인 요약 |
| `report/02_*` | 추론용 prep 스크립트 사용법 |
| `report/03_*` | breast 2장 lung-가중치 추론 결과 사례 |
| `report/04_*` | WSI→분석 end-to-end 쿡북 |
| `report/05_*` | WSI 타일링 프레임워크 224 타일 how-to |

## 진행 중 작업
- **`lung_pilot/`** — TCGA-LUAD 3장 pilot (Hist2Cell vs HEX 모델 비교).
  단계별 상태·다음 단계·산출물 위치는 **`lung_pilot/README.md`** 참조
