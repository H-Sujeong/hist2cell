# graph_output — TCGA-LUAD PyG 그래프 (.pt)

`tilitng_output/` 의 타일 좌표 HDF5 + 원본 SVS 로부터 만든 **PyG `Data` 그래프**.
패치 단위 그래프 모델(Hist2Cell 등)의 입력으로 바로 쓸 수 있다.

빌드 스크립트: `prep/build_graph_from_tiles.py`

## 두 세트 — 224 vs 112

| 폴더 | 입력 h5 | 패치 추출 | 물리 크기 | 용도 |
|---|---|---|---|---|
| `224/` | `tilitng_output/224/TCGA-LUAD/*.h5` | 224 px @ level 0, resize 없음 | 112 µm | **Hist2Cell** (20× 학습) |
| `112/` | `tilitng_output/112/*.h5` | 112 px @ level 0 → **224 로 ×2 resize** | 56 µm | **HEX 모델** (40× 224 학습) |

TCGA-LUAD 3장은 native 20× (mpp ≈ 0.502). 40× 학습 모델은 56 µm 시야를 기대하므로
20× 슬라이드에서 112 px 를 떠서 224 로 키운다. `x` 는 두 세트 모두 `[N,3,224,224]`.

## 파일

```
graph_output/
├── 224/
│   ├── <slide>.pt            PyG Data(x, edge_index, pos, spot_id)
│   └── <slide>_spots.csv     spot_id, X, Y, tile_x/y_topleft, tile_size
└── 112/
    └── (동일 구조)
```

| 슬라이드 | 224 nodes / edges / .pt | 112 nodes / edges / .pt |
|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 / 23,047 / 1.6 GB | 11,470 / 92,142 / 6.4 GB |
| TCGA-05-4245-01A-01-TS1 | 1,871 / 15,201 / 1.0 GB | 7,499 / 60,407 / 4.2 GB |
| TCGA-05-4390-01A-01-BS1 | 10,661 / 85,959 / 6.0 GB | 42,615 / 344,355 / 23.9 GB |

edges/node ≈ 8 (kNN k=6 + self-loop + 대칭 union).

## `.pt` 내용

`torch_geometric.data.Data`:

| 필드 | shape / type | 설명 |
|---|---|---|
| `x` | `[N, 3, 224, 224]` float32 | 패치, **ImageNet 정규화** (mean .485/.456/.406, std .229/.224/.225) |
| `edge_index` | `[2, E]` int64, contiguous | tile 중심 kNN(k=6) 그래프, 양방향 + self-loop |
| `pos` | `[N, 2]` float32 | tile **중심** (x, y), level-0 px |
| `spot_id` | list[str] 길이 N | `<slide>_x<cx>y<cy>` |

```python
import torch
data = torch.load("graph_output/224/TCGA-05-4245-01A-01-BS1.pt", weights_only=False)
# Data(x=[2869,3,224,224], edge_index=[2,23047], pos=[2869,2], spot_id=[2869])
```

> `spot_id` 가 python list 라 `weights_only=False` 필수. NeighborLoader 에 넣기 전
> `del data.spot_id` 가 필요할 수 있다 (`report/04` §3.8 참조).

## 주의

- **ImageNet 정규화**: repo 표준(Hist2Cell ResNet18 · DINOv2 모두 ImageNet 통계).
  정규화는 채널별 affine 이라 `x*std+mean` 으로 [0,1] 패치 복원 가능 — HEX 모델이
  다른 전처리를 기대하면 이 역변환 후 재정규화하면 된다.
- **112 세트 해상도 OOD**: TCGA-LUAD 는 native 20× — 진짜 40× 디테일이 없다.
  112→224 업샘플 패치는 면적(56 µm)은 HEX 학습과 맞지만 해상도는 더 흐리다.
  HEX 출력·UMAP 해석 시 명시할 것.
- **224 ↔ 112 그리드 비정렬**: 두 세트의 tile 중심이 일치하지 않는다
  (224 중심 = coord+112, 112 중심 = coord+56). spot 단위 paired 비교가 필요하면
  공유 center 에서 두 크기 패치를 추출하는 별도 단계가 필요하다.

## 재생성

```bash
python prep/build_graph_from_tiles.py \
    --tiles-dir lung_pilot/tilitng_output/224/TCGA-LUAD \
    --wsi-dir   /mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi \
    --output    lung_pilot/graph_output/224
# 112 는 --tiles-dir .../tilitng_output/112  --output .../graph_output/112
```
