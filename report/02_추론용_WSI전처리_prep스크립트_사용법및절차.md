# 추론(inference) 전용 WSI 전처리 — `prep/prepare_wsi_for_inference.py` 사용법 및 절차

## 0. 배경: training prep 와의 차이

학습용 데이터(`data_preparation_tutorial.ipynb`) 는 `Data(x, edge_index, y, pos, spot_id)` 5종을 모두 만든다. 이 중 `y` 는 **250개 gene 발현 + 80개 cell type abundance** 라벨로, 원래 cell2location 과 SpaceRanger 출력에서 만들어진다.

추론 시나리오에서는:

- 우리는 **저장된 weights** (`model_weights/*.pth`) 로 cell type abundance 를 **예측** 하려는 것이므로 `y` 자체가 출력 대상 → **입력에 필요 없음**.
- `model.forward(x, edge_index)` 시그니처에서 보듯, 모델이 실제로 사용하는 입력은 **`x` 와 `edge_index` 뿐**이다 (training tutorial Cell 5).
- `pos` 와 `spot_id` 는 모델 forward 에는 안 쓰이지만 시각화/후처리에 필요해서 같이 저장한다.

따라서 추론용 전처리는 다음 4가지만 만들면 된다.

| 필요 산출물 | 형태 | 용도 |
|---|---|---|
| `x` | `[N, 3, 224, 224]` float32, ImageNet 정규화 | ResNet18 입력 |
| `edge_index` | `[2, E]` long | GAT 의 인접 정보 |
| `pos` | `[N, 2]` float32 | 시각화 좌표 |
| `spot_id` | `list[str]` 길이 N | 결과 추적 |

ST 시퀀싱이 없는 일반 H&E WSI 한 장만 가지고도 만들 수 있다는 점이 핵심.

---

## 1. 스크립트 한 줄 사용법

```bash
python prep/prepare_wsi_for_inference.py \
    --input  /path/to/MY_SLIDE.svs \
    --output ./prep_out/MY_SLIDE
```

### CLI 옵션

| 옵션 | 기본값 | 설명 |
|---|---|---|
| `--input` (필수) | — | WSI 경로. `.svs/.ndpi/.tif/.tiff/.mrxs/.jpg/.png` 등 |
| `--output` (필수) | — | 출력 디렉터리 (없으면 생성) |
| `--slide-name` | 입력 파일 stem | spot ID prefix 로 사용. `<slide>_r<row>c<col>` 형태가 됨 |
| `--patch-size` | `224` | crop 한 변 길이. **모델 가중치가 224 기준이므로 변경 비권장** |
| `--spot-distance` | `200.0` | hex grid 의 spot 중심 간격(px). 자기 WSI 의 mpp 에 맞춰 조정 (§3 참고) |
| `--min-tissue-frac` | `0.30` | 패치 영역의 tissue 비율 하한. 빈 배경 spot 제거 |
| `--thumb-max-side` | `4000` | tissue mask 계산용 썸네일 한 변 최대 |
| `--save-patches` | off | 켜면 `<output>/patches/<spot_id>.jpg` 도 저장 (QC/시각화용) |

---

## 2. 스크립트 내부 단계 (6단계)

### [1/6] WSI 열기
- openslide 가 설치되어 있으면 `openslide.OpenSlide` 로 lazy 하게 read_region.
- 없거나 jpg/png 확장자면 **PIL** 로 전체 이미지 로드 (`Image.MAX_IMAGE_PIXELS = None`).
- 큰 SVS/NDPI 가 메인 사용처라면 `pip install openslide-python` 권장 — 메모리 사용량이 훨씬 작다.

### [2/6] tissue mask
1. WSI 를 `--thumb-max-side` 로 다운샘플.
2. RGB → **HSV saturation** 채널 (배경/유리는 채도가 낮고 조직 H&E 는 높음).
3. **Otsu** 임계값 적용 + 너무 밝은 픽셀(`max(R,G,B) >= 230`) 제외.
4. `scipy.ndimage` 로 closing → fill_holes → opening 후처리.
- 결과: `tissue_mask.png` 로도 저장되어 눈으로 확인 가능.

### [3/6] hex spot grid
Visium 와 동일한 격자 규약을 사용:
- 짝수 행 → `array_col ∈ {0, 2, 4, ...}`
- 홀수 행 → `array_col ∈ {1, 3, 5, ...}` (반칸 offset)
- 행 간격 `dy = spot_distance × √3 / 2` ≈ `0.866 × spot_distance`
- 행 안 간격 `dx = spot_distance`

이 규약 덕에 학습 데이터와 **똑같은 |Δcol|<3, |Δrow|<2 이웃 규칙**을 그대로 쓸 수 있고, 내부 spot 은 자동으로 6 이웃이 된다.

### [4/6] tissue 필터
각 후보 spot 에 대해 224×224 crop 영역을 썸네일 좌표로 매핑하고, 그 영역 안에서 mask 의 평균값(=tissue 비율)이 `--min-tissue-frac` 이상이면 keep.

### [5/6] 패치 추출 + 정규화
- openslide 핸들이면 `read_region((x0,y0), level=0, (224,224))`,
- PIL 핸들이면 `img.crop(...)`.
- 모든 패치에 학습과 동일한 ImageNet 정규화 적용:
  - `mean=[0.485, 0.456, 0.406]`, `std=[0.229, 0.224, 0.225]`
- `[N, 3, 224, 224] float32` 텐서로 누적.

### [6/6] 그래프 + 저장
- `array_col, array_row` 로부터 N×N 인접행렬 생성 (`|Δcol|<3 & |Δrow|<2`, self-loop 포함) → `dense_to_sparse` → `edge_index`.
- `Data(x, edge_index, pos, spot_id)` → `<slide>.pt` 로 `torch.save`.
- 부가 산출물: `spots.csv`, `tissue_mask.png`, `spot_view.jpg` (썸네일 위에 spot circle overlay).

---

## 3. `--spot-distance` 어떻게 정할까

학습 데이터(Visium) 는 spot 중심 간격이 **약 100 μm**, 즉 H&E 20× WSI(0.5 μm/px)에서 **~200 px**. 이 값이 기본값.

다른 배율이라면:

| WSI 배율 | mpp (μm/px) | 권장 `--spot-distance` (≈100μm) |
|---|---|---|
| 40× | 0.25 | `400` |
| 20× | 0.50 | `200` (기본) |
| 10× | 1.00 | `100` |

너무 촘촘하면 패치가 거의 같은 영역이 되고 그래프만 무거워진다. 너무 성기면 GAT 가 받는 spatial context 가 학습 분포와 어긋난다. **학습 시 패치가 약간 겹치는 정도(spot_distance ≲ patch_size)** 가 가장 분포에 가깝다.

---

## 4. 출력 디렉터리 구조

```
<output>/
├── <slide>.pt          # PyG Data(x, edge_index, pos, spot_id) — 모델 입력
├── spots.csv           # spot_id, X, Y, array_col, array_row
├── tissue_mask.png     # QC: 검출된 tissue 영역 이진화
├── spot_view.jpg       # QC: 썸네일 위에 spot 원 overlay
└── patches/            # (--save-patches 시) 개별 jpg
    └── <spot_id>.jpg
```

---

## 5. 만들어진 `.pt` 로 추론 돌리기

`Data` 객체의 필드 이름이 학습 prep 결과와 동일하므로 기존 추론/시각화 노트북에 그대로 끼워 넣을 수 있다. 최소 예시:

```python
import torch
from torch_geometric.data import Data
from torch_geometric.loader import NeighborLoader

# 1) 모델 정의는 tutorial_training/training_tutorial.ipynb 의 Hist2Cell 클래스 그대로
from your_module import Hist2Cell

device = "cuda" if torch.cuda.is_available() else "cpu"
model = Hist2Cell(cell_dim=80).to(device).eval()
model.load_state_dict(torch.load("model_weights/humanlung_cell2location_leave_A50_out.pth",
                                 map_location=device))

# 2) prep 출력 로드
data = torch.load("prep_out/MY_SLIDE/MY_SLIDE.pt", weights_only=False)

# 3) NeighborLoader 로 subgraph batching (GPU 메모리에 맞게)
loader = NeighborLoader(data, num_neighbors=[-1, -1], batch_size=16,
                        directed=False, shuffle=False)

# 4) 추론
preds = torch.zeros((data.num_nodes, 80))
with torch.no_grad():
    for sub in loader:
        sub = sub.to(device)
        out = model(sub.x, sub.edge_index)          # 모델 출력은 학습 시 정의에 따름
        # 학습 코드에서 반환한 cell_pred 가 (B, cell_dim) 인 spot 단위 출력이라고 할 때:
        center_idx = sub.input_id                    # 이번 배치의 center node 글로벌 인덱스
        preds[center_idx] = out[: sub.batch_size].cpu()

# 5) 좌표와 함께 저장
import pandas as pd
df = pd.DataFrame(preds.numpy())
df.insert(0, "spot_id", data.spot_id)
df["X"] = data.pos[:, 0].numpy()
df["Y"] = data.pos[:, 1].numpy()
df.to_csv("prep_out/MY_SLIDE/predictions.csv", index=False)
```

> 학습 노트북에서 `Hist2Cell.forward` 가 **여러 head** (spot/local/fused) 를 반환하므로 실제 추론 시 어떤 head 를 쓸지는 학습 코드 정의를 그대로 따르면 된다. tutorial 의 fused head 가 최종 보고용 출력이다.

---

## 6. 한계와 주의사항

- **세포 종류는 학습된 80종 그대로**. 학습에 쓴 cell type 어휘(humanlung_cell2location 의 `cell_types.pkl`) 와 동일한 출력 차원만 의미가 있다. 다른 조직(예: skin, breast) 가중치를 쓰면 80종이 그 조직 기준으로 해석된다.
- **gene 출력 250종은 학습 데이터의 슬라이드들에서 정해진 상위 발현 gene 집합**. 다른 데이터셋에 그대로 쓰면 의미 해석이 달라진다.
- **mpp 일치가 중요**. 학습 분포와 한참 다른 배율이라면 결과 신뢰도가 떨어진다. 가능하면 WSI 를 학습 분포에 맞춰 리샘플 후 입력.
- **tissue mask 가 너무 빈약하면** `--min-tissue-frac` 을 낮추거나 `--thumb-max-side` 를 키워 더 정밀한 mask 를 얻는다. 반대로 배경에 spot 이 들어오면 `--min-tissue-frac` 을 0.5–0.7 까지 올린다.
- **PIL 로 거대 WSI 를 열면 RAM 폭주** — 30k×30k RGB ≈ 2.7GB. SVS/NDPI 사용자는 `pip install openslide-python` (시스템에 `libopenslide` 도 필요) 를 강력 권장.

---

## 7. 검증

본 repo 의 예시 슬라이드 `WSA_LngSP9258467_low_res.jpg` (8208×4608) 로 smoke test:

```
[6/6] building graph and writing outputs
Done. Wrote:
  /tmp/h2c_prep_test/WSA_LngSP9258467_low_res.pt   (169 nodes, 1051 edges)
```

저장된 Data:
```
Data(x=[169, 3, 224, 224], edge_index=[2, 1051], pos=[169, 2], spot_id=[169])
edges/node = 6.22  → 헥사그리드 위상 정상
x dtype = float32, ImageNet 분포 (mean ≈ 1.26, std ≈ 0.76 of normalized values)
```
스크립트가 모델 입력 사양과 학습 prep 노트북의 출력 형태를 그대로 따르고 있음을 확인.

---

## 8. 파일 위치 정리

- 코드: `prep/prepare_wsi_for_inference.py`
- 본 보고서: `report/02_추론용_WSI전처리_prep스크립트_사용법및절차.md`
- 관련 선행 문서: `report/01_원본_WSI전처리파이프라인_튜토리얼이전단계_요약.md` (학습 prep 가 가정하는 외부 파이프라인 정리)
