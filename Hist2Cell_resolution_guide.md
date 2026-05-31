# Hist2Cell를 다른 Resolution / Magnification / 플랫폼에 적용하기

> Hist2Cell 공식 구현( `Weiqin-Zhao/Hist2Cell` )의 실제 코드를 기준으로,
> 입력 해상도·배율·spot geometry가 다른 데이터를 정확히 처리하기 위한 실무 가이드.
> 핵심 결론부터: **Hist2Cell 모델 자체에는 resolution/magnification 파라미터가 없다.**
> 배율 대응은 전적으로 **데이터 전처리 단계(patch crop + graph 구성)**에서 사용자가 맞춰야 한다.

---

## 0. 왜 "how-to"가 없어 보이는가

Hist2Cell이 학습/추론에서 실제로 보는 입력은 **물리적 해상도(µm)가 아니라 고정된 224×224 픽셀 텐서**다.
- H&E 이미지는 20× magnification으로 촬영되었고, Visium ST 자체는 10× 으로 수행되었다.
- 모델 입력 노드는 각 ST spot을 중심으로 crop한 224×224 RGB 패치이며, 모델은 그 패치 안에 몇 µm가 담겼는지 **인식하지도, 보정하지도 않는다.**

따라서 다른 배율/MPP/플랫폼에 쓰려면 모델이 아니라 **전처리에서 두 가지를 맞춰야 한다.**
1. **패치의 물리적 시야(FOV, µm 기준)** — magnification / MPP 대응
2. **그래프 edge 정의** — spot geometry(간격·격자 구조) 대응

논문·코드에 "다른 배율은 이렇게 하라"가 일반화된 절차로 안 적혀 있을 뿐, 코드를 뜯어보면 어디를 수정해야 하는지는 명확하다. 아래가 그 지점이다.

---

## 1. 데이터 객체 구조 (수정 대상의 전체 그림)

전처리 결과물은 PyTorch Geometric `Data` 객체이며, 슬라이드 1장 = 그래프 1개다.

| 필드 | shape (예시) | 의미 | 배율/해상도 영향 |
|------|------|------|------|
| `x` | `[N, 3, 224, 224]` | spot별 H&E 패치 (ResNet18 입력) | **★ 직접 영향** — crop FOV가 여기에 들어감 |
| `edge_index` | `[2, E]` | spot 간 spatial neighbor (COO) | **★ 직접 영향** — spot 간격/격자로 정의 |
| `y` | `[N, 330]` | label = 250 genes + 80 cell types | 영향 없음(라벨) |
| `pos` | `[N, 2]` | spot의 **픽셀 좌표** (시각화·분석용) | 스케일만 일치시키면 됨 |

> `N` = spot 수, `E` = edge 수. 모델 가중치는 224×224·채널3·label차원에만 묶여 있고, **µm나 배율에는 묶여 있지 않다.** → 입력 정규화만 맞추면 재학습 없이 추론 가능, 도메인 차이가 크면 fine-tuning.

---

## 2. 표준 파이프라인에서 배율과 직결되는 두 코드 지점

### 2-1. 패치 생성 (`x`) — magnification / MPP 대응 지점

공식 `STDataset`은 디스크에 미리 잘려있는 `patches/*.jpg`를 읽어 **무조건 224로 resize**한다:

```python
# data_preparation_tutorial.ipynb - STDataset.__getitem__
patch = Image.open(patch_path).convert('RGB')
data  = transforms.Resize((224, 224))(patch)   # ← 물리적 크기 보정이 전혀 없음
```

즉 **"패치 1장이 실제 몇 µm를 담는가"는 patch를 만드는 단계(WSI에서 crop할 때)에서 이미 결정**되며,
공식 코드는 그 crop을 외부 도구(**DSMIL-WSI 파이프라인**)에 맡긴다. resize는 단지 224 픽셀로 맞추는 것뿐이다.

> **여기가 핵심.** 다른 magnification/MPP를 쓸 때 바꿔야 할 것은 `Resize((224,224))`가 아니라 **crop 단계의 crop 크기(픽셀)**다. 학습 도메인과 동일한 **물리적 FOV(µm)**가 224 픽셀 안에 담기도록 crop 픽셀 수를 환산해야 한다.

### 2-2. 그래프 edge 생성 (`edge_index`) — spot geometry 대응 지점

공식 코드는 픽셀 거리가 아니라 **Visium의 array 좌표(`array_col`, `array_row`) 격자 인접성**으로 edge를 만든다:

```python
# data_preparation_tutorial.ipynb - 인접행렬 구성
for i in range(num_spots):
    for j in range(num_spots):
        if i == j:
            adj[i][j] = 1.0                      # self-loop
        else:
            x1, y1 = spot_array_x_y[i]           # (array_col, array_row)
            x2, y2 = spot_array_x_y[j]
            col_dist = abs(x2 - x1)
            row_dist = abs(y2 - y1)
            if col_dist < 3 and row_dist < 2:    # ← Visium 육각격자 6-이웃 규칙
                adj[i][j] = 1.0
# adj → dense_to_sparse → edge_index
```

`col_dist < 3 and row_dist < 2`는 **Visium 육각 격자에서 6개 최근접 이웃**을 잡는 규칙이다(Visium은 한 칸 건너 열이 한 줄이라 col이 ±2까지 이웃).

> **여기가 두 번째 핵심.** 다른 플랫폼(다른 spot 간격, 사각 격자, 비격자 좌표, super-resolution 가상 spot)에서는 이 규칙이 깨진다. array 좌표가 없거나 격자가 다르면 **픽셀 좌표 기반 거리/KNN으로 edge 정의를 바꿔야** 한다.

### 2-3. 픽셀 좌표 (`pos`) — 스케일 일치만 주의

`pos`는 시각화·colocalization 분석용 픽셀 좌표다. super-resolution 튜토리얼에서 보듯, 좌표가 저장된 해상도와 표시용 이미지 해상도가 다르면 스케일 보정이 필요하다(예: 좌표가 4× 이미지 기준이면 저해상도 위에 그릴 때 ÷4). 모델 예측 자체에는 영향 없지만, 분석 단계에서 어긋나면 결과 해석이 틀어진다.

---

## 3. 케이스별 How-To

### Case A. magnification만 다름 (예: 학습=20×, 신규=40×)

배율이 다르면 같은 224 픽셀이 담는 µm가 달라진다 → **동일 FOV가 되도록 crop 크기 환산 후 224로 resize**.

1. 학습 도메인의 FOV 확인: 20× H&E의 MPP가 약 0.5 µm/px라면 224 px ≈ **112 µm**.
2. 신규 슬라이드의 MPP로 같은 112 µm에 해당하는 픽셀 수 계산.
   `crop_px = target_FOV_µm / MPP_new`
   예: 40×(≈0.25 µm/px) → `112 / 0.25 = 448 px`로 crop → 224로 downsize.
3. crop은 **각 spot의 픽셀 좌표를 중심**으로 수행(아래 §4 코드).

> 직관: 모델은 "조직이 화면을 얼마나 채우는지"를 학습했다. 배율이 올라가면 더 넓은 픽셀 영역을 잘라 같은 시야를 만들어줘야 학습 분포와 맞는다. (GHIST가 112로 center-crop 후 256으로 resize해 단일세포 모델의 해상도를 맞춘 것과 같은 발상.)

### Case B. MPP만 다름 (배율 라벨은 같지만 스캐너가 다름)

magnification 라벨(20×, 40×)은 스캐너마다 실제 MPP가 다르다. **항상 라벨이 아니라 MPP(µm/px)를 신뢰**하라. WSI 메타데이터에서 MPP를 읽어 Case A의 공식으로 crop 픽셀을 정한다.

```python
import openslide
slide = openslide.OpenSlide("sample.svs")
mpp_x = float(slide.properties.get('openslide.mpp-x'))   # µm/px
```

### Case C. spot 간격/플랫폼이 다름 (Visium HD, Xenium, Slide-seq, 비격자 등)

`array_col/array_row` 격자 규칙(§2-2)이 안 맞으므로 **edge를 픽셀 좌표 기반으로 재정의**한다. 두 가지 방법:

- **반경(radius) 기반**: spot 간 중심거리 `d`의 약 1.5배를 임계로.
- **KNN 기반**: 학습 도메인 평균 차수(약 6, self-loop 제외)에 맞춰 k=6.

```python
import numpy as np
from sklearn.neighbors import NearestNeighbors

def build_edges_from_pixels(pos_xy, k=6, add_self_loop=True):
    """pos_xy: [N,2] 픽셀(또는 µm) 좌표. Visium 격자 규칙을 대체."""
    nn = NearestNeighbors(n_neighbors=k+1).fit(pos_xy)   # +1: 자기 자신 포함
    _, idx = nn.kneighbors(pos_xy)
    src, dst = [], []
    for i in range(len(pos_xy)):
        for j in idx[i]:
            if i == j and not add_self_loop:
                continue
            src.append(i); dst.append(j)
    edge_index = np.vstack([src, dst])
    # 무방향 대칭화
    edge_index = np.hstack([edge_index, edge_index[::-1]])
    edge_index = np.unique(edge_index, axis=1)
    return edge_index   # [2, E]
```

> 목표는 **학습 그래프의 평균 차수(~6)와 비슷한 국소 연결성**을 재현하는 것. 차수가 크게 달라지면 GAT/Transformer가 받는 이웃 통계가 바뀌어 성능이 흔들린다.

### Case D. super-resolution (더 촘촘한 가상 spot)

공식 2× 튜토리얼은 학습된 모델을 **2× 해상도 패치 그래프**에 그대로 추론한다. 핵심은:
- 가상 spot 좌표(원 spot 사이 보간 위치)를 만들고, 각 좌표 중심으로 **동일 FOV** 패치를 crop(§3 Case A 논리 동일).
- edge는 새 좌표 집합 위에서 다시 구성(§2-2 또는 §3-C).
- 추론은 2-hop subgraph 단위라 메모리 부담이 작아 4×/8×도 시간만 더 들 뿐 가능.

---

## 4. 좌표 기반 패치 crop 레퍼런스 코드 (배율 환산 포함)

공식 파이프라인의 "미리 잘린 patches/" 의존을 없애고, **WSI에서 직접 동일 FOV로 crop**하는 범용 함수.

```python
import openslide
import numpy as np
from PIL import Image

def crop_patches_fixed_fov(
    wsi_path,
    spot_pixel_coords,        # [N,2] level-0 픽셀 좌표 (각 spot 중심)
    target_fov_um=112.0,      # 학습 도메인 FOV. 20×·0.5µm/px·224px ≈ 112µm
    out_px=224,               # 모델 입력 크기(고정)
):
    slide = openslide.OpenSlide(wsi_path)
    mpp = float(slide.properties.get('openslide.mpp-x'))   # µm/px (level 0)
    crop_px = int(round(target_fov_um / mpp))              # 동일 FOV에 해당하는 crop 픽셀
    half = crop_px // 2

    patches = []
    for (cx, cy) in spot_pixel_coords:
        x0, y0 = int(cx - half), int(cy - half)
        region = slide.read_region((x0, y0), 0, (crop_px, crop_px)).convert('RGB')
        region = region.resize((out_px, out_px), Image.BILINEAR)   # 224로 통일
        patches.append(np.asarray(region))
    return np.stack(patches)   # [N, 224, 224, 3]
```

이후 표준 정규화(ImageNet 통계)는 공식 코드와 동일하게 적용:

```python
from torchvision import transforms
norm = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485,0.456,0.406], std=[0.229,0.224,0.225]),
])
```

---

## 5. 전체 흐름 (다른 해상도 데이터 → Hist2Cell 입력)

```
WSI(.svs/.tif) + spot 좌표(픽셀) + (선택)gene/cell 라벨
        │
        ├─ [§4] MPP 읽기 → 동일 FOV(µm)로 crop_px 환산 → crop → 224 resize → 정규화 ──► x [N,3,224,224]
        │
        ├─ [§2-2 / §3-C] spot geometry 판단
        │        ├─ Visium 격자면      → array_col/row 인접 규칙
        │        └─ 그 외/비격자/SR면   → 픽셀좌표 KNN(k≈6) 또는 radius
        │                                                              ──► edge_index [2,E]
        │
        ├─ 라벨 있으면 250 genes + 80 cell types 병합 (없으면 추론 전용)  ──► y [N,330]
        │
        └─ spot 픽셀 좌표(시각화 해상도와 스케일 일치)                    ──► pos [N,2]
                                   │
                                   ▼
        torch_geometric.data.Data(x, edge_index, y, pos)  →  .pt 저장
                                   │
                                   ▼
        NeighborLoader(hop=2, num_neighbors=[-1]*hop)  →  Hist2Cell 추론/학습
```

---

## 6. 체크리스트 (배율을 바꿀 때 반드시 확인)

- [ ] **FOV 일치**: 신규 패치의 µm 시야 ≈ 학습 도메인(예 112µm)인가? (magnification 라벨 말고 MPP로 계산)
- [ ] **224 통일**: crop 픽셀 수와 무관하게 최종 입력은 224×224인가?
- [ ] **정규화 동일**: ImageNet mean/std로 정규화했는가?
- [ ] **edge 차수**: 평균 차수가 학습 그래프(~6, self-loop 제외)와 비슷한가?
- [ ] **self-loop**: 학습과 동일하게 self-loop 포함 규칙을 따랐는가?
- [ ] **무방향 대칭**: `edge_index`가 대칭(undirected)인가?
- [ ] **pos 스케일**: 시각화 이미지 해상도와 `pos` 좌표 스케일이 맞는가?
- [ ] **라벨 차원**: 라벨을 쓸 경우 250+80=330 순서/차원이 맞는가?
- [ ] **도메인 차이**: stain/조직/장기 차이가 크면 zero-shot 대신 fine-tuning 고려.

---

## 7. 자주 틀리는 지점

1. **magnification 라벨을 그대로 믿음** → 스캐너마다 실제 MPP가 다르다. 항상 MPP로 환산.
2. **224 resize만 하고 FOV를 안 맞춤** → 40× 패치를 그냥 224로 줄이면 학습보다 좁은 시야가 들어가 분포가 어긋난다. crop 픽셀부터 환산해야 한다.
3. **다른 플랫폼인데 Visium 격자 규칙(`col<3, row<2`)을 그대로 사용** → edge가 거의 안 생기거나 과하게 생긴다. 픽셀 KNN/radius로 교체.
4. **`pos` 스케일 불일치** → 예측은 맞는데 시각화에서 점이 엉뚱한 곳에 찍힌다.
5. **도메인 갭 무시** → 다른 장기/염색 프로토콜은 zero-shot 성능이 떨어질 수 있다. 소량 fine-tuning이 안전.

---

### 참고
- 코드 근거: `Weiqin-Zhao/Hist2Cell` — `tutorial_data_preparation/data_preparation_tutorial.ipynb`(STDataset, 인접행렬), `tutorial_analysis_evaluation/super_resovled_cell_abundance_tutorial.ipynb`(2× 추론·좌표 스케일).
- 패치 추출 권장 도구: DSMIL-WSI (`binli123/dsmil-wsi`).
- 모델: ResNet18 인코더 + GATv2Conv + Transformer, 입력 224×224 고정, 2-hop subgraph 추론.
