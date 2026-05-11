# WSI 한 장에서 cell-type 공간 분석까지 — Cookbook

> **무엇을 위한 문서인가**
>
> 새 H&E WSI (`.svs / .ndpi / .tif / .jpg / .png` …) 한 장이 들어왔을 때 **prep → inference → analysis** 까지 전부 돌리고, 결과를 동료에게 공유 가능한 산출물 (`predictions.csv`, spatial heatmap PNG, Moran's R clustermap, abundance group table) 까지 만드는 흐름. 각 명령어의 모든 옵션을 정리해서, 사용자가 자기 케이스에 맞춰 조합해서 돌릴 수 있도록 한다.

```
WSI ─ prep_v2 ─► <slide>.pt ─ infer ─► predictions.csv ─ analyze ─► spatial PNG + Moran R + abundance CSV
                 (PyG Data)         (multi-GPU)        (.npy + .h5)
```

---

## 0. 먼저 결정할 4 가지

| 결정 | 영향 | 빠른 결정 가이드 |
|---|---|---|
| (a) 사용할 가중치 | 출력 cell type 분포 | breast/일반 sanity-check → `model_weights/humanlung_cell2location_leave_A50_out.pth`. 폐 데이터면 그대로 의미. 다른 조직이면 결과는 패턴/그룹 단위로만 해석 |
| (b) `--tile-size` (prep) | spot 격자 크기 | Visium 학습 분포 ~100 μm 에 맞추기 → `--tile-size = round(100 / mpp)`. 0.26 mpp 면 ~400 px |
| (c) `--min-tissue-frac` | spot 수와 background false positive trade-off | 깔끔한 슬라이드 0.30, 라벨/inkstain 많은 슬라이드 0.40, 신호 빠지면 0.10 |
| (d) GPU 수 | inference 시간 | 디폴트 = 모든 가용 GPU. 4 GPU 면 35k spot 약 3 분, 단일 GPU 면 12 분 |

---

## 1. Prerequisites (1회 셋업)

가상환경은 **uv** 로 관리된 `.venv` (`/home/sjhong/hist2cell/.venv`). 첫 실행 전 다음을 확인.

```bash
# 가상환경 활성 (또는 직접 path 호출)
source /home/sjhong/hist2cell/.venv/bin/activate     # (대안) 직접 .venv/bin/python3 사용

# 핵심 모듈 체크
python -c "import torch, torch_geometric, openslide, scipy, h5py, seaborn; \
           print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'devices', torch.cuda.device_count())"
# 기대 출력 예시: torch 2.5.1+cu121 cuda True devices 4

# 누락 시 설치 (uv 사용)
VIRTUAL_ENV=/home/sjhong/hist2cell/.venv uv pip install \
    openslide-bin openslide-python opencv-python-headless h5py natsort \
    scikit-image scikit-learn seaborn pymupdf
```

> openslide 는 시스템 lib 없이 `openslide-bin` (pip wheel) 만으로 동작. `apt` 권한 없는 컨테이너에서도 OK.

---

## 2. STEP 1 — Prep (`prep/prepare_wsi_for_inference_v2.py`)

WSI 한 장 → tissue mask + tile coords + 224×224 patch tensor + kNN graph → PyG `Data(.pt)`.

### 2.1 풀 옵션

```text
--input INPUT                  (필수) WSI 경로 (.svs/.ndpi/.tif/.tiff/.mrxs/.jpg/.png)
--output OUTPUT                (필수) 출력 디렉토리
--slide-name SLIDE_NAME        spot_id prefix. 디폴트 = 입력 stem 의 공백·쉼표 제거
--tile-size TILE_SIZE          격자 spacing (level-0 px). 디폴트 400 (≈100 μm at 0.26 mpp)
--patch-size PATCH_SIZE        모델 입력 crop 크기. 디폴트 224 — 변경 시 가중치와 mismatch
--min-tissue-frac FRAC         tile 영역 안 tissue mask 평균 ≥ FRAC 만 keep. 디폴트 0.10
--thumb-max-side N             tissue mask 계산용 thumbnail 한 변 최대. 디폴트 4000
--knn K                        spatial graph 의 이웃 수 (self 제외). 디폴트 6
--save-patches                 (옵션) 개별 jpg 패치도 <output>/patches/ 에 저장 (QC용, 디스크 큼)
```

### 2.2 mpp 별 `--tile-size` 권장

```bash
# 슬라이드의 mpp 확인 (한 번만)
python -c "
import openslide
s = openslide.OpenSlide('YOUR_SLIDE.svs')
print('dims:', s.dimensions, 'levels:', s.level_count)
print('mpp_x:', s.properties.get('openslide.mpp-x'),
      'objective:', s.properties.get('openslide.objective-power'))
"
```

| mpp | objective | `--tile-size` 권장 (≈100 μm 격자) |
|---|---|---|
| 0.25 | 40× | `400` |
| 0.50 | 20× | `200` |
| 1.00 | 10× | `100` |

학습 분포 (Visium 20×, ~150 μm spacing) 에 더 가깝게 가려면 위 표 ×1.5. 학습보다 **공간 패턴 우선** 이면 ×0.5–1.0 도 무방.

### 2.3 `--min-tissue-frac` 가이드

| 슬라이드 상태 | 권장 | 이유 |
|---|---|---|
| 깔끔한 H&E, 라벨/잉크 적음 | `0.30` | tissue 위주만 잡힘, false positive 적음 |
| 라벨/inkstain 큰 슬라이드 | `0.40` | 가장자리 artifact 컷 |
| 조직 자체가 sparse (skin, stroma) | `0.10` | 누락 방지 |
| spot 수가 너무 많아 inference 시간 부담 | `0.30~0.50` 으로 ↑ | 50% spot 감소 가능 |

`spot_view.jpg` 보고 조정한 뒤 다시 돌리는 것이 정석.

### 2.4 일반 케이스

```bash
python prep/prepare_wsi_for_inference_v2.py \
    --input  /mnt/path/to/MY_SLIDE.svs \
    --output ./inference/my_slide_v2 \
    --slide-name my_slide \
    --tile-size 400 \
    --min-tissue-frac 0.40
```

### 2.5 출력 구조

```
inference/my_slide_v2/
├── my_slide.pt                # PyG Data(x, edge_index, pos, spot_id) — STEP 2 입력
├── my_slide_coords.h5         # tile coords + metadata (다른 모달리티 매칭 용)
├── spots.csv                  # spot_id, X, Y, tile_x/y_topleft
├── tissue_mask.png            # mask QC (thumbnail 해상도)
└── spot_view.jpg              # thumbnail 위 tile 박스 overlay (정확 스케일) ← QC 1순위
```

### 2.6 QC 체크리스트

`spot_view.jpg` 를 열어 다음 확인:
- [ ] tissue 가 빨간 사각형으로 거의 다 덮였는가? (덮였다 = pass)
- [ ] 슬라이드 라벨/inkstain 도 빨간 사각형으로 잡혔는가? (5–10% 미만 = 수용, 그 이상이면 `--min-tissue-frac` ↑)
- [ ] tissue 가 누락된 영역이 있는가? (있으면 `--min-tissue-frac` ↓)
- [ ] 콘솔의 `kept N tiles` 가 합리적인가? (대형 슬라이드 ~10k–50k 정상, 100k 이상이면 `--tile-size` ↑ 또는 `--min-tissue-frac` ↑)

---

## 3. STEP 2 — Inference (`inference/infer.py`)

PyG `Data(.pt)` → `predictions.csv` + `predictions.npy`. 4-GPU 분산.

### 3.1 풀 옵션

```text
--data DATA                    (필수) prep 산출 .pt 파일 경로
--weights WEIGHTS              (필수) Hist2Cell 가중치 .pth (model_weights/*.pth)
--output OUTPUT                (필수) 출력 디렉토리 (보통 prep 와 같은 폴더 가리키면 됨)
--cell-types CELL_TYPES        cell type 이름 list pickle. 디폴트 example_data/humanlung_cell2location/cell_types.pkl
--batch-size BATCH_SIZE        NeighborLoader batch (center node 수 / GPU). 디폴트 16 — 24GB GPU 기준
--hop HOP                      NeighborLoader hop. 디폴트 2 (학습과 동일)
--world-size WORLD_SIZE        사용할 GPU 수. 디폴트 = `torch.cuda.device_count()`
```

### 3.2 가중치 선택

| 파일 | 학습 cohort | 권장 |
|---|---|---|
| `humanlung_cell2location_leave_A50_out.pth` | Madissoon human lung, donor A50 leave-out | 일반 sanity-check 1순위 |
| `humanlung_cell2location_leave_A37_out.pth` | 동일, donor A37 leave-out | 비교 / 앙상블 |
| `demo_ckpt.pth` | training tutorial 의 데모 (5 epoch) | full-trained 보다 약함, 가벼운 테스트만 |

### 3.3 GPU memory ↔ batch size

| GPU VRAM | `--batch-size` 권장 |
|---|---:|
| 8 GB | `4` |
| 12 GB | `8` |
| 16 GB | `12` |
| 24 GB+ | `16` (디폴트) |

OOM 나면 `--batch-size` 절반으로 ↓.

### 3.4 일반 케이스

```bash
python inference/infer.py \
    --data    inference/my_slide_v2/my_slide.pt \
    --weights model_weights/humanlung_cell2location_leave_A50_out.pth \
    --output  inference/my_slide_v2
```

### 3.5 단일 GPU / CPU 폴백

```bash
# GPU 1장만 강제
python inference/infer.py --world-size 1 --data ... --weights ... --output ...

# CPU 만 사용해야 한다면 — 매우 느림 (35k spot ≈ 1시간 이상). 권장 안함.
CUDA_VISIBLE_DEVICES= python inference/infer.py --data ... --weights ... --output ...
# 단, infer.py 는 CUDA 필수 가정 — CPU-only 는 코드 수정 필요
```

### 3.6 출력 (output 디렉토리에 추가됨)

```
inference/my_slide_v2/
├── ... (prep 산출물 그대로)
├── predictions.csv        # spot_id, X, Y, AT1, AT2, …, gdT  (N × 83)
└── predictions.npy        # float32 [N, 80] (csv 와 동일 row 순서)
```

### 3.7 진행 모니터링

`infer.py` 는 GPU 0 에서 약 800 spot 마다 진행률을 stdout 으로 찍는다. 백그라운드 실행 후 `tail -f` 로 추적:
```bash
python inference/infer.py ... 2>&1 | tee /tmp/infer.log &
tail -f /tmp/infer.log | grep -E "gpu0|done|Saved"
```

### 3.8 자주 발생하는 버그 대응 (이미 코드에 반영됨)

- **`Encountered invalid feature tensor type (got 'list')`** → worker 가 자동으로 `data.spot_id` 를 제거 (코드 반영 완료)
- **`Input should be contiguous`** → prep 의 edge_index 가 contiguous 로 저장 (코드 반영 완료). 이전 .pt 면 한 번 `torch.load` → `data.edge_index = data.edge_index.contiguous()` → `torch.save` 로 patch 가능
- **`PyG 2.7 input_id is local`** → infer 가 `shard[input_id]` 로 글로벌 인덱스 매핑 (코드 반영 완료)

---

## 4. STEP 3 — Analysis (`inference/analysis/analyze.py`)

`predictions.csv` + `coords.h5` + `cell_type_groups.csv` → 그룹별 abundance + spatial heatmap + 80×80 Moran's R clustermap.

### 4.1 풀 옵션

```text
--predictions PREDICTIONS      (필수) STEP 2 의 predictions.csv
--coords COORDS                (필수) STEP 1 의 *_coords.h5
--groups GROUPS                (필수) cell type → lineage 분류 CSV
--output OUTPUT                (필수) 출력 디렉토리
--knn K                        Moran's R weight matrix 의 k. 디폴트 20
```

### 4.2 `--groups` (cell type 분류 CSV) 형식

`inference/analysis/cell_type_groups.csv` 가 표준. 헤더:

```
cell_type,group,is_strict_proxy,is_broad_proxy,note
Basal,Epithelial-airway,1,1,airway basal stem cell (strongest cross-tissue marker)
AT2,Epithelial-alveolar,0,1,alveolar type 2 (broad-only — verification hypothesis)
Dividing_AT2,Epithelial-alveolar,1,1,explicitly dividing AT2 (proliferative-like signal)
…
```

자세한 설계 근거는 `inference/analysis/EPITHELIAL_PROXY_METHODOLOGY.md` (strict 3종 + broad 5종 의 cross-tissue 신뢰도 + reference 13 개).

총 80 row, 모든 cell type 의 `cell_type` 이 정확히 `cell_types.pkl` 의 80 이름과 set 일치 해야 함 (`analyze.py` 가 시작 시 검증).

**다른 조직 / 다른 가중치를 쓰면** 본 80 cell type 분류는 의미가 달라진다. 자기 케이스에 맞게 group / epithelial-activity proxy flags 직접 정의 권장.

### 4.3 일반 케이스

```bash
python inference/analysis/analyze.py \
    --predictions inference/my_slide_v2/predictions.csv \
    --coords      inference/my_slide_v2/my_slide_coords.h5 \
    --groups      inference/analysis/cell_type_groups.csv \
    --output      inference/analysis/my_slide_v2
```

### 4.4 `--knn` 가이드 (Moran's R weight)

| spot 수 | 권장 `--knn` |
|---|---:|
| < 5,000 | `8` |
| 5,000 – 20,000 | `12` |
| 20,000 – 50,000 | `20` (디폴트) |
| 50,000+ | `30` |

너무 작으면 spatial scale 이 짧아져 noise. 너무 크면 distant spot 까지 weight 가 가서 R 이 평탄.

### 4.5 출력

```
inference/analysis/my_slide_v2/
├── abundance_by_celltype.csv    # 80 type 별 mean/median/max/fraction-nonzero
├── abundance_by_group.csv       # 그룹 합 + strict / broad epithelial-activity proxy
├── spatial_top10_celltypes.png  # 평균 상위 10 type spot scatter
├── spatial_group_heatmaps.png   # 10 그룹 spatial sum panel
├── spatial_immune_vs_epithelial.png # 1×3 panel: immune / strict proxy / broad proxy
├── moran_r_pairs.csv            # cell-pair (3,240) Moran's R, z, p
└── moran_r_clustermap.png       # 80×80 hierarchical clustermap
```

소요 시간: 슬라이드 당 30–60초 (CPU). 가장 무거운 부분은 80×80 Moran R sparse 곱.

---

## 5. End-to-end one-shot

세 명령을 chain. 새 슬라이드가 도착하면 다음 셀을 슬라이드 path / mpp 만 바꿔 실행.

```bash
SLIDE=/mnt/path/to/MY_SLIDE.svs
NAME=my_slide
OUT=inference/${NAME}_v2
ANA=inference/analysis/${NAME}_v2

# (선택) mpp 확인
python -c "import openslide; s=openslide.OpenSlide('$SLIDE'); \
  print('dims', s.dimensions, 'mpp_x', s.properties.get('openslide.mpp-x'))"

# 1. PREP
python prep/prepare_wsi_for_inference_v2.py \
    --input "$SLIDE" --output "$OUT" --slide-name "$NAME" \
    --tile-size 400 --min-tissue-frac 0.40

# 2. spot_view.jpg 확인 후 진행 — 라벨 잡혔거나 tissue 빠졌으면 위 명령 다시
echo "open $OUT/spot_view.jpg ; press ENTER if OK" && read

# 3. INFERENCE (4-GPU)
python inference/infer.py \
    --data    "$OUT/${NAME}.pt" \
    --weights model_weights/humanlung_cell2location_leave_A50_out.pth \
    --output  "$OUT"

# 4. ANALYSIS
python inference/analysis/analyze.py \
    --predictions "$OUT/predictions.csv" \
    --coords      "$OUT/${NAME}_coords.h5" \
    --groups      inference/analysis/cell_type_groups.csv \
    --output      "$ANA"

# 5. 결과 확인
ls "$OUT" "$ANA"
```

---

## 6. 케이스 별 옵션 조합 (recipes)

### 6.1 일반 H&E 슬라이드 (40×, 깔끔)

```bash
prep:    --tile-size 400 --min-tissue-frac 0.30
infer:   default
analyze: default
```

### 6.2 라벨 sticker / inkstain 많은 슬라이드

```bash
prep:    --tile-size 400 --min-tissue-frac 0.50 --thumb-max-side 6000
# 더 정밀한 mask 를 위해 thumb 키움. spot 줄어드는 만큼 false positive 도 줄어듦
infer:   default
analyze: default
```

### 6.3 20× 슬라이드 (mpp ≈ 0.5)

```bash
prep:    --tile-size 200 --min-tissue-frac 0.30
# 학습 분포에 가까움. spot 수도 절반 수준으로 줄어 inference 빠름
```

### 6.4 매우 큰 WSI (level-0 > 200k × 100k)

```bash
prep:    --tile-size 600 --min-tissue-frac 0.40 --thumb-max-side 6000
# spot 수를 적정 수준 (10k–30k) 로 제어. memory 안전
infer:   --batch-size 16 (24GB GPU 기준)
```

### 6.5 작은 GPU (8GB) 환경

```bash
infer:   --batch-size 4 --world-size 1
# 단일 GPU + 작은 batch. 시간 5-10× 길어짐
```

### 6.6 sanity-check (작게 빨리)

```bash
prep:    --tile-size 800 --min-tissue-frac 0.50
# spot 1-2k 정도로 줄여 모든 단계 1-2분 안에 완료
```

### 6.7 super-resolution (고해상도 spatial map)

```bash
prep:    --tile-size 200 --min-tissue-frac 0.20
# 격자 절반, threshold 낮춰 spot 4× 증가. inference 시간/메모리 4× 부담
infer:   default (GPU 충분히)
```

---

## 7. 자기 분석 코드 작성용 cheat sheet

`predictions.csv / .npy / coords.h5` 만 있으면 본 `analyze.py` 없이 임의 분석 가능:

```python
import pandas as pd, numpy as np, h5py

# 1) 핵심 입력 4개
preds_df = pd.read_csv("inference/my_slide_v2/predictions.csv")
preds = np.load("inference/my_slide_v2/predictions.npy")          # 같은 데이터 npy 형식
with h5py.File("inference/my_slide_v2/my_slide_coords.h5") as f:
    coords = f["coords"][:]                                          # (N, 2) tile top-left
    meta = dict(f["metadata"].attrs.items())                         # mpp, tile_size, etc.
groups = pd.read_csv("inference/analysis/cell_type_groups.csv")

cell_cols = [c for c in preds_df.columns if c not in ("spot_id", "X", "Y")]
assert len(cell_cols) == 80
assert preds.shape == (len(preds_df), 80)

# 2) 인덱스 정렬 — row i 가 같은 spot
#    preds_df.iloc[i] ↔ preds[i] ↔ coords[i]

# 3) 그룹 합산 예시
def group_sum(g_name):
    members = groups.loc[groups["group"] == g_name, "cell_type"].tolist()
    return preds_df[members].sum(axis=1)

immune = group_sum("Immune-lymphoid") + group_sum("Immune-myeloid")
strict_proxy = preds_df[groups.loc[groups["is_strict_proxy"]==1, "cell_type"]].sum(axis=1)
broad_proxy  = preds_df[groups.loc[groups["is_broad_proxy"]==1,  "cell_type"]].sum(axis=1)

# 4) spatial 시각화 (matplotlib)
import matplotlib.pyplot as plt
plt.scatter(preds_df["X"], preds_df["Y"], c=immune, s=1, cmap="viridis")
plt.gca().invert_yaxis(); plt.colorbar(label="Immune total"); plt.show()

# 5) 224x224 패치 재추출 (모델 입력 그대로)
import openslide
sl = openslide.OpenSlide("/mnt/path/MY_SLIDE.svs")
i = 100
X, Y = int(preds_df.loc[i, "X"]), int(preds_df.loc[i, "Y"])
patch = sl.read_region((X-112, Y-112), 0, (224, 224)).convert("RGB")
sl.close()
```

---

## 8. Troubleshooting

| 증상 | 원인 / 조치 |
|---|---|
| `ModuleNotFoundError: openslide` | `VIRTUAL_ENV=.venv uv pip install openslide-bin openslide-python` |
| prep 결과 spot 수가 비정상적 (수십 / 수십만) | `spot_view.jpg` 확인 → `--tile-size` 또는 `--min-tissue-frac` 조정 |
| `cuda: True` 인데 inference 시작 직후 OOM | `--batch-size` 절반으로 ↓ 또는 `--world-size 1` 로 변경 후 batch 작게 |
| `Encountered invalid feature tensor type` | 옛 prep 결과 사용 중 — 최신 prep_v2 로 다시 돌리거나 .pt 에서 `data.spot_id` 삭제 |
| `Input should be contiguous` | `data.edge_index = data.edge_index.contiguous(); torch.save(data, …)` 로 한 번 patch |
| `analyze.py` 가 "groups CSV is missing N cell types" 에러 | `cell_type_groups.csv` 의 cell_type 컬럼이 cell_types.pkl 과 set 불일치 — 80개 모두 정확히 |
| Moran R clustermap 이 거의 균일 | `--knn` 너무 큼 — 절반으로 |
| Moran R clustermap 너무 noisy | `--knn` 너무 작음 — 두 배로 |
| GitHub push 가 50MB warning | 분석 산출물 (.pt) 은 `.gitignore` 의 `**/*.pt` 로 차단되어 있음. 다른 큰 파일이면 LFS 고려 |
| ssh push 시 DNS 일시 실패 | `GIT_SSH_COMMAND="ssh -o ConnectTimeout=10" git push origin main` |

---

## 9. 산출물 → 동료 공유 체크리스트

| 받는 쪽 | 필요 파일 | 위치 |
|---|---|---|
| spatial 분석가 | `predictions.csv`, `*_coords.h5`, `spots.csv`, `cell_type_groups.csv`, `analyze.py` | `inference/<slide>_v2/`, `inference/analysis/` |
| proteomics 매칭 | 위 + `tissue_mask.png`, `spot_view.jpg` (좌표계 검증), `inference/analysis/README.md` (mpp/매칭 워크플로) | 같음 |
| 빠른 시각 검토 | `spatial_top10_celltypes.png`, `spatial_group_heatmaps.png`, `spatial_immune_vs_epithelial.png`, `moran_r_clustermap.png` | `inference/analysis/<slide>_v2/` |
| 임상 review | `findings.md` + 위 PNG 4장 | 같은 폴더 |

git push 시 `.pt` (수십 GB) 는 자동 차단, 나머지 (predictions / heatmap / Moran) 만 올라간다 (`.gitignore` 규칙 참고).

---

## 10. 관련 문서

- 이 cookbook: `report/04_WSI에서_분석까지_쿡북.md`
- prep v1 → v2 배경 / 의도: `report/02_추론용_WSI전처리_prep스크립트_사용법및절차.md`, `report/03_breast슬라이드2장_lung가중치_추론결과_v2framework.md`
- 학습 prep 의 외부 파이프라인 (cell2location 등): `report/01_원본_WSI전처리파이프라인_튜토리얼이전단계_요약.md`
- 분석 결과 해석 가이드 (proteomics 매칭, caveat, KBSMC cohort, TCGA 검증): `inference/analysis/README.md`
- 슬라이드 별 통합 소견 예시: `inference/analysis/slide{1,2}_*_v2/findings.md`
- 실제 코드: `prep/prepare_wsi_for_inference_v2.py`, `inference/infer.py`, `inference/analysis/analyze.py`

---

## 11. 한 번 더 caveat

본 파이프라인은 **lung-trained Hist2Cell 가중치** 를 기준으로 한다. breast / 다른 조직에 적용 시:
- cell type 이름은 lung 분류 — 그룹 단위 / 공간 패턴 / 상대 비교만 신뢰
- epithelial-activity proxy: strict (Basal+Dividing_AT2+Dividing_Basal, 3 종) vs broad (위 + AT2+Suprabasal, 5 종). 종양 직접 검출 아님 — `inference/analysis/EPITHELIAL_PROXY_METHODOLOGY.md` 의 cross-tissue 신뢰도 표 필독
- Visium 학습 분포 ~0.5 mpp vs Aperio 40× 0.26 mpp: 모델 시야 절반 수준 — 절대값보다 패턴
- proteomics 등 다른 modality 와 **공간 일치 검증 후** 정량 결론 도출 권장
