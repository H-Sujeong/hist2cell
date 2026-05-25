# lung_pilot — TCGA-LUAD 3장 cell-type 표현 비교 pilot

## 목적
TCGA-LUAD H&E 슬라이드 3장에 **두 모델**을 돌려 cell-type 표현을 UMAP 으로 비교한다.

- **Hist2Cell** (20× 학습) — cell-type abundance 직접 예측
- **HEX 모델** (40×·224 학습, repo 외부 모델) — expression 예측 → DINO 벡터와 concat
- 비교: `[HEX expression ⊕ DINO]` UMAP  vs  `Hist2Cell cell-type` UMAP

## 슬라이드 (3장, 모두 20×, mpp ≈ 0.502)
원본 위치: `/mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi/`

| short name | SVS 파일 |
|---|---|
| TCGA-05-4245-01A-01-BS1 | `...-BS1.41d3cf23-4e36-4e42-9e08-adfea139f37e.svs` |
| TCGA-05-4245-01A-01-TS1 | `...-TS1.bf71c76b-e802-4a7a-b6c3-c5f46212fab0.svs` |
| TCGA-05-4390-01A-01-BS1 | `...-BS1.38f2a7ef-442a-4fa6-acad-6e5d567bdcfd.svs` |

## 224 vs 112 — 배율 매칭
20× 슬라이드에서 224 px = **112 µm**, 112 px = **56 µm**.

| 타일 | 물리 크기 | 대상 모델 | 패치 처리 |
|---|---|---|---|
| **224** | 112 µm | Hist2Cell (20× 학습) | 224 px 그대로 |
| **112** | 56 µm | HEX (40×·224 학습) | 112 px crop → 224 로 ×2 resize |

## 파이프라인 상태

| 단계 | 상태 | 위치 |
|---|---|---|
| ① Tiling | ✅ 완료 | `tilitng_output/224/TCGA-LUAD/`, `tilitng_output/112/` |
| ② Graph (`.pt`) | ✅ 완료 | `graph_output/224/`, `graph_output/112/` |
| ③ Inference | ✅ 완료 (2026-05-24, predictions + features) | `inference_output/<slide>/predictions.{csv,npy}` + `features_resnet.npy` + `features_fused.npy` |
| ④ Hist2Cell UMAP baseline | ✅ 완료 (2026-05-24) | `umap_output/` (초기 4 PNG + `summary.md`) |
| ⑤ DINOv2 ViT-B/14 추론 | ✅ 완료 (2026-05-25) | `dino_output/<slide>/features_dinov2.npy` [N,768] |
| ⑥ UMAP 4 rep 비교 (Hist2Cell × 3 + DINOv2) | ✅ 완료 (2026-05-25) | `umap_output/` PNG 재생성 + `summary.md` 갱신 |
| ⑦ 정량 비교 (Hist2Cell vs DINOv2) | ✅ 완료 (2026-05-25) | `compare_output/` (`metrics.csv` + `metrics_bars.png` + `summary.md`) |

## 폴더 구조
```
lung_pilot/
├── tilitng_output/
│   ├── 224/TCGA-LUAD/   # 224 타일: <slide>.h5 + Thumbnails/Masks/Overlays + tiling_summary.md
│   └── 112/             # 112 타일: 동일 구조 + tiling_summary.md
├── graph_output/
│   ├── 224/             # <slide>.pt (Hist2Cell 입력) + <slide>_spots.csv
│   ├── 112/             # <slide>.pt (HEX 입력) + <slide>_spots.csv
│   └── README.md        # .pt 포맷·로드법·주의
├── inference_output/    # Hist2Cell 추론 결과 — <slide>/{predictions.{csv,npy}, features_resnet.npy [N,512], features_fused.npy [N,256]} + _logs/
├── dino_infer.py        # DINOv2 ViT-B/14 추론 (외부 /home/sjhong/dinov2 import + 가중치 절대경로)
├── dino_output/         # DINOv2 추론 결과 — <slide>/features_dinov2.npy [N,768] + _logs/
├── umap_compare.py      # 4 rep × 3 slide UMAP 시각화 스크립트 (Hist2Cell 3 + DINOv2)
├── umap_output/         # UMAP PNG 4장 + summary.md (해석)
├── compare_hist2cell_vs_dinov2.py  # 정량 metric (1-NN purity / kNN overlap / silhouette)
└── compare_output/      # metrics.csv + metrics_bars.png + summary.md
```
세부 문서: 각 `tilitng_output/*/tiling_summary.md`, `graph_output/README.md`.

## 타일 / 노드 수

| 슬라이드 | 224 | 112 |
|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 | 11,470 |
| TCGA-05-4245-01A-01-TS1 | 1,871 | 7,499 |
| TCGA-05-4390-01A-01-BS1 | 10,661 | 42,615 |

## 추론 결과 (2026-05-24 완료, features 포함 재추론)

| slide | spots | predictions (80-d) | features_resnet (512-d) | features_fused (256-d) | shard 시간 (4× A5000) |
|---|---|---|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 | row_sum 1.32–62.81 | NaN 0, range [0, 4.5] | NaN 0, range ≈[-5, 5] | 17.4 s |
| TCGA-05-4245-01A-01-TS1 | 1,871 | row_sum 1.85–49.43 | NaN 0, range [0, 4.2] | NaN 0, range ≈[-5, 4] | 21.6 s |
| TCGA-05-4390-01A-01-BS1 | 10,661 | row_sum 1.70–49.28 | NaN 0, range [0, 4.5] | NaN 0, range ≈[-5, 5] | 55.4 s |

세 가지 산출:
- **`predictions.{csv,npy}`** — 80 cell-type abundance (cell2location scale, 확률 아님). CSV = `spot_id, X, Y, <80 cell types>` (83 col), NPY = `(N, 80)`.
- **`features_resnet.npy`** — `(N, 512)`. ResNet18 backbone 의 spot 단위 visual feature (graph context 미포함, ReLU 후 비음수). HEX 의 DINO 와 직접 대응되는 모폴로지 representation.
- **`features_fused.npy`** — `(N, 256)`. `(x_spot_e + x_local + x_global)/3` — Hist2Cell 의 fused_head 직전 통합 representation (visual + GAT graph + transformer). prediction 의 직접 precursor.

추론 명령 (참고):
```bash
cd /home/sjhong/hist2cell
for s in TCGA-05-4245-01A-01-BS1 TCGA-05-4245-01A-01-TS1 TCGA-05-4390-01A-01-BS1; do
  .venv/bin/python inference/infer.py \
    --data    lung_pilot/graph_output/224/$s.pt \
    --weights model_weights/humanlung_cell2location_leave_A50_out.pth \
    --output  lung_pilot/inference_output/$s
done
```
(`infer.py` 의 `Hist2Cell.forward(..., return_features=True)` 가 두 feature 를 함께 반환,
worker 가 shard 에 저장 → main 이 `features_{resnet,fused}.npy` 로 머지.)

## ⑤ DINOv2 ViT-B/14 — 외부 self-supervised baseline (2026-05-25)

같은 224 patch 에 DINOv2 ViT-B/14 (CLS, 768-d) 를 통과시켜 외부 모델의
visual representation 도 비교 대상에 포함. 외부 코드/가중치 (**git ignore
처리, 직접 다운로드 필요**):

- 코드: `external/dinov2/` (= facebookresearch/dinov2 clone; `dino_infer.py` 의 default)
  ```bash
  git clone https://github.com/facebookresearch/dinov2.git external/dinov2
  # 또는 zip:
  # wget https://github.com/facebookresearch/dinov2/archive/refs/heads/main.zip -O /tmp/dinov2.zip
  # unzip -q /tmp/dinov2.zip -d external/ && mv external/dinov2-main external/dinov2
  ```
- 가중치: `/home/sjhong/dinov2_vitb14_pretrain.pth` (~330 MB, 절대경로 default)
  ```bash
  wget https://dl.fbaipublicfiles.com/dinov2/dinov2_vitb14/dinov2_vitb14_pretrain.pth \
    -O /home/sjhong/dinov2_vitb14_pretrain.pth
  ```
- xformers 없는 환경에서 fallback attention 으로 동작 (warning 만)
- 4× A5000 DataParallel, 3 슬라이드 합쳐 ~25 s

추론 명령 (참고):
```bash
cd /home/sjhong/hist2cell
for s in TCGA-05-4245-01A-01-BS1 TCGA-05-4245-01A-01-TS1 TCGA-05-4390-01A-01-BS1; do
  .venv/bin/python lung_pilot/dino_infer.py \
    --data   lung_pilot/graph_output/224/$s.pt \
    --output lung_pilot/dino_output/$s
done
```

## ⑥ UMAP 4 rep 비교 (2026-05-25 완료)

Hist2Cell 3 rep (`prediction_log1p` / `features_fused` / `features_resnet`)
+ DINOv2 (`features_dinov2`) = **총 4 representation** 으로 cross-slide
+ per-slide UMAP 재생성. 해석: **`umap_output/summary.md`**.

주요 관찰 (요약):

- **prediction_log1p** (Hist2Cell, 80-d) — 3 슬라이드 거의 완전 섞임 → tissue-generic.
- **features_fused** (Hist2Cell, 256-d) — 부분 batch.
- **features_resnet** (Hist2Cell, 512-d) — 슬라이드별 강한 분리. raw visual 의 batch.
- **features_dinov2** (DINOv2, 768-d) — Hist2Cell ResNet 보다는 약하지만 여전히
  부분 slide-separation. self-supervised 라 cell-type cluster 는 약하고,
  manifold 가 morphology/stain 축으로 정리됨.

→ HEX 결과가 도착하면 `prediction_log1p ↔ HEX expression`, `features_dinov2 ↔
HEX 측 DINO` 짝이 가장 자연스러운 대응.

## ⑦ 정량 비교 결과 요약 (2026-05-25)

| metric | prediction_log1p | features_fused | features_resnet | features_dinov2 |
|---|---|---|---|---|
| slide 1-NN purity (낮을수록 mix 좋음, chance≈0.33) | **0.775** | 0.947 | 0.813 | 0.908 |
| kNN overlap vs DINOv2 (k=10 / k=50, 높을수록 유사) | 0.005 / 0.012 | **0.014 / 0.021** | 0.007 / 0.016 | — |
| silhouette by dominant lineage (높을수록 cell-type 분리) | **0.017** | -0.000 | 0.011 | 0.013 |

**핵심 교훈** — UMAP 의 시각적 batch 판단과 raw 1-NN purity 가 *순위가
다르다*. 특히 `features_fused` 가 UMAP 에서는 중간 batch 였지만 raw
에서는 가장 batch-confined (0.947). 이유: GAT graph aggregation 이
같은 슬라이드 이웃 spot feature 를 평균내기 때문 — cell-type 신호와
무관한 구조적 강제. 자세한 해석: **`compare_output/summary.md`**.

**Hist2Cell ↔ DINOv2 representation 의 이웃 구조** 는 chance 보다 10–30배
높지만 절대값으로는 매우 낮음 (Jaccard ≤ 0.021) — 두 모델이 같은 spot
의 "유사 spot" 을 거의 다른 신호로 정의한다.

## 다음 단계

1. **HEX expression** (외부, 사용자 측) — `graph_output/112/*.pt` 입력.
   도착하면 4 rep UMAP 와 정량 비교 framework 에 5번째 rep 으로 추가.
   본 figure 의 `features_dinov2` 가 HEX 의 DINO 블록 baseline 역할.
2. (선택) **정확한 chance baseline** — slide-size 가중 1-NN purity baseline
   계산해 비교를 정량적으로 sharper 하게.
3. (선택) **DINOv2 다른 사이즈** (ViT-S, ViT-L) 비교.
4. (선택) **abundance vector 직접 거리** — argmax label silhouette 대신
   abundance 자체의 inter-spot Spearman correlation 등.

## 주의
- TCGA-LUAD 는 native 20× — 112→224 업샘플 패치는 면적(56µm)만 맞고 해상도는 OOD.
  HEX 출력·UMAP 해석 시 명시할 것.
- `.pt` 로드: `torch.load(..., weights_only=False)` (spot_id 가 python list).
  `graph_output/112/TCGA-05-4390-01A-01-BS1.pt` 는 24 GB — 메모리 여유 확보.
- 224 ↔ 112 그리드 center 불일치 (224 중심 = coord+112, 112 중심 = coord+56).
  spot 단위 paired 비교가 필요하면 공유 center 에서 두 크기 패치를 추출해야 한다.
- concat→UMAP 시 HEX·DINO 블록은 스케일·차원이 다르므로 **블록별 정규화/PCA** 후 concat.
- 코드: 타일링 `WSI_tile_sampling_framework/run_tiling_tcga_luad.py`,
  그래프 `prep/build_graph_from_tiles.py`, 추론 `inference/infer.py`.
  타일링 how-to 는 `report/05_*`.
