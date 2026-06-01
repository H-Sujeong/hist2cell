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

## 224 vs 146 — FOV(배율) 매칭
20× 슬라이드(mpp 0.5015)에서 224 px = **112 µm**.
HEX 학습 FOV = 224 px × **0.325 µm/px = 72.8 µm** → 0.5015 슬라이드에선 **146 px**(=73.2µm) crop 이 맞다.
(동료 타일링이 146 px·overlap 0 으로 이미 일치.)

| 타일 | 물리 FOV | resize 후 eff. mpp | 대상 모델 | 패치 처리 |
|---|---|---|---|---|
| **224** | 112 µm | 0.5015 (native) | Hist2Cell (20× 학습) | 224 px 그대로 |
| **146** | 73.2 µm | **0.327 ≈ HEX 0.325** | HEX (학습 mpp 0.325) | 146 px crop → 224 로 ×1.53 resize |

> **정정 (2026-05-29)**: HEX 학습 mpp 가 **0.325** 로 확인 → 기존 112px(56µm·40×0.25 가정) 폐기.
> 기존 `graph_output/112`·`tilitng_output/112` 삭제, 동료 146 타일링으로 `graph_output/146` 재생성.
> 슬라이드 native 0.5015 > 0.325 라 146→224 는 보간 업샘플 — FOV 는 맞지만 해상도는 여전히 OOD.

## 파이프라인 상태

| 단계 | 상태 | 위치 |
|---|---|---|
| ① Tiling | ✅ 완료 | `tilitng_output/224/TCGA-LUAD/`, `tilitng_output/146/` (HEX, 동료 제공) |
| ② Graph (`.pt`) | ✅ 완료 | `graph_output/224/`, `graph_output/146/` (2026-05-29 FOV 정정) |
| ③ Inference | ✅ 완료 (2026-05-24, predictions + features) | `inference_output/<slide>/predictions.{csv,npy}` + `features_resnet.npy` + `features_fused.npy` |
| ④ Hist2Cell UMAP baseline | ✅ 완료 (2026-05-24) | `umap_output/` (초기 4 PNG + `summary.md`) |
| ⑤ DINOv2 ViT-B/14 추론 | ✅ 완료 (2026-05-25) | `dino_output/<slide>/features_dinov2.npy` [N,768] |
| ⑥ UMAP 4 rep 비교 (Hist2Cell × 3 + DINOv2) | ✅ 완료 (2026-05-25) | `umap_output/` PNG 재생성 + `summary.md` 갱신 |
| ⑦ 정량 비교 (Hist2Cell vs DINOv2) | ✅ 완료 (2026-05-25) | `compare_output/` (`metrics.csv` + `metrics_bars.png` + `summary.md`) |
| ⑧ 3×4 UMAP grid + Epi/Stro subtype 색칠 | ✅ 완료 (2026-05-26) | `umap_output/per_slide_grid_{lineage,epithelial,stromal}.png` + `embeddings/` cache |
| ⑨ Batch metric 정정 (size-weighted chance + balanced subsample UMAP) | ✅ 완료 (2026-05-26) | `compare_output/metrics_corrected.csv` + `umap_output/cross_slide_balanced.png` |
| ⑩ Slide-별 TOP10 cell type 통계 + UMAP overlay | ✅ 완료 (2026-05-26) | `top10_output/` (`top10_stats.csv` + `top10_union.csv` + per-slide PNG × 3 + `summary.md`) |
| ⑪ HEX FOV 정정(146px) + 146-grid Hist2Cell·DINOv2·UMAP·TOP10 | ✅ 완료 (2026-05-29) | `graph_output/146/` + `inference_output_146/` + `dino_output_146/` + `umap_output_146/` (4 PNG + `summary.md`) + `top10_output_146/` (stats/union CSV + 3 PNG + `summary.md`) |
| ⑫ DINO cluster(=dominant ct) 패치 grid + 4-rep UMAP overlay (224·146) | ✅ 완료 (2026-05-29) | `dino_cluster_output/{224,146}/` (umap_4rep_by_dominant_ct.png + cluster_NN_*.png + dino_clusters_*.csv) + `summary.md` |
| ⑬ hex+dino vs dino vs prediction 3×3 UMAP + kNN purity (224·146) | ✅ 완료 (2026-06-01) | `hex_compare_{224,146}/` (umap_3x3 + knn_purity{,_weighting}.csv + summary.md). 결과 **해상도 의존**: 224 미지지, 146(OOD) 2/3 슬라이드 부분 지지. TS1 agg 버그 정정 확인 |
| ⑭ dominant cell type 별 대표 패치 montage + 조직학 가이드 (224·146) | ✅ 완료 (2026-06-01) | `celltype_examples_{224,146}/` (montage PNG + csv + summary). prediction centroid 최근접 예시 |
| ⑮ Visium human lung GT 셋 3장 + DINO 추론 (hex 결정평가용) | ✅ 완료 (2026-06-01) | `visium_gt/` (dino_output/ + hist2cell_output/ + gt/ 80celltype + README) |
| ⑯ GT vs Hist2Cell vs DINO 비교 (실제 GT 기준, circular 아님) | ✅ 완료 (2026-06-01) | `visium_gt/compare/` (purity/accuracy csv + umap + summary). Hist2Cell≈GT, DINO 는 gap −0.12~−0.29 → hex 채울 여지 실재. HEX feature(FOV 재crop) 대기 |
| ⑰ dino→19dim(PCA/중요도) 축소 후 hex 비교 (224·146) | ✅ 완료 (2026-06-01) | `hex_dino19_{224,146}/` (purity csv + umap + summary). 차원 맞추니 hex 가 dino 와 대등(224)/우세(146) — 이전 "hex 무의미" 는 차원지배 아티팩트 |
| ⑱ hex19 패치별 하위10%→0(denoise) 후 비교 (224·146) | ✅ 완료 (2026-06-01) | `hex_dino19_denoise_{224,146}/`. purity 무영향(≤0.004), UMAP 모양만 호→분절(모양≠정보 실증) |
| ⑲ hex19 패치당 상위50%만 유지(하위50%→0) 비교 (224·146) | ✅ 완료 (2026-06-01) | `hex_dino19_top50_{224,146}/`. hex 단독 약손해, 146 pca+hex 만 소폭↑(≤0.009) — 자를수록 손해 쪽 |

## 폴더 구조
```
lung_pilot/
├── tilitng_output/
│   ├── 224/TCGA-LUAD/   # 224 타일: <slide>.h5 + Thumbnails/Masks/Overlays + tiling_summary.md
│   └── 146/             # 146 타일 (HEX·73.2µm, 동료 제공): <slide>.h5 + Overlays
├── graph_output/
│   ├── 224/             # <slide>.pt (Hist2Cell 입력) + <slide>_spots.csv
│   ├── 146/             # <slide>.pt (HEX 입력, 146→224 resize) + <slide>_spots.csv
│   └── README.md        # .pt 포맷·로드법·주의
├── inference_output/    # Hist2Cell 추론 결과 (224 grid) — <slide>/{predictions.{csv,npy}, features_resnet.npy [N,512], features_fused.npy [N,256]} + _logs/
├── inference_output_146/ # Hist2Cell 추론 (146 grid, HEX FOV 73.2µm·OOD) — 동일 구조
├── dino_output_146/     # DINOv2 (146 grid) — <slide>/features_dinov2.npy [N,768]
├── umap_output_146/     # 146-grid 4-rep UMAP (per_slide ×3 + cross_slide) + embeddings/ + summary.md
├── top10_output_146/    # 146-grid Slide-별 TOP10 (stats/union CSV + 3 PNG + summary.md)
├── dino_cluster_patches.py # DINO cluster(=dominant ct) centroid 최근접 패치 grid + 4-rep UMAP overlay
├── dino_cluster_output/  # {224,146}/ umap_4rep_by_dominant_ct.png + cluster_NN_*.png + dino_clusters_*.csv + summary.md
├── hex_compare.py        # prediction vs dino vs hex+dino 3×3 UMAP + kNN purity (경로 인자형)
├── hex_compare_{224,146}/ # 224·146 비교 (umap_3x3 + knn_purity{,_weighting}.csv + summary.md)
│   # 입력 dino/agg = /mnt/fileserver/lung_pilot/{dino_output,dino_hex_agg}{,_146} (agg=dino768⊕hex19)
├── celltype_examples.py  # dominant cell type 별 prediction centroid 최근접 대표 패치 (경로 인자형)
├── celltype_examples_224/ # 18 cell type 예시 montage + csv + summary(폐 조직학 가이드)
├── celltype_examples_146/ # 15 cell type 예시 montage + csv + summary(224 대비 FOV 차이)
├── visium_gt_compare.py  # Visium GT vs Hist2Cell vs DINO 비교 (실제 GT 기준)
├── visium_gt/            # Visium human lung GT 평가셋 (README + gt/ + dino_output/ + hist2cell_output/ + compare/)
├── hex_dino19_compare.py # dino 768→19(PCA/중요도) 축소 후 hex 비교 (차원 지배 보수)
├── hex_dino19_{224,146}/ # 224·146 결과 (knn_purity{,_pivot}.csv + umap + summary)
├── hex_dino19_denoise.py # hex19 패치별 하위 N%→0 denoise 변형 (--pct)
├── hex_dino19_denoise_{224,146}/ # 하위10%→0 결과 (purity 무영향, 모양만 변화)
├── hex_dino19_top50_{224,146}/   # 상위50%만 유지(--pct 50) 결과 (hex 단독 약손해)
│   # HEX 추출 코드: hex_inference.ipynb + model_hex_compgat_clpg_cv.py (입력=Optimus 1536-d → 19 marker)
├── dino_infer.py        # DINOv2 ViT-B/14 추론 (외부 /home/sjhong/dinov2 import + 가중치 절대경로)
├── dino_output/         # DINOv2 추론 결과 — <slide>/features_dinov2.npy [N,768] + _logs/
├── umap_compare.py      # 4 rep × 3 slide UMAP 시각화 (1×4 per-slide + cross-slide)
├── umap_subtype_grid.py # 3×4 grid + Epi/Stro subtype 색칠 (cached embedding 재사용)
├── umap_output/         # UMAP PNG (per_slide ×3, cross_slide, grid ×3) + embeddings/ + summary.md
├── compare_hist2cell_vs_dinov2.py  # 정량 metric (1-NN purity / kNN overlap / silhouette)
├── compare_output/      # metrics.csv + metrics_bars.png + summary.md + metrics_corrected.csv
├── batch_recheck.py     # size-weighted chance + balanced subsample UMAP
├── top10_umap.py        # slide-별 TOP10 (mean abundance desc) 통계 + UMAP overlay
└── top10_output/        # top10_stats.csv + top10_union.csv + per-slide PNG × 3 + summary.md
```
세부 문서: 각 `tilitng_output/*/tiling_summary.md`, `graph_output/README.md`.

## 타일 / 노드 수

| 슬라이드 | 224 (Hist2Cell) | 146 (HEX, 동료 grid) |
|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 | 6,020 |
| TCGA-05-4245-01A-01-TS1 | 1,871 | 4,257 |
| TCGA-05-4390-01A-01-BS1 | 10,661 | 24,462 |

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
| 1-NN purity (전체 15,401, excess over weighted chance 0.529) | 0.775 (+0.523) | 0.947 (+0.888) | 0.813 (+0.604) | 0.908 (+0.805) |
| balanced 1-NN purity (1,871×3, excess over 0.333) | **0.663 (+0.495)** | 0.897 (+0.846) | 0.711 (+0.567) | 0.836 (+0.754) |
| kNN overlap vs DINOv2 (k=10 / k=50, 높을수록 유사) | 0.005 / 0.012 | **0.014 / 0.021** | 0.007 / 0.016 | — |
| silhouette by dominant lineage (높을수록 cell-type 분리) | **0.017** | -0.000 | 0.011 | 0.013 |

**핵심 교훈** (2026-05-26 정정 후) — chance baseline 을 정확히 잡고
balanced subsample 으로 spot-수 효과를 제거하면, **prediction 이 명확히
가장 batch-free** (balanced excess +0.495, 다른 rep 의 절반 정도). 다른
rep 의 순위 (fused > dinov2 > resnet > prediction) 는 유지. 특히
`features_fused` 의 매우 높은 purity 는 GAT graph aggregation 이 같은
슬라이드 spot 의 feature 를 평균내는 구조적 강제 — cell-type 신호와
무관. 자세한 해석: **`compare_output/summary.md`** + balanced UMAP:
**`umap_output/cross_slide_balanced.png`**.

**Hist2Cell ↔ DINOv2 representation 의 이웃 구조** 는 chance 보다 10–30배
높지만 절대값으로는 매우 낮음 (Jaccard ≤ 0.021) — 두 모델이 같은 spot
의 "유사 spot" 을 거의 다른 신호로 정의한다.

## 다음 단계

1. **HEX expression** (외부, 사용자 측) — `graph_output/146/*.pt` 입력 (FOV 73.2µm, eff. mpp 0.327 ≈ HEX 0.325).
   도착하면 4 rep UMAP 와 정량 비교 framework 에 5번째 rep 으로 추가.
   본 figure 의 `features_dinov2` 가 HEX 의 DINO 블록 baseline 역할.
2. (선택) **정확한 chance baseline** — slide-size 가중 1-NN purity baseline
   계산해 비교를 정량적으로 sharper 하게.
3. (선택) **DINOv2 다른 사이즈** (ViT-S, ViT-L) 비교.
4. (선택) **abundance vector 직접 거리** — argmax label silhouette 대신
   abundance 자체의 inter-spot Spearman correlation 등.

## 주의
- TCGA-LUAD 는 native 20× (mpp 0.5015) — 146→224 업샘플 패치는 FOV(73.2µm)는 HEX 학습과 맞으나
  해상도는 여전히 OOD (0.5015→0.327 보간). HEX 출력·UMAP 해석 시 명시할 것.
- `.pt` 로드: `torch.load(..., weights_only=False)` (spot_id 가 python list).
  `graph_output/146/TCGA-05-4390-01A-01-BS1.pt` 는 14.7 GB — 메모리 여유 확보.
- 224 ↔ 146 그리드 center 불일치 (224 중심 = coord+112, 146 중심 = coord+73).
  146 grid 는 동료 타일링과 1:1 정합(centers == coords+73 검증) → HEX expression spot 대응.
  224↔146 spot 단위 paired 비교가 필요하면 공유 center 에서 두 크기 패치를 추출해야 한다.
- concat→UMAP 시 HEX·DINO 블록은 스케일·차원이 다르므로 **블록별 정규화/PCA** 후 concat.
  실측: agg=[dino768 ⊕ hex19], hex 블록 per-dim std 가 dino 의 ~100배(157~5625) → raw 면 hex 가 거리 지배.
  `hex_compare.py` 는 per-dim z-score 적용. (hex_compare_224/summary.md 참고)
- **fileserver agg 버그 (정정됨, 2026-06-01)**: `dino_hex_agg_146/TCGA-05-4245-01A-01-TS1/features_agg.npy`
  가 N=24462(4390 것)였으나 4257 로 재생성 완료 → 146 hex+dino 비교 수행됨. (224 는 처음부터 정상.)
- 코드: 타일링 `WSI_tile_sampling_framework/run_tiling_tcga_luad.py`,
  그래프 `prep/build_graph_from_tiles.py`, 추론 `inference/infer.py`.
  타일링 how-to 는 `report/05_*`.
