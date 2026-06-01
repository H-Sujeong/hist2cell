# dino→19dim 축소 후 hex 비교 (224) — 차원 지배 문제 보수

생성: 2026-06-01. 스크립트: `lung_pilot/hex_dino19_compare.py`.
이전 `hex_compare_224` 에서 agg=dino768⊕hex19 의 per-dim concat 은 dino 768차원이 hex(19=2.4%)를
압도해 hex 효과가 묻혔다. 여기선 **dino 를 19-dim 으로 줄여 hex 와 동등 차원**으로 맞춰 재비교.

## dino 19-dim 축소 2가지
- **PCA**: dino 768 → 19 주성분 (unsupervised). 19 PC 가 dino 분산의 **68.2%** 포착.
- **IMP**: dominant cell type 에 대한 ANOVA F-score 상위 19 원본 dim (supervised).
  ⚠️ 라벨=Hist2Cell prediction 으로 고른 것이라 dino 에 유리·circular.

라벨 = dominant cell type = argmax(Hist2Cell prediction). 모든 rep per-dim z-score 후 kNN purity(k=10).

## kNN purity (k=10, dominant cell type)

| rep | dim | 4245-BS1 | 4245-TS1 | 4390-BS1 |
|---|---|---|---|---|
| prediction_log1p (ref) | 80 | 0.612 | 0.644 | 0.732 |
| dino768 (원본) | 768 | 0.373 | 0.452 | 0.509 |
| dino19_pca | 19 | 0.361 | 0.434 | 0.487 |
| dino19_imp | 19 | 0.353 | 0.430 | 0.486 |
| hex19 | 19 | 0.360 | 0.424 | 0.477 |
| **dino19pca+hex19** | 38 | **0.381** | **0.449** | **0.498** |
| dino19imp+hex19 | 38 | 0.379 | 0.442 | 0.497 |

![umap](umap_dino19_vs_hex.png)

## 핵심 결론 — 차원 맞추니 **hex 기여가 드러난다**

1. **dino 768→19 는 정보 거의 안 잃음**: dino19(pca) purity 가 dino768 의 96~98% (0.373→0.361 등).
   cell-type 관련 구조는 19-dim 에 대부분 담김 (PCA EVR 0.68).
2. **dino19 ≈ hex19**: 동일 19-dim 에서 둘이 비슷 (hex 가 근소하게 낮음). 한쪽이 압도하지 않음.
3. **dino19+hex19 가 dino19·hex19 둘 다보다 일관되게 높음** (3 슬라이드 전부):
   - vs dino19_pca: **+0.020 / +0.015 / +0.011**
   - vs hex19:     **+0.021 / +0.025 / +0.021**
   → 이전 768+19 concat 에서 **묻혀 있던 hex 의 보완 효과가, 차원을 맞추니 나타남.**
4. **dino19+hex19 ≈ dino768**: 38-dim 결합이 768-dim dino 성능에 근접/동률(BS1 은 +0.008 로 추월).
   즉 hex 19-dim 이 dino 차원 대폭 축소분을 메워준다.
5. PCA ≈ importance (PCA 가 근소 우세). supervised 선택(IMP)이 라벨을 썼는데도 PCA 보다 낫지 않음.

→ 사용자 원래 가설("hex 가 cell-type 정보를 더한다")이 **부분 지지**됨 — 단 효과는 **작다(~+0.02)**.
이전 hex_compare 의 "무의미" 결론은 **차원 지배 아티팩트**였고, 보수 후엔 작지만 **재현성 있는(3/3 슬라이드,
두 축소법 모두) 양의 기여**가 보인다.

> **146 과 종합 (`../hex_dino19_146/summary.md`)**: 224(dino 강함)에선 dino19≈hex19, 결합이 둘 다보다 +0.02.
> 146(OOD-FOV, dino 약함)에선 **hex19 가 dino 를 추월**(dino768 까지), 결합은 dino19 보단 낫지만 hex19 단독엔
> 못 미침. → 차원 맞추면 hex 가 dino 와 대등/우세이고, 결합 이득은 **dino 품질에 의존**.

## 정직한 한계
- 라벨이 **argmax(Hist2Cell prediction)** = H&E-형태 유래 → 지표가 morphology(dino)에 유리.
  따라서 hex 기여는 **과소평가** 가능 (실제 cell-type GT 면 더 클 수도). 결정적 판정엔 Visium GT 필요.
- dino19_imp 는 라벨로 dim 선택(circular)인데도 PCA 대비 이득 없음 → 특정 19 dim 이 특출나진 않음.
- 효과 +0.02 는 moderate 아님, 작음. 방향 일치·재현성으로 읽을 것(수치 절대치 보조).

## 산출물
- `knn_purity.csv`, `knn_purity_pivot.csv`, `umap_dino19_vs_hex.png`, `embeddings/`
- 146 도 동일 가능: `--infer-dir inference_output_146 --dino-dir .../dino_output_146 --agg-dir .../dino_hex_agg_146 --label 146`
