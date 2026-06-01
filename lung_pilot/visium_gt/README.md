# visium_gt — Visium human lung GT 평가 셋 (dino vs dino+hex 결정적 비교용)

목적: TCGA-LUAD 에선 cell-type "정답" 이 Hist2Cell 자기 예측이라 circular 했음
(`../hex_compare_224/summary.md`). 여기선 **실제 cell2location cell-type abundance(GT)** 가 있는
Visium human lung 으로 dino vs dino+hex 를 공정 평가한다.

## 데이터 (repo 내 기존 데이터에서 선택, 2026-06-01)
원본: `example_data/humanlung_cell2location/<sample>.pt`
= `Data(x[N,3,224,224] H&E 224px·112µm, y[N,330]=250유전자+80celltype, edge_index, pos)`.
**cell type GT = y[:,250:] (80-d, cell2location abundance).** cell type 이름은 `cell_types.pkl`(80, AT1 부터).

3장 (서로 다른 배치 + 조직 조성 다양):

| 샘플 | N | GT dominant 상위 | 조직 |
|---|---|---|---|
| WSA_LngSP10193347 | 1937 | AT2 894 / Fibro_alveolar 456 / AT1 352 | 폐포 실질 |
| WSA_LngSP8759313 | 2001 | Ciliated 1654 / Muscle_smooth_syst_arterial 167 | 기도 |
| WSA_LngSP9258468 | 2285 | Fibro_adventitial 973 / SMG_Duct 418 / Ciliated 372 | 혼합 |

## 산출물
- `dino_output/<sample>/features_dinov2.npy` — DINOv2 ViT-B/14 CLS [N,768] (`dino_infer.py` 로 추출, 2026-06-01)
- `gt/<sample>_celltype_gt.npy` — cell-type GT abundance [N,80]
- `gt/cell_types.npy` — 80 cell type 이름 (열 순서)

## 다음 단계 (평가)
1. **HEX feature 를 이 3 샘플에 추출** (동료 HEX 모델) → `[N,19]` 또는 agg `[N,787]`.
   ⚠️ 여기 `x` 는 **20×·224px=112µm FOV**(Hist2Cell 학습). HEX(학습 72.8µm)에 공정히 넣으려면
   원본 Visium WSI+spot 좌표에서 **72.8µm 로 재crop** 필요 — pre-crop `.pt` 만으론 FOV 불일치.
   (DINO 는 generic 이라 그대로 사용.)
2. dino vs dino+hex 를 **실제 GT(y 80 celltype)** 와 비교: dominant-type purity, abundance 회귀(R²),
   neighborhood 정렬 등. TCGA 와 달리 정답이 독립 → 결정적.
