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
| ③ Inference | ✅ 완료 (2026-05-24) | `inference_output/<slide>/predictions.{csv,npy}` |

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
└── inference_output/    # Hist2Cell 추론 결과 — <slide>/predictions.{csv,npy} + _logs/
```
세부 문서: 각 `tilitng_output/*/tiling_summary.md`, `graph_output/README.md`.

## 타일 / 노드 수

| 슬라이드 | 224 | 112 |
|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 | 11,470 |
| TCGA-05-4245-01A-01-TS1 | 1,871 | 7,499 |
| TCGA-05-4390-01A-01-BS1 | 10,661 | 42,615 |

## 추론 결과 (2026-05-24 완료)

| slide | spots | cell types | row_sum (abundance) | shard 시간 (4× A5000) |
|---|---|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 2,869 | 80 | 1.32 – 62.81 | 23.6 s |
| TCGA-05-4245-01A-01-TS1 | 1,871 | 80 | 1.85 – 49.43 | 7.2 s |
| TCGA-05-4390-01A-01-BS1 | 10,661 | 80 | 1.70 – 49.28 | 60.0 s |

NaN 없음. CSV 컬럼 = `spot_id, X, Y, <80 cell types>` (총 83), `.npy` 는 (N, 80).
row_sum 은 cell2location abundance scale 그대로 (확률 아님).

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

## 다음 단계

1. **HEX / DINO** (repo 외부 모델, 사용자 측) — `graph_output/112/*.pt` 사용.
   HEX expression + DINO 벡터 concat → UMAP.
2. **비교 UMAP** — Hist2Cell cell-type UMAP (위 `predictions.npy`) vs `[HEX⊕DINO]` UMAP.

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
