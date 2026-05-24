# TCGA-LUAD WSI 타일링 결과 (112×112)

`WSI_tile_sampling_framework` 의 `WSITileSampler` 로 TCGA-LUAD 슬라이드 3장을
**112×112 px** 타일로 분할한 결과입니다. (224×224 버전은 `../TCGA-LUAD/` 참조)

## 왜 112 인가 — HEX 모델(40× 학습) 경로용

| | tile 한 변 | 20× 슬라이드에서 본 물리 크기 |
|---|---|---|
| 224 타일 (`../TCGA-LUAD/`) | 224 px | 112 µm — Hist2Cell(20× 학습)용 |
| **112 타일 (본 폴더)** | **112 px** | **56 µm** — HEX 모델(40×·224 학습)용 |

TCGA-LUAD 3장은 native **20× (mpp ≈ 0.502)**. 40×·224 로 학습된 HEX 모델은
56 µm(= 40× 의 224 px) 시야를 기대하므로, 20× 슬라이드에서는 **112 px** crop 이
그 물리 크기에 해당합니다.

> **HEX 추론 시**: 112×112 패치를 그대로 넣지 말고 **224×224 로 ×2 resize** 해서
> 입력해야 합니다 (HEX 입력층도 224 고정). 본 타일링은 112 px 좌표만 제공합니다.

## 실행 방법

```bash
cd WSI_tile_sampling_framework
/home/sjhong/hist2cell/.venv/bin/python run_tiling_tcga_luad.py \
    --tile-size 112 \
    --output /home/sjhong/hist2cell/lung_pilot/tilitng_output/112
```

파라미터: `tile_size=112`, `overlap=0`(stride 112), `min_tiles=5`,
`is_normalized=False`. centroid 버그 수정본(`Tiling.py`) 적용 상태.

## 슬라이드별 결과

| 슬라이드 | 원본 크기 (px) | 썸네일 | scaler | 조직 비율 | 타일 수 (112) | 참고: 224 타일 수 |
|---|---|---|---|---|---|---|
| TCGA-05-4245-01A-01-BS1 | 18000 × 19921 | 2250 × 2490 | 8 | 40.2% | 11,470 | 2,869 |
| TCGA-05-4245-01A-01-TS1 | 14000 × 13527 | 1750 × 1690 | 8 | 49.7% | 7,499 | 1,871 |
| TCGA-05-4390-01A-01-BS1 | 36001 × 21855 | 2250 × 1365 | 16 | 68.1% | 42,615 | 10,661 |

stride 를 양축 모두 절반으로 줄였으므로 타일 수는 224 버전의 약 4배입니다.

## 산출물 구조

```
lung_pilot/tilitng_output/112/
├── <name>.h5                          # full-res 타일 좌표 + metadata (tile_size=112)
├── Thumbnails/<basename>.png          # 슬라이드 썸네일 (실제 배열 해상도)
├── Masks/<name>_tissue_mask.png       # 조직(foreground) mask, 0/255
├── Overlays/<name>_tiles.png          # 썸네일 위 타일 위치 overlay (QA)
└── tiling_summary.md                  # 본 문서
```

`coords` = 원본 해상도(level 0) 기준 타일 좌상단 `(x, y)`, 모두 112 배수.
패치 추출: `openslide.OpenSlide(svs).read_region((x, y), 0, (112, 112))` → 224 로 resize.

## 주의

- **그리드 비정렬**: 112 타일 중심(`coord+56`)은 224 타일 중심(`coord+112`)과
  일치하지 않습니다. HEX(112)와 Hist2Cell(224)을 **spot 단위로 paired 비교**하려면
  같은 center 를 공유해야 하므로, 별도 정렬(공유 center 에서 두 크기 패치 추출)이
  필요합니다. 본 폴더는 독립 112 그리드입니다.
- **해상도 OOD**: TCGA-LUAD 는 native 20× — 진짜 40× 디테일이 없습니다.
  112→224 업샘플 패치는 면적(56µm)은 HEX 학습과 맞지만 해상도는 더 흐립니다.
  HEX 출력 해석·비교 시 이 점을 명시할 것.
