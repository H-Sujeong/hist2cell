# HEX 분석 재실행 가이드 (교정 HEX 도착 시)

## 배경
2026-06-01 기준 동료 HEX/optimus 추출에 **치명적 결함**이 발견되어, 기존 HEX 기반 분석은
`_archive_invalid_hex/` 로 아카이브(무효). 교정된 HEX(agg) 파일이 오면 **동일 분석을 그대로** 재실행한다.
**224·146 둘 다** 올 예정.

## 입력 가정 (중요)
- `agg = [dino 768 ⊕ hex 19] = 787-d` 형식 (hex = `agg[:,768:]`). 동료가 같은 포맷으로 주면 수정 불필요.
- 포맷이 바뀌면(hex 별도 파일, 차원 변경 등) `hex_compare.py` / `hex_dino19_compare.py` /
  `hex_dino19_denoise.py` 의 `agg[:,768:]` 슬라이스만 고치면 됨.
- spot 순서는 우리 `inference_output{,_146}` (= graph_output{224,146} 순서)와 정합되어야 함
  (146 은 동료 146 타일링 공유라 정합 OK. 224 는 우리 타일링이라 동료가 그 좌표로 떠야 함).

## 재실행 (한 줄씩)
```bash
cd /home/sjhong/hist2cell
# 교정 파일이 같은 fileserver 경로로 다시 오면 그대로:
lung_pilot/run_hex_analysis.sh 224 lung_pilot/inference_output     /mnt/fileserver/lung_pilot/dino_output     /mnt/fileserver/lung_pilot/dino_hex_agg
lung_pilot/run_hex_analysis.sh 146 lung_pilot/inference_output_146 /mnt/fileserver/lung_pilot/dino_output_146 /mnt/fileserver/lung_pilot/dino_hex_agg_146
```
(교정본이 다른 경로면 3·4번째 인자만 교체.)

## 산출 (lung_pilot/ 루트에 재생성)
| 폴더 | 내용 |
|---|---|
| `hex_compare_<L>/` | 3×3 UMAP(pred/dino/hex+dino) + knn_purity{,_weighting}.csv + summary |
| `hex_dino19_<L>/` | dino→19(PCA/ANOVA-F 중요도) vs hex19 vs dino19+hex19, **Q1 chance/excess 포함** |
| `hex_dino19_denoise_<L>/` | hex 하위10%→0 |
| `hex_dino19_top50_<L>/` | hex 상위50%만 |

## 분석 정의 (요지 — 자세히는 각 summary)
- **Q1** = representation *자체* 가 cell type 으로 뭉치나 = rep 단독 purity 의 **chance 대비 excess**.
- **Q2** = hex 가 dino 에 *추가* 보탬 = `purity(dino+hex) − purity(dino)`.
- **chance = size-weighted ∑ᵢpᵢ²** (이웃을 라벨분포대로 무작위 추출 시 같은 cell type 일 기대확률).
- 핵심: dino768 차원이 hex19 를 압도 → **dino 를 19-d 로 축소**(hex_dino19)해야 공정.

## 재실행 후 할 일
1. 새 결과의 Q1(excess)·Q2 가 무효본과 어떻게 달라졌는지 비교.
2. summary 작성(한국어 + 이미지 임베드 + Q1/Q2/chance 정의 명시).
3. (선택) flow_overview 재생성 — `flow_overview.py` 의 폴더 경로를 새 결과로 맞춘 뒤 실행.
4. **Visium GT(`visium_gt/`)로 확장** — 교정 HEX 를 Visium 3장에 뽑으면 실제 cell-type GT 로 결정 평가
   (단 그 전에 optimus/HEX 추출이 제대로 되는지부터 검증할 것: H-optimus-0 HF 가중치 + HF 전처리).

## 무효본 위치
`_archive_invalid_hex/` (결론 인용 금지, 참고용).
