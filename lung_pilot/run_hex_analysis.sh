#!/usr/bin/env bash
# 교정 HEX/agg 도착 시 동일 HEX 분석 일괄 재실행.
#
# 사용:
#   ./run_hex_analysis.sh <label:224|146> <infer_dir> <dino_dir> <agg_dir>
# 예 (교정 파일이 fileserver 같은 경로로 다시 오면 그대로):
#   ./run_hex_analysis.sh 224 lung_pilot/inference_output     /mnt/fileserver/lung_pilot/dino_output     /mnt/fileserver/lung_pilot/dino_hex_agg
#   ./run_hex_analysis.sh 146 lung_pilot/inference_output_146 /mnt/fileserver/lung_pilot/dino_output_146 /mnt/fileserver/lung_pilot/dino_hex_agg_146
#
# 가정: agg = [dino768 ⊕ hex19] = 787-d (hex = agg[:,768:]). 포맷 다르면 스크립트 수정 필요.
# 산출(lung_pilot/ 루트, 덮어씀):
#   hex_compare_<label>/ hex_dino19_<label>/ hex_dino19_denoise_<label>/ hex_dino19_top50_<label>/
set -euo pipefail
LABEL="${1:?label(224|146)}"; INFER="${2:?infer_dir}"; DINO="${3:?dino_dir}"; AGG="${4:?agg_dir}"
ROOT=/home/sjhong/hist2cell
PY="$ROOT/.venv/bin/python"
cd "$ROOT"

echo "### [$LABEL] hex_compare (3×3 UMAP + purity + 가중분석)"
"$PY" lung_pilot/hex_compare.py --infer-dir "$INFER" --dino-dir "$DINO" --agg-dir "$AGG" \
  --out-dir "lung_pilot/hex_compare_$LABEL" --label "$LABEL"

echo "### [$LABEL] hex_dino19 (dino→19 PCA/중요도 + Q1 chance/excess)"
"$PY" lung_pilot/hex_dino19_compare.py --infer-dir "$INFER" --dino-dir "$DINO" --agg-dir "$AGG" \
  --out-dir "lung_pilot/hex_dino19_$LABEL" --label "$LABEL"

echo "### [$LABEL] denoise (hex 하위10%→0)"
"$PY" lung_pilot/hex_dino19_denoise.py --infer-dir "$INFER" --dino-dir "$DINO" --agg-dir "$AGG" \
  --out-dir "lung_pilot/hex_dino19_denoise_$LABEL" --label "$LABEL" --pct 10

echo "### [$LABEL] top50 (hex 상위50%만)"
"$PY" lung_pilot/hex_dino19_denoise.py --infer-dir "$INFER" --dino-dir "$DINO" --agg-dir "$AGG" \
  --out-dir "lung_pilot/hex_dino19_top50_$LABEL" --label "$LABEL" --pct 50

echo "### [$LABEL] DONE → lung_pilot/hex_*_$LABEL/"
