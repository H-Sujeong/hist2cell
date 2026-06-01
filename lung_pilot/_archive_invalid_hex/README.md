# _archive_invalid_hex — 무효 처리된 HEX 기반 분석 (2026-06-01)

## 왜 무효인가
동료의 **HEX feature 추출에 치명적 결함**이 있었음(HEX 잘못 뽑힘). 이 폴더의 모든 분석은
`/mnt/fileserver/lung_pilot/dino_hex_agg{,_146}` 의 **잘못된 hex19 블록**(agg[:,768:])에
의존하므로 결과·결론을 **신뢰할 수 없음**.

(참고: 우리 H-optimus-0 재현 검증에서도 동료 optimus 와 cosine ~0.56 으로 안 맞았음 →
optimus/HEX 파이프라인 자체 문제. `_optimus_validate.log`.)

## 보관 내용 (참고용, 결론 인용 금지)
- hex_compare_{224,146}/        : 3×3 UMAP + kNN purity (per-dim/block-EQ) + overlap
- hex_dino19_{224,146}/         : dino→19dim 축소 비교
- hex_dino19_denoise_{224,146}/ : hex 하위10%→0
- hex_dino19_top50_{224,146}/   : hex 상위50%만
- *.log                         : 실행 로그

## 교정 HEX 도착 시
`lung_pilot/run_hex_analysis.sh` 로 동일 분석을 새 입력에 재실행 → 새 결과는 lung_pilot/ 루트에 생성.
스크립트(hex_compare.py / hex_dino19_compare.py / hex_dino19_denoise.py)는 그대로 재사용.
자세한 절차: `lung_pilot/HEX_RERUN.md`.
