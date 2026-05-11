# Notion 보고서 — "proteomics 결과 비교 결론" 섹션 수정안

> 본 섹션 전체를 아래로 교체. 본문 (슬라이드별 findings) 의 새 용어 (broad / strict epithelial-activity proxy) 와 정렬, slide1/slide2 비교의 strict vs broad 갭 명시.

---

## proteomics 결과 비교 결론

### 슬라이드 1 (1_085_12) — 깔끔한 일치 ✓

지표: 활성도
Proteomics: low-risk Tumor (b:21) ≫ high-risk (a:10) → 대체로 quiescent
Hist2Cell: broad epithelial-activity proxy-dominant 10.87% (원본) / 13.0% (필터), strict 0.35% / 0.59% — 양쪽 모두 소수
────────────────────────────────────────
지표: proliferative 신호
Proteomics: high-risk 마커 KIF20A/KIF22/INCENP (mitosis)
Hist2Cell: Dividing_AT2 / AT2 Moran I 원본 0.749 / 0.745, 필터 후 0.679 / 0.722 — 양쪽 모두 robust blob
────────────────────────────────────────
지표: stromal context
Proteomics: high-risk 마커에 MYH11/TAGLN (smooth muscle)
Hist2Cell: Stromal-muscle μ=2.23 (원본 1위), Fibro 합 #2 — 필터 후 μ=3.58 로 상승하지만 순위 유지

→ "stromal-rich + 일부 영역에 mitotic epithelial" 의 같은 그림. **strict / broad 양쪽 + 원본 / 필터 모든 조합에서 결론 보존**.

---

### 슬라이드 2 (2_152_19) — 일치하되 nuance 큼 ⚠️

지표: 활성도
Proteomics: high/low Tumor 균형 (e:13, f:15) — 환자 1 보다 균형적
Hist2Cell: **broad-proxy-dominant 17.7% (원본) → 3.64% (필터) — 5 배 감소**, strict 0.04% / 0.03% (양쪽 ~0). → broad 17.7% 결론이 측부 덩어리 (전체 30%) 의 AT2 / Suprabasal 신호에 강하게 의존, strict 으론 본 슬라이드 전체에서 dominant 영역 사실상 없음.
────────────────────────────────────────
지표: high-risk Tumor 마커
Proteomics: GZMH/LCK (immune mixed!) + MAPK12/MARK3 (cytoskeleton)
Hist2Cell (필터): 큰 덩어리 안 immune dominant 96.4%, top Moran R 페어가 DC/Macro myeloid 중심 (B_memory 도 top5)
────────────────────────────────────────
지표: high-risk T-cell 마커
Proteomics: TFAP2C (mammary epithelial transcription factor!)
Hist2Cell (필터): SMG_Duct ↔ SMG_Serous spatial co-occurrence Moran R=0.604 — ductal-glandular 신호 후보

→ 두 modality 가 동일한 "active + immune-mixed + ductal-epithelial 동반" 그림. 단 원본의 "Goblet ↔ immune mutual exclusion" 결론은 측부 덩어리 의존 (필터 후 사라짐) — Hist2Cell 단독 신호이며 proteomics 가 직접 검증한 것이 아님 (mucin marker 미측정). **broad-proxy 17.7% 표현은 외부 reader 에게 단순화시켜 전달 금지** — strict 으로 재검증 시 사실상 0% 이며, 결론의 robust 부분은 AT2 hot-spot blob + B-cell TLS + ductal-glandular signal 정도.

---

### 두 슬라이드 비교 — modality-cross 일관성

- proteomics 는 환자 2 에서 high/low Tumor 비율이 환자 1 보다 균형적 (13:15 vs 10:21), Hist2Cell 의 broad epithelial-activity proxy 도 환자 2 가 1.6배 우세 (broad 17.7% vs 10.9%, 원본 기준). 단 **strict (방어 가능한 3 종) 으로 보면 양 슬라이드 모두 dominant 영역 ~0% 로 동일** — "active" 의 정의는 broad-proxy 의 AT2/Suprabasal 가정에 의존.
- 두 modality 모두 immune compartment 의 spatial 존재 detect (slide1 은 NK + Monocyte / B-cell 응집형, slide2 는 myeloid 분산 + B-cell 응집부 측부).
- proliferative-like epithelial signal: slide1 = AT2 / Dividing_AT2 Moran I 0.7+ blob + proteomics KIF20A/22/INCENP mitosis 일치, slide2 = SMG_Duct / SMG_Serous co-occurrence + TFAP2C 일치 (단 ductal-glandular 신호의 cross-tissue mapping 은 가설 수준).

---

### ⚠️ "유사한 경향" 의 한계

1. **정량 검증 불가** — `.tmpprotocol` ROI 좌표 없어서 ROI-spot 매칭 정량 (Wilcoxon, ROI-level Pearson) 미수행. 정성 cross-check 까지만.
2. **lung → breast cross-tissue proxy 의 한계** — Hist2Cell 의 5 lung label 합 (strict 3종 / broad 5종) 은 *epithelial-activity 의 spatial proxy* 이지 breast cell type 의 직접 측정이 아님. 이전 명명법 "cancer-proxy" 는 외부 reader 의 tumor-detection 오해를 유발할 수 있어 본 분석은 **strict / broad** 두 score 로 분리하여 사용. broad-only 라벨 (AT2 / Suprabasal) 의 cross-tissue 매핑은 가설 수준 — slide2 의 "broad-proxy 17.7%" 가 측부 의존인 점에서 직접 확인됨. methodology 문서 (`EPITHELIAL_PROXY_METHODOLOGY.md`) §1, §3 참조.

⇒ 정량 검증의 *상한* 은 (a) ROI 좌표 + (b) breast-trained 모델 (her2st CUCA) 두 가지가 모두 갖춰진 후 가능. 본 보고서는 그 전 단계의 정성 cross-check.

⇒ spatial proteomics 에 사용된 coords information 제공받은 후, 그리고 CUCA her2st 가중치 도착 후, 정확한 ROI-level 정량 분석 추가 예정.
