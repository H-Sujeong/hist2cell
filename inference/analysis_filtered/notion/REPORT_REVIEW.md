# Notion 보고서 표현 리뷰

> 사용자가 작성한 Notion 보고서 (TNBC 2개 슬라이드 cell typing + ROI proteomics 연관성) 의 표현 / 일관성 / 정확성 검토. 본 분석 repo (EPITHELIAL_PROXY_METHODOLOGY.md + 4 findings.md + COMPARISON.md) 와의 일치 관점에서 평가.

---

## 🔴 critical (반드시 수정)

### 1. "TNBC 2개 슬라이드" — 확인 필요

목적 첫 줄에 "TNBC 2개 슬라이드" 라고 적혀있는데, 본 분석 문서들 어디에도 slide1/slide2 가 TNBC subtype 으로 확정된 근거가 없다. KBSMC 96-sample cohort 가 breast cancer 이긴 하지만 TNBC-only 인지 confirm 안 됨 — `inference/analysis/README.md` 의 TCGA TNBC 검증은 **EMT/IMMUNE 축의 외부 validation** 용도일 뿐, 본 슬라이드 자체의 분자 분류는 아니다.

→ KBSMC 측에서 두 슬라이드의 ER/PR/HER2 status 확인 후 표시. 안 됐으면 **"KBSMC breast 2개 슬라이드"** 로 중립화. 잘못된 분류로 외부 reader 가 read 하면 신뢰도 타격.

### 2. "tumor region 분석에 이용" — methodology 와 모순

전제 첫 bullet:

> hist2cell open-source trained weight는 lung 조직 밖에 없음, 이 예측 결과를 바탕으로 breast 에 매칭하여 **tumor region 분석에 이용**

`EPITHELIAL_PROXY_METHODOLOGY.md` 의 결론 문장과 정면 충돌. methodology 는 "tumor detector 가 아니다" 를 핵심 메시지로 명시 + 이전 critique 에서도 가장 큰 지적이었던 부분.

**수정 권장**:

> hist2cell open-source trained weight 는 lung 조직 밖에 없음. 이 예측 결과를 **lung-derived epithelial-activity spatial proxy (strict/broad)** 로 변환하여 breast 의 ROI proteomics 와 *정성적* 일치성을 탐구.

### 3. "proteomics 결과 비교 결론" 섹션 — cancer-proxy 용어 잔존

본문 (슬라이드별 findings) 은 epithelial-activity proxy (strict/broad) 로 통일되어 있는데, 마지막 결론 비교 섹션에서만 "cancer-proxy" 가 남아 있어 일관성이 깨진다.

| 위치 | 현재 | 수정 권장 |
|---|---|---|
| 슬라이드 1 활성도 row | "cancer-proxy 우세 spot 10.9%" | "broad-proxy-dominant 10.9% / strict 0.35%" |
| 슬라이드 2 활성도 row | "cancer-proxy 17.7%" | "broad-proxy-dominant 17.7% (원본) / 3.64% (필터) / strict 0.04% — broad-only 신호의 측부 의존이 큼" |
| "유사한 경향" 한계 §2 | "Hist2Cell 의 cancer-proxy 5 type" | "Hist2Cell 의 epithelial-activity proxy (strict 3종 / broad 5종)" |

methodology 와 본문이 새 용어로 통일된 만큼 결론 섹션도 정렬 필요.

### 4. 슬라이드 2 "더 active" 표현 — broad 기준임을 명시

현재:

> 두 modality 모두 slide2 가 slide1 보다 active 하다고 보고

이 표현은 **broad-proxy 기준** 에서만 참 (broad 17.7% vs 10.9%). **strict 기준** 으로는 양쪽 모두 ~0% 로 차이 없음 → "active 의 정의에 따라 결론이 달라진다" 가 정확.

**수정 권장**:

> 두 modality 가 같은 방향성: proteomics 는 환자 2 에서 high/low Tumor 비율이 환자 1 보다 균형적 (13:15 vs 10:21), Hist2Cell 의 broad epithelial-activity proxy 도 환자 2 가 1.6배 우세. 단 strict (방어 가능한 3 종) 으로 보면 양 슬라이드 모두 dominant 영역 ~0% 로 동일 — "active" 의 정의는 broad-proxy 가정에 의존.

### 5. 슬라이드 1 Moran I 수치 — 원본/필터 명시

현재:

> Hist2Cell: AT2/Dividing_AT2 Moran I 0.74-0.75, blob 응집

이건 **원본 (un-filtered) 값**. 필터 후엔 Dividing_AT2=0.679, AT2=0.722. 어느 버전 인용 중인지 모호.

**수정 권장**: "원본 AT2 0.745 / Dividing_AT2 0.749, 필터 후 0.722 / 0.679 — 양쪽 robust" 또는 일관되게 한 버전 골라 표기.

---

## 🟡 권장 (있으면 더 좋음)

### 6. 한계 §1 "정량 검증 불가" + §2 "lung→breast" 의 관계 정리

두 항목이 서로 의존: 좌표 받아도 lung-proxy 한계는 남고, her2st-trained 모델 와도 좌표 없으면 ROI 매칭 안 됨. → 둘이 합쳐져야 진짜 정량 검증.

**추가 권장 문구**:

> ⇒ 정량 검증의 *상한* 은 (a) ROI 좌표 + (b) breast-trained 모델 (her2st CUCA) 두 가지가 모두 갖춰진 후. 본 보고서는 그 전 단계의 정성 cross-check.

### 7. cohort context 누락 — KBSMC 96 sample + TCGA TNBC 검증의 위치 한 줄로

전제에 (B) bulk cohort 만 언급되고, 두 슬라이드가 cohort 의 column 30 / column 3 (중간 ~ 중상위 위치, outlier 아님) 라는 점 + TCGA TNBC 의 EMT/IMMUNE axis 가 본 cohort 에서 재현된다는 점이 빠져 있음 — 본 결과를 "cohort 안에서 일반화 가능한 패턴" 으로 위치시키는 데 도움.

**전제 끝에 한 줄 추가 권장**:

> KBSMC 96-sample 안에서 두 슬라이드는 column 30 (slide1) / column 3 (slide2) 의 평범 ~ 중상위 sample. cohort 의 EMT-high/IMMUNE-low 축은 TCGA TNBC 에서도 재현됨 (외부 검증).

### 8. "사실상 비약적인 명명법" 표현 강화

현재:

> cancer-proxy는 사실 비약적인 명명법 으로 볼 수 있음

→ 좀 더 actionable 하게:

> "cancer-proxy" 라는 이전 명명법은 외부 reader 의 tumor-detection 오해를 유발할 수 있어 본 분석은 **strict / broad epithelial-activity proxy** 두 score 로 분리. methodology 문서 §1 참조.

---

## 🟢 small (선택)

### 9. 첨부 zip 위치

본문 중간 `[slide2_152_19_v2.zip]` 가 "두 슬라이드 비교" 섹션 안에 있어 흐름이 어색. 결론 끝 또는 별도 "첨부" 섹션으로 모으는 게 깔끔.

### 10. slide2 의 strict/broad 갭

"slide2 의 가장 큰 발견 — broad/strict 차이" 박스를 결론 섹션에도 한 줄 정리. 본문에선 잘 드러나는데 결론 비교에선 nuance 가 약함.

---

## 한 줄 요약

본문 (슬라이드별 findings) 은 이미 새 용어로 정리되어 있으니, **❶ 첫머리 "TNBC" 확인 + ❷ "tumor region 분석" 표현 다운 + ❸ 마지막 비교 섹션의 cancer-proxy 표현을 새 용어로 통일 + ❹ "더 active" 같은 비교 claim 의 strict/broad 기준 명시** 이 4 가지가 critical. 나머지는 polish.

이 4 가지 반영 시점에 외부 reviewer 대응도 깔끔해진다.

---

## 관련 파일

- `inference/analysis/EPITHELIAL_PROXY_METHODOLOGY.md` — 본 리뷰의 기준이 되는 methodology
- `inference/analysis/{slide1,slide2}_v2/findings.md` — 본문 (원본 분석)
- `inference/analysis_filtered/{slide1,slide2}_v2/findings.md` — 본문 (필터 분석)
- `inference/analysis_filtered/COMPARISON.md` — 원본 vs 필터 비교
- `inference/analysis_filtered/notion/{slide1,slide2}_v2/findings.md` — Notion 패키지 (사용자가 보고서 작성에 활용)
