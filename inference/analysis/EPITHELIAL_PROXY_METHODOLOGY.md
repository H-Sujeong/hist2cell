# Epithelial-activity proxy (strict / broad) — 5 lung label 의 선정 근거 (lung-trained Hist2Cell → breast 적용)

> **이 문서가 필요한 이유**
>
> 우리가 KBSMC breast 슬라이드 (slide1, slide2) 에 적용한 Hist2Cell 가중치 `humanlung_cell2location_leave_A50_out.pth` 는 **healthy human lung** 데이터로 학습된 모델이다. 80개 출력 cell type 은 모두 lung 분류 라벨이며 breast 에 직접 대응하는 type 이 없다. 그럼에도 우리는 `cell_type_groups.csv` 에서 5 개 lung label (`AT2`, `Basal`, `Suprabasal`, `Dividing_AT2`, `Dividing_Basal`) 을 **합산 score 로** 사용해왔다. 본 문서는 이 선정을 정직한 한계와 함께 명시한다.
>
> ⚠️ **결정적 caveat (먼저 읽기)**
>
> 이 5 lung label 합은 **breast tumor detector 가 아니다.** lung-trained 모델이 breast H&E 에서 *상대적으로 epithelial / progenitor / proliferative morphology 와 유사하다고 판단한 spatial signal* 을 모은 것이다. 따라서 "cancer-proxy" 라는 표현 (이전 버전 문서에서 사용) 은 외부 reader 의 오해를 유발할 수 있어 본 문서부터는 **`epithelial-activity proxy`** 라는 neutral 용어로 통일한다. 정량적 tumor abundance 측정이 아니라 **상대적 spatial pattern 비교** 에 한정해 사용하며, 최종 타당성은 breast-trained 모델 (CUCA her2st 가중치) 또는 pathologist ROI / IHC 와의 spatial overlap 으로 검증되어야 한다.

---

## 1. 두 단계 score 설계 (strict / broad)

비판적 검토 결과, 5 lung label 의 cross-tissue 매핑 신뢰도가 균일하지 않다는 점이 명백해졌다. 따라서 본 분석은 **두 개의 score 를 병렬로 산출** 한다.

| score | 포함 라벨 | type 수 | 해석 |
|---|---|---:|---|
| **strict** epithelial-proliferative proxy | `Basal`, `Dividing_AT2`, `Dividing_Basal` | 3 | 가장 방어 가능. basal stem cell marker (KRT5/TP63) 공유 + 직접 cell-cycle 표현형 |
| **broad** epithelial-activity proxy | 위 + `AT2`, `Suprabasal` | 5 | 이전 버전과 호환. AT2 는 cross-tissue 유사 가설, Suprabasal 은 sensitivity 항목 |

`cell_type_groups.csv` 의 컬럼:
```
is_strict_proxy=1  → strict 셋 (3종) 에 포함
is_broad_proxy=1   → broad 셋 (5종) 에 포함 (strict 포함)
```

분석 산출 (`analyze.py`) 은 두 score 모두 별도 row 로 보고 — 두 score 의 일치는 결론의 robustness 를, 차이는 Suprabasal/AT2 에 의존하는 부분을 드러낸다.

---

## 2. 개별 라벨 별 근거 (tone-down)

### 2.1 Basal (airway basal cell) — **가장 강한 cross-tissue 근거**
- 기도 (airway) 상피의 **basal layer stem cell** (TP63+/KRT5+/KRT14+). 자기 갱신 + ciliated/secretory 로 분화.
- **Rock et al. 2009** [4]: 마우스 trachea + 인간 기도 상피에서 basal cell 이 multipotent stem cell 임을 lineage tracing 으로 보고. **단 본 연구는 postnatal growth + repair 맥락**이며, lung cancer cell-of-origin 의 단정적 증거가 아니다.
- **Sutherland & Berns 2010** [5]: lung cancer origin 의 review — Basal 이 LUSC 후보로 거론되지만 context-dependent.

**breast 맥락**: Lung airway basal 과 **mammary basal/myoepithelial** 은 동일 marker (KRT5+/TP63+/KRT14+) 를 공유 → **marker-level analogy 가 가장 강한 lung label**. 다만 **basal-like / TNBC 의 cell-of-origin 을 basal cell 로 단정하지 않는다** — BRCA1-associated basal tumor 의 candidate target population 으로 **luminal progenitor** 가 보고되어 있어 (Lim et al. 2009 [7]) "basal hot-spot = TNBC origin" 의 단순 매핑은 반박 가능. **본 분석에서는 spatial marker analogy 의 사후 검증 대상으로 둔다.**

→ 분류: strict + broad 양쪽.

### 2.2 Dividing_Basal — **직접 cell-cycle 표현형, 단 H&E 기반 모델 한계 명시**
- Basal 의 부분집합 중 명시적으로 mitotic (Ki67+/cyclin-positive) 인 cell. Madissoon atlas [3] 에서 cell-cycle gene 발현으로 sub-cluster 분리.
- 정상에서도 항상 일정 비율 존재. 활성 (regenerative / 종양 미세환경) 환경에서 abundance 증가.

⚠️ **중요한 caveat**: Hist2Cell 은 H&E 형태학 기반 모델이다. **단일 spot 의 Dividing_Basal abundance 예측이 실제 Ki67/MKI67 단백질 / mRNA 발현을 직접 검출했다는 뜻이 아니다.** 모델이 H&E 형태에서 cell atlas 의 cell-cycle cluster 와 유사한 morphology 를 식별했음을 의미. 따라서 *model-derived proliferative-like signal* 로 해석.

→ 분류: strict + broad 양쪽.

### 2.3 Dividing_AT2 — **strict 의 alveolar branch**
- AT2 의 부분집합 중 mitotic cell. lung distal alveolar 의 활성 분열 세포.
- **Barkauskas et al. 2013** [1]: AT2 가 polyclonal 으로 분열하며 alveolar epithelium 재생을 담당함을 보고 — **AT2 stem cell 기능의 직접 lineage tracing 증거**. 단 lung 맥락.
- **Kim et al. 2005** [2]: AT2 (BASC 포함) 가 K-Ras 변이 시 LUAD 의 cell-of-origin 후보. 단 마우스 실험.

**breast 맥락**: breast 에 직접 대응하는 alveolar progenitor 라벨은 없으나, model-derived proliferative-like signal 로서 strict set 에 포함. *mammary luminal progenitor 와의 직접 매핑은 별도 검증 대상이며 본 score 가 이를 의미하지 않는다*.

→ 분류: strict + broad 양쪽.

### 2.4 AT2 (alveolar type 2 cell) — **broad 의 alveolar branch, cross-tissue 매핑은 가설 수준**
- 폐포의 **distal lung stem cell**, SP-C / SP-B / ABCA3 발현. 자기 갱신 + AT1 분화.
- Barkauskas 2013 [1] / Kim 2005 [2] 의 lung 맥락 stem 활성은 견고.

⚠️ **breast 맥락의 어려움**: AT2 는 lung-specific 라벨이며 breast 에는 직접 대응하는 cell type 이 없다. 이전 버전 문서에서 "mammary luminal progenitor 의 신호로 read 가능" 으로 적었으나 이는 **검증 가설** 수준의 진술이며 본 버전에서는 톤다운한다:

> AT2 label 은 breast cell type 에 직접 대응하지 않는다. epithelial progenitor / secretory epithelial morphology 에 대한 cross-tissue 유사 signal 의 *후보* 로만 해석하며, breast 맥락에서는 luminal/progenitor-like epithelial compartment 와의 spatial overlap 을 검증 대상으로 삼는다.

mammary luminal progenitor 연구 [7][9] 는 본 score 가 *대응한다* 는 근거가 아니라 *검증 가설의 background* 로 인용된다.

→ 분류: broad 만. strict 에는 미포함 (lung-specific 한 lineage 라 cross-tissue 신뢰도 중간).

### 2.5 Suprabasal — **가장 약한 근거, sensitivity 항목**
- Basal layer 위, 분화 transitional state. KRT5 감소 + KRT8 증가. cell cycle 음성이지만 분화-중간 표현형.
- 정상에서는 short-lived intermediate. 염증/수복 환경에서 abundance 증가.

⚠️ **본 분석에서 가장 보수적으로 다룬다**:
- breast 에 1:1 대응 label 없음.
- 이전 버전 문서의 "DCIS → invasive carcinoma transitional epithelial state 와 functional analog" 는 **직접 문헌 근거가 없는 과해석** 으로 판단, 본 버전에서 삭제.
- 본 분석에서는 broad set 의 **auxiliary / sensitivity 항목** 으로만 사용. 단독 해석하지 않으며, strict vs broad 결과 비교 시 Suprabasal 의 기여를 분리해 평가.

→ 분류: broad 만. strict 미포함.

---

## 3. lung → breast 매핑 신뢰도 정리

| lung label | breast 매핑 신뢰도 | 본 분석에서의 역할 |
|---|---|---|
| `Basal` | **높음** (KRT5/TP63/KRT14 marker 공유) | strict + broad |
| `Dividing_Basal` | **높음** (basal stem + cell-cycle marker, H&E 한계 명시 후) | strict + broad |
| `Dividing_AT2` | 중 (lung-specific lineage, cell-cycle 신호 직접성 일부) | strict + broad |
| `AT2` | 중-낮 (cross-tissue 유사성 가설 수준) | broad only |
| `Suprabasal` | 낮 (direct breast counterpart 없음, transitional state) | broad only (sensitivity) |

따라서 **strict score 의 해석 신뢰도 ≫ broad score**. 두 score 의 일치 → 결론 robust. 차이 → AT2/Suprabasal 에 의존하는 부분.

> 본 5 label 합은 her2st 의 mammary epithelial 라벨 (`basal cell`, `luminal epithelial cell of mammary gland`, `mammary gland epithelial cell`) 과 **비교 가능한 후보 epithelial-activity score** 로 볼 수 있다. 그러나 1:1 대응은 아니므로, CUCA / her2st 가중치 도착 후 spatial overlap 검증을 통해 사후적으로 타당성을 평가한다.

---

## 4. 명시적 한계 (외부 reader 용)

1. **proxy ≠ detector** — 이 score 들은 종양 detection 이 아니라 *epithelial activity* 의 spatial reference. 진짜 종양 영역 검출은 별도 모델 (tiatoolbox 위험도, CUCA her2st, IHC) 가 담당.
2. **모델 feature space mismatch** — lung 조직으로 학습한 backbone (ResNet18) 의 representation 이 breast H&E 와 어떤 차이를 보이는지는 정량 검증 없음. *상대 spatial pattern* 으로만 해석.
3. **H&E 기반 모델 한계** — `Dividing_*` 라벨이라 해도 실제 Ki67/MKI67 abundance 의 직접 측정이 아니라 **morphology-similar signal**. 정량 검증에는 IHC / multiplex IF 필요.
4. **mpp / tile_size mismatch** — Visium 20× (~0.5 μm/px) 학습 vs Aperio 40× (0.2615 μm/px) 적용. 모델 시야가 학습 분포의 절반 — 절대값 비교 금지.
5. **분리된 tissue blob 의 영향** — `inference/analysis_filtered/COMPARISON.md` 에서 보았듯, slide2 의 경우 epithelial-activity 우세 spot 중 상당 부분이 가장 큰 덩어리 바깥에 위치. abundance 의 절대값 / spatial 분포 해석 시 슬라이드 내 다중 compartment 가능성 명심.
6. **5 외 다른 후보 type 제외 사유**: lung 80 type 중 `Macro_dividing` (면역세포), `B_plasmablast` (면역세포), `Deuterosomal` (ciliated-progenitor intermediate) 는 epithelial cancer compartment 의 spatial proxy 로 부적합하여 제외. Suprabasal 는 transitional epithelial 이라 broad 에만 포함.

---

## 5. CUCA / her2st 가중치 도착 후 권장 흐름

| 단계 | 본 lung-기반 분석 | her2st 가중치 도착 후 |
|---|---|---|
| epithelial-activity 정의 | strict 3종 + broad 5종 | her2st 의 3 mammary epithelial type (`basal cell`, `luminal epithelial cell of mammary gland`, `mammary gland epithelial cell`) 합 |
| `is_*` flag | `is_strict_proxy` / `is_broad_proxy` | `is_mammary_epithelial=1` (3 종) → analyze.py 1줄 패치 |
| 해석 | cross-tissue 유사성 기반 spatial proxy | 직접 라벨 측정 |
| 신뢰도 | 정성 cross-check 까지 | 정량 ROI-level 분석 가능 (ROI 좌표 받은 후) |
| 검증 | 본 score 의 spatial overlap 가설 | strict / broad / mammary 세 score 모두 같은 hot-spot 잡는지 일치성 검증 |

→ 두 결과의 spatial overlap (Pearson ρ / Jaccard) 측정 시, 본 strict / broad score 의 사후 타당성이 직접 검증된다.

---

## 6. 참고문헌

[1] **Barkauskas CE**, Cronce MJ, Rackley CR, et al. *Type 2 alveolar cells are stem cells in adult lung*. **J Clin Invest** 2013;123(7):3025–3036. doi:10.1172/JCI68782 (PMID: 23921127)

[2] **Kim CFB**, Jackson EL, Woolfenden AE, et al. *Identification of bronchioalveolar stem cells in normal lung and lung cancer*. **Cell** 2005;121(6):823–835. doi:10.1016/j.cell.2005.03.032

[3] **Madissoon E**, Oliver AJ, Kleshchevnikov V, et al. *A spatially resolved atlas of the human lung characterizes a gland-associated immune niche*. **Nat Genet** 2023;55:66–77. doi:10.1038/s41588-022-01243-4 — 본 모델 학습 데이터의 source atlas, 80 lung cell types via Cell2location

[4] **Rock JR**, Onaitis MW, Rawlins EL, et al. *Basal cells as stem cells of the mouse trachea and human airway epithelium*. **PNAS** 2009;106(31):12771–12775. doi:10.1073/pnas.0906850106 — postnatal growth + repair 맥락의 airway basal stem cell 증거

[5] **Sutherland KD**, Berns A. *Cell of origin of lung cancer*. **Mol Oncol** 2010;4(5):397–403. doi:10.1016/j.molonc.2010.05.002 — context-dependent lung cancer origin review

[6] **Wuidart A**, Sifrim A, Fioramonti M, et al. *Early lineage segregation of multipotent embryonic mammary gland progenitors*. **Nat Cell Biol** 2018;20:666–676. doi:10.1038/s41556-018-0095-2 — **embryonic mammary gland progenitor 의 lineage segregation** (p63/basal fate 측면). adult mammary basal stem cell 의 직접 증거가 아니라 lineage commitment 의 background.

[7] **Lim E**, Vaillant F, Wu D, et al. *Aberrant luminal progenitors as the candidate target population for basal tumor development in BRCA1 mutation carriers*. **Nat Med** 2009;15(8):907–913. doi:10.1038/nm.2000 — **basal-like tumor 의 cell-of-origin 으로 luminal progenitor 가 후보** (basal 단정 완화의 근거)

[8] **Visvader JE**. *Keeping abreast of the mammary epithelial hierarchy and breast tumorigenesis*. **Genes Dev** 2009;23(22):2563–2577. doi:10.1101/gad.1849509 — mammary cell hierarchy review

[9] **Pal B**, Chen Y, Vaillant F, et al. *A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast*. **EMBO J** 2021;40(11):e107333. doi:10.15252/embj.2020107333 — breast cell atlas

[10] **Sørlie T**, Perou CM, Tibshirani R, et al. *Gene expression patterns of breast carcinomas distinguish tumor subclasses with clinical implications*. **PNAS** 2001;98(19):10869–10874. doi:10.1073/pnas.191367098 — luminal A/B / HER2 / basal-like 분자 분류 background

[11] **Andersson A**, Larsson L, Stenbeck L, et al. *Spatial deconvolution of HER2-positive breast cancer delineates tumor-associated cell type interactions*. **Nat Commun** 2021;12:6012. doi:10.1038/s41467-021-26271-2 — her2st dataset 의 원본 (CUCA 학습 데이터)

[12] **Travaglini KJ**, Nabhan AN, Penland L, et al. *A molecular cell atlas of the human lung from single-cell RNA sequencing*. **Nature** 2020;587:619–625. doi:10.1038/s41586-020-2922-4

[13] **Hist2Cell** (used model). Zhao W. et al. *Hist2Cell* GitHub: https://github.com/Weiqin-Zhao/Hist2Cell — paired histology + spatial transcriptomics → 80 lung cell type prediction. 본 분석 가중치 `model_weights/humanlung_cell2location_leave_A50_out.pth` 의 source.

---

## 7. 결론 (외부 reader 안전 버전)

> 본 분석에서 사용한 5 개 lung-derived Hist2Cell label (`AT2`, `Basal`, `Suprabasal`, `Dividing_AT2`, `Dividing_Basal`) 은 breast cancer cell type 을 직접 예측하는 label set 이 아니다. 이들은 lung atlas 내 epithelial stem / progenitor, transitional, 또는 cell-cycle-associated state 로 구성되어 있으며, breast H&E 적용 시에는 tumor detector 가 아니라 **epithelial / proliferative-like spatial activity 를 나타내는 exploratory proxy** 로 해석한다.
>
> 본 분석은 두 단계 score 를 병렬로 산출한다:
> - **strict** (`Basal`, `Dividing_AT2`, `Dividing_Basal`): basal stem cell marker 공유 + 직접 cell-cycle 표현형 → 가장 방어 가능
> - **broad** (위 + `AT2`, `Suprabasal`): 이전 버전과 호환, AT2/Suprabasal 의 cross-tissue 매핑은 검증 가설
>
> 따라서 정량적 tumor abundance 가 아닌 **상대적 spatial pattern 비교** 에 한정해 사용하며, 최종 타당성은 breast-trained her2st / CUCA output 또는 pathologist ROI / IHC 와의 spatial overlap 으로 검증해야 한다. 두 score (strict, broad) 의 일치는 결론의 robustness 의 척도이며, 차이는 AT2 / Suprabasal 에 의존하는 부분을 드러낸다.

---

## 8. 본 문서와 연결된 파일

- `inference/analysis/cell_type_groups.csv` — `is_strict_proxy` / `is_broad_proxy` 컬럼
- `inference/analysis/analyze.py` — `per_group_stats()` 의 strict / broad pseudo-group, `plot_immune_vs_epithelial()` (3-panel: immune / strict / broad)
- `inference/analysis/slide{1,2}_*_v2/findings.md` — 각 슬라이드별 두 score 의 해석
- `inference/analysis_filtered/COMPARISON.md` — filter 적용 시 두 score 비율 변화
- `inference/analysis_filtered/slide{1,2}_*_v2/findings.md` — filter 후 두 score 의 재해석
- `inference/analysis/README.md` — 분석 디렉토리 전체 caveat
- `report/04_WSI에서_분석까지_쿡북.md` — pipeline cookbook 의 score 언급
- `/home/sjhong/CUCA/HER2ST_VS_LUNG_MAPPING.md` — her2st 39 type 도착 후 직접 mammary epithelial 측정으로 전환하는 길

본 문서는 위 모든 파일의 **methodological 기반** 으로서, 후속 보고서/논문 작성 시 인용 단위로 사용 가능.

---

## 9. 이전 버전과의 차이 (transparency)

이전 버전 (`CANCER_PROXY_METHODOLOGY.md`, 삭제됨) 으로부터의 주요 변경:

| 항목 | 이전 버전 | 본 버전 |
|---|---|---|
| 명칭 | "cancer-proxy" | "epithelial-activity proxy" (strict / broad) |
| 5 label 의 위계 | 단일 score (균일 취급) | 두 score 병렬 (strict 3종 / broad 5종) |
| Suprabasal | core member (DCIS analog 주장) | broad only, sensitivity 항목 (단독 해석 금지) |
| AT2 → mammary luminal progenitor | "read 가능" 단정 | **검증 가설** 로 톤다운 |
| Basal → TNBC origin | "후보" 단정 | luminal progenitor 보고와의 충돌 명시, 단정 금지 |
| Wuidart 2018 인용 | adult mammary basal stem | **embryonic mammary lineage segregation** (정확한 frame) |
| self-renewal 일반화 | 5 label 전체 | strict 3 종으로 한정, Suprabasal 제외 |
| 외부 공유 가능성 | 표현 수위 과함 | 안전한 외부 reader 톤 |

해당 변경의 사유는 외부 검토 의견 (이 conversation 의 critique session) 을 반영. critique 의 핵심: "5 label 선정 자체는 가설-기반 타당하나 *cancer-proxy* 명칭과 mapping 표현 강도가 tumor detection 으로 오해될 수 있음". 본 버전은 이를 직접 반영.
