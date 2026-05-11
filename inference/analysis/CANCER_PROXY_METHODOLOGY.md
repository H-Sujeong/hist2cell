# Cancer-proxy 5 type 선정의 근거 (lung-trained → breast 적용)

> **이 문서가 필요한 이유**
>
> 우리가 KBSMC breast 슬라이드 (slide1, slide2) 에 적용한 Hist2Cell 가중치 `humanlung_cell2location_leave_A50_out.pth` 는 **healthy human lung** 데이터로 학습된 모델이다. 80개 출력 cell type 은 모두 lung 분류 라벨이며 breast 에 직접 대응하는 type 이 없다. 그럼에도 우리는 `cell_type_groups.csv` 에서 **5 개 type** (`AT2`, `Basal`, `Suprabasal`, `Dividing_AT2`, `Dividing_Basal`) 을 `is_cancer_proxy=1` 로 marking 하여 "cancer-proxy" abundance 합으로 사용해왔다 (`inference/analysis/analyze.py`, `findings.md`, `cookbook.md`).
>
> 이 선정은 **임의 판단** 이 아니라 다음 두 가지 생물학적 가정에 기반하며, 본 문서는 각 type 별로 그 근거를 참고문헌과 함께 명시한다.
>
> 1. **Proliferative / progenitor epithelial signal 은 tissue-agnostic** — 분열·분화 중인 상피세포는 tissue-of-origin 과 무관하게 "활성 epithelial compartment" 의 spatial marker 로 기능한다.
> 2. **Lung 의 basal stem cell / alveolar progenitor 와 breast 의 basal / luminal progenitor 는 stem cell hierarchy 구조가 유사** — 두 조직 모두 basal layer (KRT5+/TP63+) + 분비/구조 layer (KRT8/18+) 의 hierarchy 를 가진다.

---

## 1. 선정된 5 type 과 lung-맥락 정의

| cell type | lineage (lung) | role in lung | 본 분석에서의 의미 |
|---|---|---|---|
| **AT2** | Epithelial-alveolar | type II pneumocyte; alveolar progenitor + surfactant producer | proliferative epithelial proxy (alveolar 계열) |
| **Dividing_AT2** | Epithelial-alveolar | 명시적으로 분열기에 있는 AT2 cell (Ki67+/MKI67+) | 분열 중 alveolar progenitor → 직접 proliferation marker |
| **Basal** | Epithelial-airway | airway basal stem cell (TP63+/KRT5+/KRT14+) | proliferative epithelial proxy (basal 계열) |
| **Suprabasal** | Epithelial-airway | basal → ciliated/secretory 분화 중간층 | 분화 시작 단계, 여전히 active epithelial compartment |
| **Dividing_Basal** | Epithelial-airway | 분열 중인 basal cell (Ki67+/cyclin+) | 분열 중 basal progenitor → 직접 proliferation marker |

5 type 의 공통점:
- 모두 **상피 (epithelial)** 계열
- 모두 **줄기/전구 (stem/progenitor) 또는 명시적 분열기** 특성
- 모두 **자기 갱신 (self-renewal) 능력 보유** — 정상 조직에서도, 종양에서도 활성

---

## 2. 개별 type 의 생물학적 근거

### 2.1 AT2 (alveolar type 2 cell)
- AT2 는 폐포 (alveoli) 의 **distal lung stem cell** 로, 정상 조직에서 자기 갱신 + AT1 분화를 모두 담당. SP-C / SP-B / ABCA3 발현.
- **Barkauskas et al. 2013** [1]: AT2 cell 이 polyclonal 으로 분열하며 AT1 으로 분화함을 lineage tracing 으로 증명. "AT2 cells are stem cells of the alveolus."
- **Kim et al. 2005** [2]: AT2 + bronchioalveolar duct junction (BADJ) 의 BASC (bronchioalveolar stem cell) 가 **lung adenocarcinoma (LUAD) 의 cell-of-origin**. K-Ras 변이로 in vivo 종양 유도.
- **Madissoon et al. 2023** [3]: 본 모델이 학습한 인체 폐 Cell Atlas 의 AT2 라벨 정의 — Cell2location deconvolution 으로 80 type 중 하나.

**breast 맥락에서의 read**: AT2 는 lung-specific 라벨이지만, 그 핵심 정의 ("자기 갱신 + 분비 능력 보유 + 분열성 progenitor") 는 **mammary luminal progenitor** [9] 와 기능적으로 평행 (둘 다 KRT8/18+, alveolar 분비, 분열성). 따라서 KBSMC 슬라이드의 AT2 spatial pattern 은 **mammary luminal progenitor 또는 일반적 분비-상피 progenitor** 의 신호로 read 가능.

### 2.2 Dividing_AT2
- AT2 의 부분집합 중 **명시적으로 mitotic** (Ki67+/cyclin-positive) 인 cell 만 별도 라벨링. Madissoon 2023 의 atlas 에서 cell-cycle gene 발현으로 sub-cluster 분리.
- AT2 의 분열은 종양 형성의 직접 전구체 — Kim 2005 [2] 가 K-Ras 활성화 시 AT2 가 분열 → BASC → LUAD 진행을 추적.

**breast 맥락**: Dividing_AT2 는 우리 모델 출력의 가장 직접적인 "spatial proliferation marker" 중 하나. breast 에서도 분열 중인 epithelial 영역과 spatial overlap 한다면 그것은 곧 **분열 활성 epithelial compartment** = potential tumor 영역.

### 2.3 Basal (airway basal cell)
- 기도 (airway) 상피의 **basal layer stem cell** (TP63+/KRT5+/KRT14+). 자기 갱신 + ciliated/secretory/goblet 으로 분화.
- **Rock et al. 2009** [4]: 마우스 trachea + 인간 기도 상피에서 basal cell 이 multipotent stem cell 임을 lineage tracing 으로 증명.
- **Sutherland & Berns 2010** [5]: Basal cell 이 **squamous cell carcinoma (LUSC, 폐 편평세포암)** 의 cell-of-origin.

**breast 맥락**: Lung airway basal cell 과 **mammary basal cell (myoepithelial-like)** 은 동일 marker (KRT5+/TP63+/KRT14+) 를 공유한다. **Wuidart et al. 2018** [6] 가 mammary basal 의 stem cell 활성을 증명. 따라서 breast 에서 lung Basal signal hot-spot 은 **mammary basal/myoepithelial compartment** (= basal-like / triple-negative breast cancer 의 cell-of-origin 후보) 의 spatial marker 로 read 가능.

### 2.4 Suprabasal
- Basal layer 위, ciliated/secretory 로 분화하기 직전의 중간 분화층. KRT5 발현 감소 + KRT8 발현 시작.
- Madissoon 2023 [3] atlas 에 별도 cluster. cell cycle gene 은 음성이지만 progenitor → differentiated 의 transitional state.
- 정상 조직에서는 short-lived intermediate. **염증/수복 (regenerative) 환경에서 Suprabasal abundance 가 늘어남** — 종양 미세환경의 wound-healing-like phenotype 신호.

**breast 맥락**: 직접 1:1 라벨 매칭은 없으나, 분화 중간 단계의 epithelial compartment 로 read. **DCIS (ductal carcinoma in situ) → invasive carcinoma 진행기** 의 transitional epithelial state 와 functional analog.

### 2.5 Dividing_Basal
- Basal 의 부분집합 중 명시적으로 mitotic (Ki67+) 인 cell. 정상에서도 항상 일정 비율 (1-5%) 존재.
- **종양 영역에서 abundance 가 증가** — basal stem cell 의 hyperproliferation 신호.

**breast 맥락**: Dividing_AT2 와 마찬가지로 직접 proliferation marker. mammary basal stem cell + cell cycle 활성 = **basal-like breast cancer 의 핵심 phenotype**.

---

## 3. Lung → Breast 매핑 요약

| lung label | 직접 매칭되는 breast label (CUCA her2st 39 type 기준 [11]) | 매핑 신뢰도 |
|---|---|---|
| AT2 / Dividing_AT2 | luminal epithelial cell of mammary gland (mammary luminal progenitor analog) | 중 (lineage 평행, type 1:1 아님) |
| Basal / Dividing_Basal | basal cell (mammary basal stem cell, KRT5+/TP63+ 공유) | 높음 (marker 공유) |
| Suprabasal | mammary gland epithelial cell (generic) | 낮 (정확한 1:1 없음) |

따라서 **본 5 type 합 ≈ her2st 의 3 종 mammary epithelial (basal + luminal + generic mammary) 합** 에 가까운 spatial signal 을 잡는다고 해석 가능 (단 lung-trained 모델의 feature space 가 breast 와 다르므로 정확한 정량 일치는 불가). 이는 우리 `findings.md` 와 `analyze.py` 의 cancer-proxy 분석이 **proliferative mammary epithelial compartment 의 spatial proxy** 라는 의미를 부여한다.

---

## 4. 명시적 한계 (re-iterated for emphasis)

1. **proxy ≠ detector**: 이 5 type 합은 종양 detection 이 아니라 *epithelial proliferation* 의 spatial reference. 진짜 종양 영역 검출은 별도 모델 (tiatoolbox 위험도, CUCA her2st, IHC) 가 담당.
2. **모델 feature space mismatch**: lung 조직으로 학습한 backbone (ResNet18) 의 representation 이 breast H&E 의 morphology 와 어떤 차이를 보이는지는 정량 검증 없음. 본 proxy 는 *상대 spatial pattern* 으로만 해석.
3. **mpp / tile_size mismatch**: Visium 20× (~0.5 μm/px) 학습 vs Aperio 40× (0.2615 μm/px) 적용. 모델 시야가 학습 분포의 절반 — 절대값 비교 금지.
4. **분리된 tissue blob 의 영향**: `inference/analysis_filtered/COMPARISON.md` 에서 보았듯, slide2 의 경우 cancer-proxy 우세 spot 17.7% 중 상당 부분이 가장 큰 덩어리 바깥에 위치. proxy abundance 의 절대값 / spatial 분포 해석 시 슬라이드 내 다중 compartment 가능성 명심.
5. **5 type 외 다른 progenitor 후보**: lung 80 type 중 본 5 종 외에도 `Macro_dividing`, `B_plasmablast`, `Deuterosomal` (ciliated-progenitor intermediate) 등 분열성 표현형이 있으나, 이들은 면역세포 또는 분비-airway 분화로 epithelial cancer compartment 의 spatial proxy 로는 부적합하여 제외.

---

## 5. 대안 — CUCA her2st 가중치 도착 후 권장 흐름

`/home/sjhong/CUCA/HER2ST_VS_LUNG_MAPPING.md` 와 함께 보면:

| 단계 | 본 lung-기반 분석 | her2st 가중치 도착 후 |
|---|---|---|
| cancer-proxy 정의 | 5 lung type 합 (AT2/Basal/Suprabasal/Dividing_*) | 3 mammary epithelial type 합 (basal/luminal/mammary gland) |
| `is_*` flag | `is_cancer_proxy=1` (5 종) | `is_mammary_epithelial=1` (3 종) → analyze.py 1줄 패치 필요 |
| 해석 | spatial proxy | 직접 측정 (proxy 가 아님) |
| 신뢰도 | 정성 cross-check 까지만 | 정량 ROI-level 분석 가능 |

→ 두 결과의 spatial overlap (Pearson ρ / Jaccard) 을 측정하면 본 lung-proxy 의 정성적 정확도가 직접 검증 가능. 후속 단계의 1순위.

---

## 6. 참고문헌

[1] **Barkauskas CE**, Cronce MJ, Rackley CR, et al. *Type 2 alveolar cells are stem cells in adult lung*. **J Clin Invest** 2013;123(7):3025–3036. doi:10.1172/JCI68782 (PMID: 23921127)

[2] **Kim CFB**, Jackson EL, Woolfenden AE, et al. *Identification of bronchioalveolar stem cells in normal lung and lung cancer*. **Cell** 2005;121(6):823–835. doi:10.1016/j.cell.2005.03.032

[3] **Madissoon E**, Oliver AJ, Kleshchevnikov V, et al. *A spatially resolved atlas of the human lung characterizes a gland-associated immune niche*. **Nat Genet** 2023;55:66–77. doi:10.1038/s41588-022-01243-4 *(본 모델 학습 데이터의 source atlas — 80 lung cell types via Cell2location)*

[4] **Rock JR**, Onaitis MW, Rawlins EL, et al. *Basal cells as stem cells of the mouse trachea and human airway epithelium*. **PNAS** 2009;106(31):12771–12775. doi:10.1073/pnas.0906850106

[5] **Sutherland KD**, Berns A. *Cell of origin of lung cancer*. **Mol Oncol** 2010;4(5):397–403. doi:10.1016/j.molonc.2010.05.002 *(review: AT2 → LUAD, Basal → LUSC)*

[6] **Wuidart A**, Sifrim A, Fioramonti M, et al. *Early lineage segregation of multipotent embryonic mammary gland progenitors*. **Nat Cell Biol** 2018;20:666–676. doi:10.1038/s41556-018-0095-2

[7] **Lim E**, Vaillant F, Wu D, et al. *Aberrant luminal progenitors as the candidate target population for basal tumor development in BRCA1 mutation carriers*. **Nat Med** 2009;15(8):907–913. doi:10.1038/nm.2000 *(breast luminal progenitor 가 basal-like TNBC 의 cell-of-origin)*

[8] **Visvader JE**. *Keeping abreast of the mammary epithelial hierarchy and breast tumorigenesis*. **Genes Dev** 2009;23(22):2563–2577. doi:10.1101/gad.1849509 *(mammary cell hierarchy review — basal vs luminal)*

[9] **Pal B**, Chen Y, Vaillant F, et al. *A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast*. **EMBO J** 2021;40(11):e107333. doi:10.15252/embj.2020107333 *(breast cell atlas — luminal progenitor + basal stem cell 정의)*

[10] **Sørlie T**, Perou CM, Tibshirani R, et al. *Gene expression patterns of breast carcinomas distinguish tumor subclasses with clinical implications*. **PNAS** 2001;98(19):10869–10874. doi:10.1073/pnas.191367098 *(luminal A/B / HER2 / basal-like 분자 분류)*

[11] **Andersson A**, Larsson L, Stenbeck L, et al. *Spatial deconvolution of HER2-positive breast cancer delineates tumor-associated cell type interactions*. **Nat Commun** 2021;12:6012. doi:10.1038/s41467-021-26271-2 *(her2st dataset 의 원본 — CUCA 가 학습에 사용)*

[12] **Travaglini KJ**, Nabhan AN, Penland L, et al. *A molecular cell atlas of the human lung from single-cell RNA sequencing*. **Nature** 2020;587:619–625. doi:10.1038/s41586-020-2922-4

[13] **Hist2Cell** (used model). Zhao W. et al. *Hist2Cell* GitHub: https://github.com/Weiqin-Zhao/Hist2Cell — paired histology + spatial transcriptomics → 80 lung cell type prediction. 본 분석의 가중치 `model_weights/humanlung_cell2location_leave_A50_out.pth` 는 이 repo 의 공식 배포본.

---

## 7. 본 문서와 연결된 파일

- `inference/analysis/cell_type_groups.csv` — `is_cancer_proxy` 컬럼 정의 위치 (line 3, 4, 7, 8, 18)
- `inference/analysis/analyze.py` — `per_group_stats()` 의 cancer-proxy pseudo-group, `plot_immune_vs_cancer()` 의 cancer_members 선택 로직
- `inference/analysis/slide{1,2}_*_v2/findings.md` — 각 슬라이드별 cancer-proxy abundance 해석
- `inference/analysis_filtered/COMPARISON.md` — filter 적용 시 cancer-proxy 비율 변화
- `inference/analysis_filtered/slide{1,2}_*_v2/findings.md` — filter 후 재해석
- `inference/analysis/README.md` — 분석 디렉토리 전체 caveat
- `report/04_WSI에서_분석까지_쿡북.md` — pipeline cookbook 의 cancer-proxy 언급 (§4.2 "cancer-proxy 5 type 자기상관")
- `/home/sjhong/CUCA/HER2ST_VS_LUNG_MAPPING.md` — her2st 39 type 도착 후 직접 mammary epithelial 측정으로 전환하는 길

본 문서는 위 모든 파일의 **methodological 기반** 으로서, 후속 보고서/논문 작성 시 인용 단위로 사용 가능.
