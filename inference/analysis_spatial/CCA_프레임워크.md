# CCA 프레임워크 — 본 프로젝트 진행 방식

`proof_ver2/` 의 cross-modality 검증에 사용한 CCA (Canonical Correlation
Analysis) 파이프라인의 단계·선택·해석 기준을 한 곳에 정리한 방법론 문서.
두 슬라이드(1_085_12, 1_152_19)에 동일하게 적용된다.

코드 본체: [`_proof_ver2_lib.py`](_proof_ver2_lib.py) (함수 `run_cca`, `permutation_null`)
호출 스크립트: 각 슬라이드 `proof_ver2/core_proofs_v2.py`

---

## 1. 무엇을 알고 싶은가

> ROI 1개당 두 종류의 측정값이 있다: Hist2Cell (80 cell-type 점수) 과
> proteomics (수천 개 단백질). **두 modality 가 ROI 들의 변동을 *같은
> 방향* 으로 설명하는가?**

대답할 수 있는 가장 자연스러운 통계 도구가 CCA.

---

## 2. CCA 의 본질 — 한 줄 설명

> 두 변수군 X (n×p), Y (n×q) 가 주어졌을 때, **각 변수군에서 선형결합을
> 한 번씩 골라 만든 두 합성치 사이의 Pearson 상관계수가 최대가 되도록**
> 결합 계수를 학습한다. 이때 만들어지는 합성치 쌍을 *canonical pair*,
> 그 상관계수를 *canonical correlation r* 이라고 한다. 첫 번째 쌍은
> 가장 강한 r, 두 번째 쌍은 첫 번째에 직교하면서 다음으로 큰 r ….

→ PCA 가 "분산 최대화" 라면 CCA 는 "**두 변수군 사이 상관 최대화**".

핵심 주의:
- CCA r 은 *측정* 이 아니라 *fitting 결과*. 자유도 (변수 수) 가 크면 신호가
  없어도 r 이 위쪽으로 inflated. → 본 문서 §6 의 inflation 항목 참조.

---

## 3. 본 프로젝트의 파이프라인 (6단계)

### 단계 1 — ROI 정렬 (`align_modalities`)

- 각 슬라이드에서 Hist2Cell tube ID 와 proteomics sample ID 가 *공통*
  으로 존재하는 ROI 만 추림.
- slide1 → 46 ROI / slide2 → 48 ROI.
- 추후 모든 분석은 이 공통 ROI 행렬에서만 수행.

산출: `H ∈ R^{n×80}`, `P ∈ R^{n×G}` (G = proteomics gene 개수 detect ≥ 50% 필터 후).

### 단계 2 — 표준화 (`StandardScaler`)

- 각 column 을 평균 0, 표준편차 1 로 정규화. modality 별로 따로.
- 이유: PCA 가 분산 큰 column 에 끌려가는 걸 막기 위해. proteomics 의
  log2 intensity 와 Hist2Cell 의 cell-type 점수는 스케일이 다르므로
  정규화 없이는 비교가 의미 없음.

### 단계 3 — PCA 차원 축소 (`PCA(n_components=10)`)

- H 와 P 각각 독립적으로 10 PC 까지 줄임.
- 이유: **CCA 의 inflation 을 줄이기 위함**. 원본 80·G 차원을 그대로
  CCA 에 넣으면 자유도가 너무 커서 ROI 가 46개뿐인 상황에서 trivial 한
  r=1 이 나옴 (perfect overfit).
- 10 이라는 숫자는 관행적 선택: `min(N-1, p, q)` 까지 가능하지만, ROI 수
  대비 너무 크지 않게 한 자릿수 후반에서 잘랐다. slide1·slide2 모두 PC10
  까지 누적 분산이 H ≈ 95%+, P ≈ 70%+ 수준이라 정보 손실은 작은 편.
- 코드 상수: `_proof_ver2_lib.py` 의 `N_PCS = 10`.

산출: `H_pcs ∈ R^{n×10}`, `P_pcs ∈ R^{n×10}`.

### 단계 4 — CCA fitting (`sklearn.cross_decomposition.CCA`)

- `CCA(n_components=3, max_iter=1000)`.
- `Hc, Pc = cca.fit_transform(H_pcs, P_pcs)` → 각 ROI 마다 3개의
  Hist2Cell canonical score 와 3개의 proteomics canonical score.
- 각 axis 의 train r 은 `pearsonr(Hc[:,i], Pc[:,i])` 로 계산.
- canonical loadings 는 원본 (PCA 이전) 공간으로 역추적:
  - `h_loadings = pca_h.components_.T @ cca.x_weights_` ∈ R^{80×3}
  - `p_loadings = pca_p.components_.T @ cca.y_weights_` ∈ R^{G×3}
  - 즉 각 cell type / gene 이 canonical axis 에 얼마나 기여하는지의 가중치.

**3개로 자른 이유는 코드 default**. PCA dim 까지 (최대 10개) 다 뽑을 수
있고, axis 별 train r 의 drop-off 와 null floor 를 보고 "유효 차원" 을
데이터 기반으로 정하는 것이 더 엄밀하나, 본 프로젝트에서는 axis 1 의
신호 유무를 보는 것이 주된 관심이라 3개에서 끊었다.

### 단계 5 — Permutation null (`permutation_null`)

- proteomics 행렬의 ROI 축을 1000회 무작위로 섞고, 매번 단계 2~4 를
  처음부터 다시 돌려서 top canonical r 의 영가설 분포를 만듦.
- empirical p (양측) = `mean(|null| ≥ |observed|)`.
- 코드 상수: `N_PERM = 1000`.
- 결과:
  - slide1: 관측 +0.936 vs null 평균 +0.778, 95% 상한 +0.863, p=0/1000
  - slide2: 관측 +0.940 vs null 평균 +0.768, 95% 상한 +0.857, p=0/1000

### 단계 6 — Loadings 해석

axis 1 의 ± top loaders 그림 (`cca_loadings_axis1.png`) 으로 *각 axis 가
어떤 조직 모듈에 해당하는지* 생물학적 의미 부여.
- 양 방향 loader 들: 모듈 A 의 구성요소.
- 음 방향 loader 들: 모듈 B (반대 극단) 의 구성요소.

slide1 axis 1 의 경우 + 방향 = 상피·선조직 (KRT8/18, B_plasma_IgA, AT2 등),
− 방향 = 면역·기질·혈관 (CD45=PTPRC, COL1A1, Muscle_smooth, Endothelia)
으로 해석.

---

## 4. 파라미터 요약

| 항목 | 값 | 코드 위치 |
|---|---|---|
| 표준화 방법 | z-score (mean 0, std 1) | `StandardScaler()` |
| PCA 차원 (modality 별) | 10 | `_proof_ver2_lib.py: N_PCS = 10` |
| CCA canonical pair 수 | 3 | `_proof_ver2_lib.py: N_CCA_COMP = 3` |
| Permutation 반복 | 1000 | `_proof_ver2_lib.py: N_PERM = 1000` |
| Random seed | 42 | `_proof_ver2_lib.py: RANDOM_SEED = 42` |

라이브러리 버전:
- scikit-learn `CCA` (`cross_decomposition`), `PCA`, `StandardScaler`
- numpy, pandas, scipy.stats.pearsonr

---

## 5. 산출물 매핑

`proof_ver2/` 안의 다음 파일들이 위 단계의 결과:

| 파일 | 단계 | 내용 |
|---|---|---|
| `cca_summary.csv` | 4·5 | axis 1~3 train r + null 95% 구간 + p |
| `cca_scatter.png` | 4 | 3 canonical pair 산점도 (section 색) |
| `permutation_null.png` | 5 | null 히스토그램 + 관측 r 선 |
| `cca_loadings_axis1.png` | 6 | axis 1 의 ± top loaders (H2C celltype / proteomics gene) |

---

## 6. ⚠ 핵심 caveat — CCA 의 inflation 문제

> **본 프레임워크가 보고하는 CCA r 의 절대값(예: +0.94)을 그대로 "두
> modality 가 94% 닮았다" 로 읽으면 안 된다.**

### 왜인가
CCA 는 *상관 최대화 fitting* 이므로, 변수 수에 비해 ROI 가 부족하면
*신호가 없어도* r 이 자동으로 위쪽으로 올라간다. 본 프로젝트에서:

- ROI N = 46~48 (작음)
- PCA dim 10 + 10 = 자유도가 N 에 비해 큼

이 조건에서 1000회 permutation 으로 측정한 영가설 평균이 **+0.78**.
즉 데이터에 진짜 신호가 없어도 자연히 그 정도까지 나온다.

### 어떻게 읽어야 하나

| 잘못된 읽기 | 올바른 읽기 |
|---|---|
| "관측 r = +0.94, 두 modality 가 강하게 일치" | "관측 r = +0.94 는 영가설 평균 +0.78 대비 +0.16 의 신호" |
| "axis 1 의 r 이 0.94 라서 거의 완전히 같다" | "관측치가 영가설 95% 구간 바깥이라 *방향성* 은 유의" |
| absolute magnitude 비교 | null 분포 대비 *상대 위치* 만 사용 |

### 다른 검정과의 비교

| 통계 | 영가설 평균 | 관측치 | 신호 = 관측 − null 평균 |
|---|---|---|---|
| **CCA top r (slide1)** | **+0.78** | +0.94 | +0.16 |
| Mantel Pearson r (slide1) | +0.00 | +0.20 | **+0.20** |
| Per-ROI cosine | (미검증, post-hoc bias) | +0.555 | — |

→ "정직한" 신호 크기는 CCA 와 Mantel 이 거의 비슷 (각각 +0.16, +0.20).
CCA 의 +0.94 라는 큰 숫자는 *통계의 inflation 가 만든 시각적 환상* 에
가깝다는 점을 본 프레임워크 사용 시 명시할 것.

---

## 7. 재현 방법

```bash
cd /home/sjhong/hist2cell/inference/analysis_spatial
# 슬라이드별
python 1_085_12/proof_ver2/core_proofs_v2.py
python 1_152_19/proof_ver2/core_proofs_v2.py
```

각 스크립트가 위 단계 1~6 을 순차 실행하고 `proof_ver2/` 폴더에
산출물을 떨군다. random seed 가 박혀있어 결과 재현 가능.

---

## 8. 같이 봐야 하는 다른 문서

- `*/proof_ver2/summary.md` — 슬라이드별 CCA 결과 + 해석
- `*/proof_ver2/연관성_분석.md` — CCA 외 다른 검정 (UMAP, silhouette,
  Mantel) 까지 종합한 보조 보고서. CCA inflation 의 맥락을 비교 검정으로
  보완하는 역할.

---

## 9. 결론적 사용 가이드

본 CCA 프레임워크는 **"두 modality 가 *어떤 공통 잠재 축* 을 공유하는지
빠르게 확인하기 위한 도구"** 로 사용한다. 강점·약점:

- ✓ 잠재 축의 *방향성* 과 *loadings* 해석을 깔끔하게 준다 (axis 1 = 상피
  vs 면역, 같은 식).
- ✓ Permutation null 과 같이 보면 통계적 유의성 검정 가능.
- ✗ 절대 r 값을 "강도" 로 해석할 수 없음 — null 대비 상대 위치만 사용.
- ✗ 정직한 cross-modality 일치 강도를 묻는다면 **Mantel test** 가 더
  적합 (`연관성_분석.md` 의 §4 참조).
