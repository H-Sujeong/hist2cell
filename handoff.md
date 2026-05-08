# Handoff — 2장 SVS 추론 작업 (Docker 재시작 후 재개)

세션 종료 컨텍스트: 2026-05-08, 사용자가 GPU 활성화를 위해 Docker 재시작 예정.

---

## 0. 이 작업이 무엇인지 한 줄

KBSMC 의 **breast** SVS 2장을 Hist2Cell 의 **lung** 가중치로 추론하여 (조직 mismatch 인 점을 알고 sanity-check 목적) `inference/` 에 결과를 정리하고 `report/` 에 전처리/추론 결과 md 작성.

---

## 1. 입력 슬라이드 (breast 조직)

```
/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-085-12,.svs
/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-152-19,,dup1.svs
```

파일명에 **공백/쉼표** 가 있으므로 모든 CLI 호출 시 큰따옴표로 감쌀 것.
파일 크기 ≈ 2.5GB / 2.8GB. SVS pyramidal TIFF.

## 2. 사용할 가중치 (lung — 의도적 mismatch)

```
model_weights/humanlung_cell2location_leave_A50_out.pth
```
- 출력 80 cell type 은 **폐 전용** (AT1/AT2/Basal/SMG_*/Schwann/Secretory_Club/Goblet 등)
- breast 슬라이드에 적용하면 **수치는 나오지만 생물학적 의미 없음** — 사용자가 sanity-check 목적임을 명시.
- 보고서 첫 단락에 이 caveat 반드시 박스로 강조.

대안 가중치: `humanlung_cell2location_leave_A37_out.pth`, `demo_ckpt.pth` (둘 다 lung 학습본). 둘 중 하나로 가도 무방.

---

## 3. Docker 재시작 직후 1회 실행할 환경 셋업

작업 디렉터리 `/home/sjhong/hist2cell` 기준.

```bash
# 1) GPU 살아있는지 먼저 확인 — 살아있어야 그 아래로 진행
python3 -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.device_count())"

# 2) openslide 시스템 라이브러리 + 파이썬 바인딩 (사용자가 OK 한 상태)
sudo apt-get update && sudo apt-get install -y openslide-tools libopenslide0
pip install openslide-python==1.3.1

# 3) 동작 확인
python3 -c "import openslide; print('openslide', openslide.__version__)"
```

GPU가 안 보이면 사용자에게 다시 알리고 멈출 것 — CPU 추론은 두 슬라이드 합쳐서 매우 길어질 수 있음.

---

## 4. 슬라이드 mpp / dimensions 확인 (spot-distance 정하기 위해)

```bash
python3 -c "
import openslide
for p in ['/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-085-12,.svs',
         '/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-152-19,,dup1.svs']:
    s = openslide.OpenSlide(p)
    print(p)
    print('  dims:', s.dimensions, 'levels:', s.level_count, 'level_dims:', s.level_dimensions)
    print('  mpp_x:', s.properties.get('openslide.mpp-x'),
          'mpp_y:', s.properties.get('openslide.mpp-y'),
          'objective:', s.properties.get('openslide.objective-power'))
    s.close()
"
```

**`--spot-distance` 매핑 룰** (학습 분포 ≈ 100μm):

| mpp (μm/px) | objective | `--spot-distance` |
|---|---|---|
| ~0.25 | 40× | `400` |
| ~0.50 | 20× | `200` |
| ~1.00 | 10× | `100` |

KBSMC SVS 는 보통 40× scan → `400` 일 가능성 큼. 실제 mpp 보고 결정.

---

## 5. 전처리 (각 슬라이드별로 1회)

`prep/prepare_wsi_for_inference.py` 는 이미 작성·검증 완료. 출력 디렉터리 구조:

```
inference/
├── slide1_085_12/
│   ├── slide1_085_12.pt          # PyG Data — 모델 입력
│   ├── spots.csv
│   ├── tissue_mask.png
│   └── spot_view.jpg
├── slide2_152_19/
│   └── ...
```

권장 명령:

```bash
SD=400   # 위에서 mpp 확인 후 결정한 값

python3 prep/prepare_wsi_for_inference.py \
  --input "/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-085-12,.svs" \
  --output ./inference/slide1_085_12 \
  --slide-name slide1_085_12 \
  --spot-distance $SD \
  --min-tissue-frac 0.4

python3 prep/prepare_wsi_for_inference.py \
  --input "/mnt/fileserver/Pathology/KBSMC/meteo_biotech_analysis_wsi/Z 2025000042,1-152-19,,dup1.svs" \
  --output ./inference/slide2_152_19 \
  --slide-name slide2_152_19 \
  --spot-distance $SD \
  --min-tissue-frac 0.4
```

> spot 수가 비정상적으로 많거나(>20k) 적으면(<200) `--spot-distance` 와 `--min-tissue-frac` 을 조정. `spot_view.jpg` 를 사용자에게 한 번 보여주고 컨펌 받는 것이 안전(전처리 caveat: §6).

---

## 6. 사용자에게 다시 컨펌 받아야 할 지점

전처리 직후·추론 직전에 한 번 멈출 것. 다음을 함께 보고:
1. 두 슬라이드의 **spot 수, edge 수** — 비정상 비대/희박 여부
2. **`spot_view.jpg`** 두 장 — tissue 위에 spot 이 잘 깔렸는지 (배경 들어왔거나 조직 빠뜨렸는지)
3. **`--spot-distance` 값** — mpp 보고 정한 값이 합당한지

문제 있으면 파라미터 조정해 prep 만 다시 돌리고, OK 받으면 §7 추론으로 넘어감.

---

## 7. 추론 코드 (작성 필요 — 아직 미작성)

`inference/infer.py` 를 신규 작성. 다음 책임:

1. `Hist2Cell` 클래스 정의 (training 노트북 cell 5 의 정의 그대로 — `cell_dim=80, vit_depth=3`)
2. 가중치 로드: `humanlung_cell2location_leave_A50_out.pth`
3. `Data(.pt)` 로드 → `NeighborLoader(num_neighbors=[-1,-1], batch_size=16, directed=False, shuffle=False)`
4. center node 단위로 fused head 출력 (cell_dim=80) 수집 → `[N, 80]`
5. `cell_types.pkl` 로 column 이름 매핑
6. 출력:
   - `inference/<slide>/predictions.csv` (`spot_id, X, Y, AT1, AT2, ..., gdT`)
   - `inference/<slide>/predictions.npy` (raw `[N, 80]` float32)
   - 옵션: 상위 5종 cell type 의 spot scatter 시각화 → `inference/<slide>/cell_<name>.png`

학습 노트북에서 모델은 (spot_pred, local_pred, fused_pred) 멀티헤드를 반환했음 — **fused 가 최종**. 단일 텐서를 반환하는 wrapper 가 있는지 model 코드 재확인 필요.

```python
# 큰 그림 (의사 코드)
model = Hist2Cell(cell_dim=80).to(device).eval()
model.load_state_dict(torch.load(weight_path, map_location=device))

data = torch.load(pt_path, weights_only=False)
loader = NeighborLoader(data, num_neighbors=[-1,-1], batch_size=16,
                        directed=False, shuffle=False, input_nodes=None)

preds = torch.zeros(data.num_nodes, 80)
with torch.no_grad():
    for sub in loader:
        sub = sub.to(device)
        out = model(sub.x, sub.edge_index)        # 학습 코드에서 fused head 가 어떤 인덱스인지 확인
        center = sub.input_id
        preds[center] = out_fused[: sub.batch_size].cpu()
```

> training tutorial 의 forward 마지막 부분을 다시 읽어서 `fused_pred` 가 단독 반환인지, tuple 의 몇 번째인지 확정 필요. (cell 5 일부만 보고 종료해서 미확정.)

---

## 8. 최종 산출물 정리 (사용자 요구)

```
inference/
├── slide1_085_12/
│   ├── slide1_085_12.pt
│   ├── spots.csv
│   ├── tissue_mask.png
│   ├── spot_view.jpg
│   ├── predictions.csv
│   ├── predictions.npy
│   └── (optional) cell_*.png
└── slide2_152_19/
    └── ...

report/
└── 03_breast슬라이드2장_lung가중치_추론결과_및_caveat.md
```

`report/03_*.md` 가 다뤄야 할 것:
- breast×lung mismatch caveat (맨 위, 강조 박스)
- 두 슬라이드의 mpp/dimensions/spot 수/edge 수 표
- `--spot-distance`, `--min-tissue-frac` 선택 근거
- spot_view 미니 첨부 (경로 링크)
- predictions 분포 요약 (각 slide 별 80 종 column means/maxes top-10 정도)
- 결론: 수치가 그럴듯해 보여도 **biological interpretation 금지**, lung-trained 출력 그대로임을 다시 명시

---

## 9. 이미 마친 것 / 안 한 것

### ✅ Done
- `prep/prepare_wsi_for_inference.py` (예시 슬라이드로 smoke test 통과: 169 nodes / 1051 edges, edges/node ≈ 6.22)
- `report/01_원본_WSI전처리파이프라인_튜토리얼이전단계_요약.md`
- `report/02_추론용_WSI전처리_prep스크립트_사용법및절차.md`

### ⏳ Pending (재개 시 할 일)
1. openslide 시스템+파이썬 설치 (§3)
2. 두 SVS 의 mpp 확인 → `--spot-distance` 결정 (§4)
3. `prep/prepare_wsi_for_inference.py` 두 번 실행 (§5)
4. **사용자 컨펌**: spot_view.jpg 보여주고 OK 받기 (§6)
5. `inference/infer.py` 작성 (§7) — Hist2Cell forward 의 fused head 인덱스 재확인
6. 두 슬라이드 추론 → predictions.csv/npy 저장 (§7)
7. `report/03_*.md` 작성 (§8)

---

## 10. 새 세션 첫 발화 권장문

> "handoff.md 따라서 이어서 진행. 지금 step 1 (환경 셋업) 시작."

이렇게 시작하면 됨. 컨텍스트 부족하면 `report/01`, `report/02`, `prep/prepare_wsi_for_inference.py` 를 우선 읽을 것.
