"""각 Hist2Cell dominant cell type 의 '대표 패치' 1장씩 — prediction_log1p centroid 최근접.

dominant cell type = argmax(prediction). 그 type 의 spot 들 중 prediction_log1p(80-d)
centroid 에 가장 가까운 spot = 가장 전형적인 예시. 그 패치를 SVS 에서 추출해 라벨 montage.

산출 (--out-dir):
  celltype_examples.png   — cell type 별 예시 패치 montage (라벨 = type, n, mean abundance)
  celltype_examples.csv   — cell_type, n_dominant, mean_abundance, slide, spot_id, x, y

Usage (224):
  .venv/bin/python lung_pilot/celltype_examples.py \
    --infer-dir lung_pilot/inference_output --graph-dir lung_pilot/graph_output/224 \
    --wsi-dir /mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi \
    --out-dir lung_pilot/celltype_examples_224 --label 224
"""
import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openslide
from PIL import Image

SLIDES = ["TCGA-05-4245-01A-01-BS1", "TCGA-05-4245-01A-01-TS1", "TCGA-05-4390-01A-01-BS1"]
DISP = 200


def short(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def find_svs(wsi, s):
    for pat in (f"{s}.*.svs", f"{s}.svs"):
        h = sorted(glob.glob(os.path.join(wsi, pat)))
        if h:
            return h[0]
    raise FileNotFoundError(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infer-dir", required=True)
    ap.add_argument("--graph-dir", required=True)
    ap.add_argument("--wsi-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    INFER, GRAPH = Path(a.infer_dir), Path(a.graph_dir)
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)

    pred_l, dom_l, sl_l, sid_l, x_l, y_l = ([] for _ in range(6))
    spots = {}
    ct_cols = None
    for s in SLIDES:
        df = pd.read_csv(INFER / s / "predictions.csv")
        ct_cols = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
        pred = np.load(INFER / s / "predictions.npy")
        spots[s] = pd.read_csv(GRAPH / f"{s}_spots.csv")
        pred_l.append(pred)
        dom_l += [ct_cols[i] for i in pred.argmax(1)]
        sl_l += [s] * len(pred); sid_l += list(df["spot_id"])
        x_l += list(df["X"]); y_l += list(df["Y"])
    pred = np.vstack(pred_l)
    logp = np.log1p(pred)
    dom = np.array(dom_l); sl = np.array(sl_l); sid = np.array(sid_l)
    x = np.array(x_l); y = np.array(y_l)
    ct_idx = {c: i for i, c in enumerate(ct_cols)}

    # dominant type 들, spot 수 desc
    from collections import Counter
    cnt = Counter(dom.tolist())
    types = [ct for ct, _ in cnt.most_common()]
    print(f"[{a.label}] distinct dominant types: {len(types)}")

    sl_h = {s: openslide.OpenSlide(find_svs(a.wsi_dir, s)) for s in SLIDES}
    rows, examples = [], []
    for ct in types:
        m = np.where(dom == ct)[0]
        cen = logp[m].mean(0)
        nearest = m[np.argmin(np.linalg.norm(logp[m] - cen, axis=1))]
        s = sl[nearest]
        sp = spots[s]
        row = sp[sp["spot_id"] == sid[nearest]].iloc[0]
        tlx, tly, ts = int(row["tile_x_topleft"]), int(row["tile_y_topleft"]), int(row["tile_size"])
        patch = sl_h[s].read_region((tlx, tly), 0, (ts, ts)).convert("RGB").resize((DISP, DISP))
        examples.append((ct, cnt[ct], float(pred[nearest, ct_idx[ct]]), np.asarray(patch)))
        rows.append({"cell_type": ct, "n_dominant": cnt[ct],
                     "abundance_at_example": round(float(pred[nearest, ct_idx[ct]]), 3),
                     "slide": short(s), "spot_id": sid[nearest],
                     "x": int(x[nearest]), "y": int(y[nearest])})
    for h in sl_h.values():
        h.close()
    pd.DataFrame(rows).to_csv(OUT / "celltype_examples.csv", index=False)

    # montage
    n = len(examples)
    ncol = 5
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.3 * nrow))
    axes = np.array(axes).reshape(-1)
    for ax, (ct, k, ab, img) in zip(axes, examples):
        ax.imshow(img)
        ax.set_title(f"{ct}\nn_dom={k}, abund={ab:.1f}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"[{a.label}] Hist2Cell dominant cell type 별 대표 패치 "
                 f"(prediction_log1p centroid 최근접) — {n} types", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = OUT / "celltype_examples.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"saved: {p}\nsaved: {OUT/'celltype_examples.csv'}\ndone.")


if __name__ == "__main__":
    main()
