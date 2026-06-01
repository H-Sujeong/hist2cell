"""3x3 per-slide UMAP — prediction_log1p vs dino vs hex+dino.

가설: hex+dino 가 dino 보다 Hist2Cell cell-type 으로 더 뭉친다.

행 = 3 슬라이드, 열 = [prediction_log1p, dino, hex+dino(agg)].
색 = Hist2Cell dominant cell type (argmax prediction).
정량: dominant-type kNN purity (k=10, 같은 type 이웃 비율) — 높을수록 cell-type clustering.

전처리: 각 representation 을 per-dim z-score 후 UMAP. agg = [dino 768 ⊕ hex 19] 인데
두 블록 스케일이 ~100배 차이라 raw 면 hex 19-d 가 거리를 지배 (README §주의). z-score 로 보정.

agg 의 N 이 prediction N 과 다르면 그 슬라이드 hex+dino 패널은 생략(오류 표기).

Usage (224):
  .venv/bin/python lung_pilot/hex_compare.py \
    --infer-dir lung_pilot/inference_output \
    --dino-dir  /mnt/fileserver/lung_pilot/dino_output \
    --agg-dir   /mnt/fileserver/lung_pilot/dino_hex_agg \
    --out-dir   lung_pilot/hex_compare_224 --label 224
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.neighbors import NearestNeighbors

SLIDES = ["TCGA-05-4245-01A-01-BS1", "TCGA-05-4245-01A-01-TS1", "TCGA-05-4390-01A-01-BS1"]
COLS = ["prediction_log1p", "dino", "hex+dino"]
UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)


def short(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def zscore(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    return (X - mu) / sd


def knn_purity(X, labels, k=10):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, idx = nn.kneighbors(X)
    lab = np.asarray(labels)
    return float((lab[idx[:, 1:]] == lab[:, None]).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infer-dir", required=True)
    ap.add_argument("--dino-dir", required=True)
    ap.add_argument("--agg-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    INFER, DINO, AGG = Path(a.infer_dir), Path(a.dino_dir), Path(a.agg_dir)
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    EMB = OUT / "embeddings"; EMB.mkdir(exist_ok=True)

    data, all_dom = {}, []
    for s in SLIDES:
        pred = np.load(INFER / s / "predictions.npy")
        df = pd.read_csv(INFER / s / "predictions.csv")
        ct_cols = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
        dom = np.array([ct_cols[i] for i in pred.argmax(1)])
        dino = np.load(DINO / s / "features_dinov2.npy")
        agg = np.load(AGG / s / "features_agg.npy")
        ok = agg.shape[0] == len(pred)
        data[s] = dict(pred=pred, dom=dom, dino=dino, agg=agg if ok else None, ok=ok)
        all_dom += dom.tolist()
        print(f"{short(s)}: N={len(pred)} dino{dino.shape} agg{agg.shape} ok={ok}"
              f"{'' if ok else '  <-- N MISMATCH'}")

    cnt = Counter(all_dom)
    major = [ct for ct, c in cnt.most_common() if c >= 50]
    cmap = plt.get_cmap("tab20", max(20, len(major)))
    color_map = {ct: cmap(i % cmap.N) for i, ct in enumerate(major)}
    OTHER = (0.7, 0.7, 0.7, 1.0)

    def reps_for(s):
        d = data[s]
        return {
            "prediction_log1p": zscore(np.log1p(d["pred"])),
            "dino": zscore(d["dino"]),
            "hex+dino": zscore(d["agg"]) if d["ok"] else None,
        }

    def get_emb(slide, col, X):
        cache = EMB / f"{slide}_{col.replace('+', '_')}_umap2d.npy"
        if cache.exists():
            return np.load(cache)
        z = umap.UMAP(**UMAP_KW).fit_transform(X)
        np.save(cache, z)
        return z

    # kNN purity
    rows = []
    for s in SLIDES:
        reps = reps_for(s)
        for col in COLS:
            X = reps[col]
            rows.append({"slide": short(s), "rep": col,
                         "knn_purity_k10": None if X is None else round(knn_purity(X, data[s]["dom"]), 4),
                         "note": "agg N mismatch" if X is None else
                                 ("circular (label=argmax)" if col == "prediction_log1p" else "")})
    pur = pd.DataFrame(rows)
    pur.to_csv(OUT / "knn_purity.csv", index=False)
    print("\n== dominant-type kNN purity (k=10) ==")
    print(pur.to_string(index=False))

    # 3x3 figure
    fig, axes = plt.subplots(len(SLIDES), len(COLS), figsize=(6 * len(COLS), 6 * len(SLIDES)))
    for r, s in enumerate(SLIDES):
        reps = reps_for(s)
        dom = data[s]["dom"]
        pcs = [color_map.get(ct, OTHER) for ct in dom]
        for c, col in enumerate(COLS):
            ax = axes[r, c]
            X = reps[col]
            pr = pur[(pur.slide == short(s)) & (pur.rep == col)]["knn_purity_k10"].iloc[0]
            if X is None:
                ax.text(0.5, 0.5, "agg 파일 N 불일치\n재생성 필요", ha="center", va="center",
                        fontsize=14, color="crimson")
            else:
                z = get_emb(s, col, X)
                ax.scatter(z[:, 0], z[:, 1], c=pcs, s=4, alpha=0.6, linewidths=0)
            ax.set_title((col if r == 0 else "") + (f"\npurity={pr}" if pr is not None else ""),
                         fontsize=12)
            ax.set_xticks([]); ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(f"{short(s)} ({len(dom)})", fontsize=12)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=color_map[ct], markersize=8,
                          label=ct) for ct in major]
    handles.append(plt.Line2D([], [], marker="o", linestyle="", color=OTHER, markersize=8, label="other"))
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5), fontsize=9,
               title="Hist2Cell dominant cell type")
    fig.suptitle(
        f"[{a.label}] per-slide UMAP — prediction_log1p vs dino vs hex+dino   "
        f"(color=dominant cell type, per-dim z-scored, title=kNN purity k10)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 0.9, 0.98])
    p = OUT / "umap_3x3_pred_dino_hexdino.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"\nsaved: {p}\ndone.")


if __name__ == "__main__":
    main()
