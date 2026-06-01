"""hex19 패치별 하위 N% 값을 0으로(노이즈 제거) 후 dino19(PCA) 비교.

hex19 = agg[:,768:] (19 marker). 각 행(패치)에서 그 행의 N-percentile 미만 값을 0 으로 치환.
(N=10 → 19개 중 ~최저 2개가 0). 목적: 패치별 저신호 marker 를 노이즈로 보고 제거.

컬럼 배열: hist2cell(prediction) · PCA(dino19_pca) · hex(denoised) · pca+hex(denoised). (imp 제외)
denoise 효과 확인용으로 purity 표엔 원본 hex19 도 포함.

Usage:
  .venv/bin/python lung_pilot/hex_dino19_denoise.py \
    --infer-dir lung_pilot/inference_output --dino-dir /mnt/fileserver/lung_pilot/dino_output \
    --agg-dir /mnt/fileserver/lung_pilot/dino_hex_agg --out-dir lung_pilot/hex_dino19_denoise_224 --label 224
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

SLIDES = ["TCGA-05-4245-01A-01-BS1", "TCGA-05-4245-01A-01-TS1", "TCGA-05-4390-01A-01-BS1"]
UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)
N19 = 19


def short(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def zsc(X):
    sd = X.std(0); sd[sd == 0] = 1.0
    return (X - X.mean(0)) / sd


def denoise_rows(H, pct):
    thr = np.percentile(H, pct, axis=1, keepdims=True)
    return np.where(H < thr, 0.0, H)


def knn_purity(X, lab, k=10):
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    lab = np.asarray(lab)
    return float((lab[idx[:, 1:]] == lab[:, None]).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infer-dir", required=True)
    ap.add_argument("--dino-dir", required=True)
    ap.add_argument("--agg-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--pct", type=float, default=10.0, help="하위 몇 % 를 0 으로")
    a = ap.parse_args()
    INFER, DINO, AGG = Path(a.infer_dir), Path(a.dino_dir), Path(a.agg_dir)
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    EMB = OUT / "embeddings"; EMB.mkdir(exist_ok=True)

    pred_l, dino_l, hex_l, dom_l, slide_l = [], [], [], [], []
    for s in SLIDES:
        df = pd.read_csv(INFER / s / "predictions.csv")
        ct = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
        pred = np.load(INFER / s / "predictions.npy")
        dino = np.load(DINO / s / "features_dinov2.npy")
        agg = np.load(AGG / s / "features_agg.npy")
        assert agg.shape[0] == len(pred)
        pred_l.append(pred); dino_l.append(dino); hex_l.append(agg[:, 768:])
        dom_l += [ct[i] for i in pred.argmax(1)]; slide_l += [s] * len(pred)
    pred = np.vstack(pred_l); dino = np.vstack(dino_l); hexb = np.vstack(hex_l)
    dom = np.array(dom_l); slide = np.array(slide_l)

    hex_dn = denoise_rows(hexb, a.pct)
    zeroed = (hex_dn == 0).sum(1).mean()
    print(f"[{a.label}] N={len(dom)} | hex19 행당 0 처리 개수 평균={zeroed:.2f}/19 (pct={a.pct})")

    dino_z = zsc(dino)
    dino19_pca = PCA(N19, random_state=42).fit_transform(dino_z)

    REPS = {
        "hist2cell": zsc(np.log1p(pred)),
        "dino768": dino_z,
        "PCA": zsc(dino19_pca),                                   # dino19_pca
        "hex_orig": zsc(hexb),
        "hex": zsc(hex_dn),                                       # denoised hex19
        "pca+hex_orig": zsc(np.hstack([dino19_pca, hexb])),
        "pca+hex": zsc(np.hstack([dino19_pca, hex_dn])),          # denoised
    }

    rows = []
    for name, X in REPS.items():
        for s in SLIDES:
            m = slide == s
            rows.append({"slide": short(s), "rep": name, "dim": X.shape[1],
                         "knn_purity_k10": round(knn_purity(X[m], dom[m]), 4)})
    pur = pd.DataFrame(rows); pur.to_csv(OUT / "knn_purity.csv", index=False)
    piv = pur.pivot(index="rep", columns="slide", values="knn_purity_k10").reindex(list(REPS.keys()))
    print("\n== kNN purity (k=10, dominant cell type) ==\n", piv.to_string())
    piv.to_csv(OUT / "knn_purity_pivot.csv")

    cols = ["hist2cell", "PCA", "hex", "pca+hex"]   # 요청 배열
    major = [c for c, n in Counter(dom).most_common() if n >= 50]
    cmap = plt.get_cmap("tab20", max(20, len(major)))
    cmd = {c: cmap(i % cmap.N) for i, c in enumerate(major)}; OTHER = (.7, .7, .7, 1.)
    fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(5.6 * len(cols), 5.6 * len(SLIDES)))
    for r, s in enumerate(SLIDES):
        m = slide == s; pc = [cmd.get(c, OTHER) for c in dom[m]]
        for c, name in enumerate(cols):
            ax = axes[r, c]
            cache = EMB / f"{short(s)}_{name.replace('+','_')}_umap.npy"
            z = np.load(cache) if cache.exists() else umap.UMAP(**UMAP_KW).fit_transform(REPS[name][m])
            if not cache.exists(): np.save(cache, z)
            ax.scatter(z[:, 0], z[:, 1], c=pc, s=4, alpha=.6, linewidths=0)
            pr = pur[(pur.slide == short(s)) & (pur.rep == name)]["knn_purity_k10"].iloc[0]
            ax.set_title((name if r == 0 else "") + f"\npurity={pr}", fontsize=11)
            ax.set_xticks([]); ax.set_yticks([])
            if c == 0: ax.set_ylabel(f"{short(s)} ({m.sum()})", fontsize=11)
    h = [plt.Line2D([], [], marker="o", ls="", color=cmd[c], ms=8, label=c) for c in major]
    h.append(plt.Line2D([], [], marker="o", ls="", color=OTHER, ms=8, label="other"))
    fig.legend(handles=h, loc="center right", bbox_to_anchor=(1, .5), fontsize=8, title="dominant cell type")
    fig.suptitle(f"[{a.label}] hex19 하위 {a.pct:.0f}% →0 (denoise) | hist2cell · PCA(dino19) · hex · pca+hex  "
                 f"(color=dominant cell type, z-scored, title=kNN purity k10)", fontsize=13)
    fig.tight_layout(rect=[0, 0, .9, .98])
    p = OUT / "umap_denoise.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"\nsaved: {p}\ndone.")


if __name__ == "__main__":
    main()
