"""dino(768) 를 hex 와 동일한 19-dim 으로 축소 후 비교 — dino 차원 지배 문제 보수.

이전 hex_compare 에서 agg=dino768⊕hex19 의 per-dim concat 은 dino 768차원에 hex(2.4%)가
묻혔다. 여기선 dino 를 19-dim 으로 줄여 hex 와 동등 차원으로 맞춘 뒤 비교.

dino 19-dim 축소 2가지:
  - PCA      : dino 768 → 19 주성분 (unsupervised, top-variance 방향)
  - IMP      : dominant cell type 에 대한 ANOVA F-score 상위 19 원본 dim (supervised
               중요도; ⚠️ 라벨=Hist2Cell prediction 이라 dino 에 유리·circular)

라벨 = dominant cell type = argmax(Hist2Cell prediction). 모든 rep per-dim z-score 후
kNN purity(k=10) + per-slide UMAP.

Usage:
  .venv/bin/python lung_pilot/hex_dino19_compare.py \
    --infer-dir lung_pilot/inference_output \
    --dino-dir  /mnt/fileserver/lung_pilot/dino_output \
    --agg-dir   /mnt/fileserver/lung_pilot/dino_hex_agg \
    --out-dir   lung_pilot/hex_dino19_224 --label 224
"""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_classif
from sklearn.neighbors import NearestNeighbors

SLIDES = ["TCGA-05-4245-01A-01-BS1", "TCGA-05-4245-01A-01-TS1", "TCGA-05-4390-01A-01-BS1"]
UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)
N19 = 19


def short(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def zsc(X):  # per-dim z-score
    sd = X.std(0); sd[sd == 0] = 1.0
    return (X - X.mean(0)) / sd


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
    a = ap.parse_args()
    INFER, DINO, AGG = Path(a.infer_dir), Path(a.dino_dir), Path(a.agg_dir)
    OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
    EMB = OUT / "embeddings"; EMB.mkdir(exist_ok=True)

    # load combined
    pred_l, dino_l, hex_l, dom_l, slide_l = [], [], [], [], []
    for s in SLIDES:
        df = pd.read_csv(INFER / s / "predictions.csv")
        ct = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
        pred = np.load(INFER / s / "predictions.npy")
        dino = np.load(DINO / s / "features_dinov2.npy")
        agg = np.load(AGG / s / "features_agg.npy")
        assert agg.shape[0] == len(pred), f"{s} N mismatch"
        pred_l.append(pred); dino_l.append(dino); hex_l.append(agg[:, 768:])
        dom_l += [ct[i] for i in pred.argmax(1)]; slide_l += [s] * len(pred)
    pred = np.vstack(pred_l); dino = np.vstack(dino_l); hexb = np.vstack(hex_l)
    dom = np.array(dom_l); slide = np.array(slide_l)
    print(f"[{a.label}] combined N={len(dom)} dino{dino.shape} hex{hexb.shape}")

    # --- dino 19-dim 축소 (combined fit) ---
    dino_z = zsc(dino)
    pca = PCA(n_components=N19, random_state=42).fit(dino_z)
    dino19_pca = pca.transform(dino_z)
    evr = pca.explained_variance_ratio_.sum()
    F, _ = f_classif(dino_z, dom)
    top19 = np.argsort(-np.nan_to_num(F))[:N19]
    dino19_imp = dino_z[:, top19]
    print(f"  PCA19 explained variance = {evr:.3f}  | IMP top19 dims = {sorted(top19.tolist())}")

    # reps (combined, per-dim z-scored)
    REPS = {
        "prediction_log1p": zsc(np.log1p(pred)),
        "dino768": dino_z,
        "dino19_pca": zsc(dino19_pca),
        "dino19_imp": zsc(dino19_imp),
        "hex19": zsc(hexb),
        "dino19pca+hex19": zsc(np.hstack([dino19_pca, hexb])),
        "dino19imp+hex19": zsc(np.hstack([dino19_imp, hexb])),
    }

    # --- kNN purity per-slide ---
    rows = []
    for name, X in REPS.items():
        for s in SLIDES:
            mask = slide == s
            rows.append({"slide": short(s), "rep": name, "dim": X.shape[1],
                         "knn_purity_k10": round(knn_purity(X[mask], dom[mask]), 4)})
    pur = pd.DataFrame(rows)
    pur.to_csv(OUT / "knn_purity.csv", index=False)
    piv = pur.pivot(index="rep", columns="slide", values="knn_purity_k10")
    piv = piv.reindex(list(REPS.keys()))
    print("\n== kNN purity (k=10, dominant cell type) ==")
    print(piv.to_string())
    piv.to_csv(OUT / "knn_purity_pivot.csv")

    # --- UMAP per-slide, 핵심 컬럼 ---
    cols = ["prediction_log1p", "dino19_pca", "dino19_imp", "hex19", "dino19pca+hex19"]
    major = [c for c, n in Counter(dom).most_common() if n >= 50]
    cmap = plt.get_cmap("tab20", max(20, len(major)))
    cmd = {c: cmap(i % cmap.N) for i, c in enumerate(major)}; OTHER = (.7, .7, .7, 1.)
    fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(5.5 * len(cols), 5.6 * len(SLIDES)))
    for r, s in enumerate(SLIDES):
        mask = slide == s; pc = [cmd.get(c, OTHER) for c in dom[mask]]
        for c, name in enumerate(cols):
            ax = axes[r, c]
            cache = EMB / f"{short(s)}_{name.replace('+','_')}_umap.npy"
            X = REPS[name][mask]
            z = np.load(cache) if cache.exists() else umap.UMAP(**UMAP_KW).fit_transform(X)
            if not cache.exists(): np.save(cache, z)
            ax.scatter(z[:, 0], z[:, 1], c=pc, s=4, alpha=.6, linewidths=0)
            pr = pur[(pur.slide == short(s)) & (pur.rep == name)]["knn_purity_k10"].iloc[0]
            ax.set_title((name if r == 0 else "") + f"\npurity={pr}", fontsize=11)
            ax.set_xticks([]); ax.set_yticks([])
            if c == 0: ax.set_ylabel(f"{short(s)} ({mask.sum()})", fontsize=11)
    h = [plt.Line2D([], [], marker="o", ls="", color=cmd[c], ms=8, label=c) for c in major]
    h.append(plt.Line2D([], [], marker="o", ls="", color=OTHER, ms=8, label="other"))
    fig.legend(handles=h, loc="center right", bbox_to_anchor=(1, .5), fontsize=8,
               title="dominant cell type")
    fig.suptitle(f"[{a.label}] dino→19dim (PCA / importance) vs hex19 vs dino19+hex19  "
                 f"(color=dominant cell type, z-scored, title=kNN purity k10)", fontsize=13)
    fig.tight_layout(rect=[0, 0, .9, .98])
    p = OUT / "umap_dino19_vs_hex.png"
    fig.savefig(p, dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"\nsaved: {p}\ndone.  (PCA19 EVR={evr:.3f})")


if __name__ == "__main__":
    main()
