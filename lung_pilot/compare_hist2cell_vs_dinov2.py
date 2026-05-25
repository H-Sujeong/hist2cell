"""Hist2Cell ↔ DINOv2 정량 비교 (cross-slide 단위, raw representation).

산출 (lung_pilot/compare_output/):
  metrics.csv        — 모든 정량 metric 표
  metrics_bars.png   — bar chart 시각화
  summary.md         — 해석

Metric 3종:
  1) slide 1-NN purity  — 각 spot 의 1-NN 이 같은 슬라이드인 비율
                          (낮을수록 batch-mix 좋음, 4 rep 비교).
  2) kNN overlap (Jaccard) @ k=10, k=50 — 같은 spot 의 top-k 이웃 집합이
                          두 representation 에서 얼마나 일치하나
                          (Hist2Cell 3 rep × DINOv2).
  3) silhouette by lineage — 각 rep 에서 dominant cell-type lineage 가
                          얼마나 분리되나 (높을수록 cell-type 신호 보유).
                          시간 절약을 위해 5,000 spot sample.

전처리:
  - prediction 은 log1p 적용 (umap_compare.py 와 동일).
  - 그 외 rep 은 그대로.
  - 모든 metric 은 cross-slide 통합 (3 슬라이드 모든 spot) 에서 계산.

Usage:
    .venv/bin/python lung_pilot/compare_hist2cell_vs_dinov2.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parent.parent
INFER_DIR = ROOT / "lung_pilot" / "inference_output"
DINO_DIR = ROOT / "lung_pilot" / "dino_output"
GROUPS_CSV = ROOT / "inference" / "analysis_spatial" / "cell_type_groups.csv"
OUT_DIR = ROOT / "lung_pilot" / "compare_output"

SLIDES = [
    "TCGA-05-4245-01A-01-BS1",
    "TCGA-05-4245-01A-01-TS1",
    "TCGA-05-4390-01A-01-BS1",
]

REPS = [
    ("prediction_log1p", lambda d: np.log1p(d["pred"]), 80),
    ("features_fused",   lambda d: d["ff"],            256),
    ("features_resnet",  lambda d: d["fr"],            512),
    ("features_dinov2",  lambda d: d["dn"],            768),
]

K_VALUES = (10, 50)
SILHOUETTE_SAMPLE = 5000
SEED = 42


def load_slide(s):
    df = pd.read_csv(INFER_DIR / s / "predictions.csv")
    ct_cols = [c for c in df.columns if c not in ("spot_id", "X", "Y")]
    return dict(
        df=df,
        ct_cols=ct_cols,
        pred=np.load(INFER_DIR / s / "predictions.npy"),
        ff=np.load(INFER_DIR / s / "features_fused.npy"),
        fr=np.load(INFER_DIR / s / "features_resnet.npy"),
        dn=np.load(DINO_DIR / s / "features_dinov2.npy"),
    )


def stack_rep(payloads, slides, fn):
    return np.vstack([fn(payloads[s]) for s in slides]).astype(np.float32)


def slide_labels(payloads, slides):
    labels = []
    for s in slides:
        labels.extend([s] * len(payloads[s]["df"]))
    return np.array(labels)


def dominant_lineage(payloads, slides, lineage_map):
    """argmax cell type → lineage group, concatenated across slides."""
    parts = []
    for s in slides:
        pred = payloads[s]["pred"]
        ct_cols = payloads[s]["ct_cols"]
        ct_arr = np.array([lineage_map.get(c, "Unknown") for c in ct_cols])
        parts.append(ct_arr[pred.argmax(axis=1)])
    return np.concatenate(parts)


def knn_indices(X, k):
    """Return [N, k] neighbor indices (excluding self)."""
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", n_jobs=-1)
    nn.fit(X)
    _, idx = nn.kneighbors(X)
    return idx[:, 1:]


def slide_1nn_purity(X, slide_lbl):
    idx = knn_indices(X, k=1).ravel()
    return float((slide_lbl[idx] == slide_lbl).mean())


def knn_jaccard(idx_a, idx_b):
    """Mean Jaccard between row-wise neighbor sets of two index matrices."""
    n, k = idx_a.shape
    jacc = np.zeros(n, dtype=np.float32)
    for i in range(n):
        a = set(idx_a[i])
        b = set(idx_b[i])
        u = len(a | b)
        jacc[i] = 0.0 if u == 0 else len(a & b) / u
    return float(jacc.mean())


def silhouette_by_lineage(X, lineage_arr, rng):
    n = X.shape[0]
    sample = rng.choice(n, size=min(SILHOUETTE_SAMPLE, n), replace=False)
    Xs = X[sample]
    Ls = lineage_arr[sample]
    uniq = np.unique(Ls)
    if len(uniq) < 2:
        return float("nan")
    return float(silhouette_score(Xs, Ls, metric="euclidean", random_state=SEED))


def plot_bars(metrics, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    rep_names = [r[0] for r in REPS]
    rep_dims = {r[0]: r[2] for r in REPS}
    x_labels = [f"{r}\n({rep_dims[r]}-d)" for r in rep_names]
    x = np.arange(len(rep_names))

    # 1) slide 1-NN purity (lower is better)
    vals = [metrics["slide_1nn_purity"][r] for r in rep_names]
    bars = axes[0].bar(x, vals, color=["#4c78a8"] * len(rep_names))
    axes[0].set_xticks(x); axes[0].set_xticklabels(x_labels, fontsize=9)
    axes[0].set_ylabel("slide 1-NN purity"); axes[0].set_ylim(0, 1)
    axes[0].axhline(1 / len(SLIDES), color="grey", linestyle="--", linewidth=0.8,
                    label=f"chance ({1/len(SLIDES):.2f})")
    axes[0].set_title("slide 1-NN purity\n(lower = better batch mix)", fontsize=10)
    axes[0].legend(fontsize=8, loc="upper right")
    for b, v in zip(bars, vals):
        axes[0].text(b.get_x() + b.get_width() / 2, v + 0.01,
                     f"{v:.3f}", ha="center", fontsize=8)

    # 2) kNN overlap (Hist2Cell 3 rep × DINOv2)
    pairs = [r for r in rep_names if r != "features_dinov2"]
    width = 0.35
    xp = np.arange(len(pairs))
    for i, k in enumerate(K_VALUES):
        vals_k = [metrics["knn_overlap"][(r, "features_dinov2", k)] for r in pairs]
        offset = (i - 0.5) * width
        bars = axes[1].bar(xp + offset, vals_k, width=width, label=f"k={k}")
        for b, v in zip(bars, vals_k):
            axes[1].text(b.get_x() + b.get_width() / 2, v + 0.005,
                         f"{v:.3f}", ha="center", fontsize=8)
    axes[1].set_xticks(xp)
    axes[1].set_xticklabels([f"{r}\nvs DINOv2" for r in pairs], fontsize=9)
    axes[1].set_ylabel("kNN overlap (Jaccard)")
    axes[1].set_ylim(0, max(0.25, max(axes[1].get_ylim()) * 1.2))
    axes[1].set_title("Hist2Cell rep vs DINOv2 — kNN neighbor overlap\n(higher = more similar neighbor structure)", fontsize=10)
    axes[1].legend(fontsize=8)

    # 3) silhouette by lineage (higher = better cell-type separation)
    vals = [metrics["silhouette_lineage"][r] for r in rep_names]
    bars = axes[2].bar(x, vals, color=["#54a24b"] * len(rep_names))
    axes[2].set_xticks(x); axes[2].set_xticklabels(x_labels, fontsize=9)
    axes[2].set_ylabel("silhouette by dominant lineage")
    axes[2].set_title("silhouette by dominant lineage\n(higher = stronger cell-type separation)", fontsize=10)
    axes[2].axhline(0, color="grey", linewidth=0.6)
    for b, v in zip(bars, vals):
        axes[2].text(b.get_x() + b.get_width() / 2, v + 0.005,
                     f"{v:.3f}", ha="center", fontsize=8)

    fig.suptitle(
        "Hist2Cell vs DINOv2 — quantitative comparison, cross-slide (3 slides, 15,401 spots)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    print("loading 3 슬라이드...")
    payloads = {s: load_slide(s) for s in SLIDES}
    ct_cols = payloads[SLIDES[0]]["ct_cols"]
    grp = pd.read_csv(GROUPS_CSV)
    lineage_map = dict(zip(grp["cell_type"], grp["group"]))
    lineage_arr = dominant_lineage(payloads, SLIDES, lineage_map)
    slide_lbl = slide_labels(payloads, SLIDES)
    n_total = len(slide_lbl)
    print(f"  total spots: {n_total}")
    print(f"  lineage distribution: {pd.Series(lineage_arr).value_counts().to_dict()}")

    # build cross-slide representations
    X_by_rep = {}
    for name, fn, dim in REPS:
        X_by_rep[name] = stack_rep(payloads, SLIDES, fn)
        print(f"  rep {name}: {X_by_rep[name].shape}")

    # 1) slide 1-NN purity
    print("\n[1/3] slide 1-NN purity...")
    purity = {}
    for name in X_by_rep:
        p = slide_1nn_purity(X_by_rep[name], slide_lbl)
        purity[name] = p
        print(f"  {name:18s}  {p:.4f}")

    # 2) kNN overlap (Hist2Cell 3 rep × DINOv2)
    print("\n[2/3] kNN overlap (Jaccard) Hist2Cell × DINOv2...")
    dino_idx = {k: knn_indices(X_by_rep["features_dinov2"], k) for k in K_VALUES}
    overlap = {}
    for name in X_by_rep:
        if name == "features_dinov2":
            continue
        for k in K_VALUES:
            a = knn_indices(X_by_rep[name], k)
            j = knn_jaccard(a, dino_idx[k])
            overlap[(name, "features_dinov2", k)] = j
            print(f"  k={k:2d}  {name:18s} vs DINOv2  {j:.4f}")

    # 3) silhouette by lineage
    print(f"\n[3/3] silhouette by lineage (sample {SILHOUETTE_SAMPLE} of {n_total})...")
    sil = {}
    for name in X_by_rep:
        s = silhouette_by_lineage(X_by_rep[name], lineage_arr, rng)
        sil[name] = s
        print(f"  {name:18s}  {s:.4f}")

    # write CSV
    rows = []
    for name in X_by_rep:
        rows.append({"representation": name, "metric": "slide_1nn_purity",
                     "k": None, "vs": None, "value": purity[name]})
    for (a, b, k), v in overlap.items():
        rows.append({"representation": a, "metric": "knn_overlap_jaccard",
                     "k": k, "vs": b, "value": v})
    for name in X_by_rep:
        rows.append({"representation": name, "metric": "silhouette_lineage",
                     "k": None, "vs": None, "value": sil[name]})
    df = pd.DataFrame(rows)
    csv_path = OUT_DIR / "metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nsaved: {csv_path}")

    # plot
    metrics = {
        "slide_1nn_purity": purity,
        "knn_overlap": overlap,
        "silhouette_lineage": sil,
    }
    png_path = OUT_DIR / "metrics_bars.png"
    plot_bars(metrics, png_path)
    print(f"saved: {png_path}")


if __name__ == "__main__":
    main()
