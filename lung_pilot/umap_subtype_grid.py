"""3×4 UMAP grid (3 슬라이드 × 4 representation) — 3 가지 색칠.

UMAP embedding 은 1회 fit + 캐시. umap_compare.py 와 같은 파라미터
(random_state=42, n_jobs=1 deterministic) 라 같은 PNG 와 동일 좌표 보장.
캐시 위치: lung_pilot/umap_output/embeddings/<slide>_<rep>_umap2d.npy.

3 출력 PNG (lung_pilot/umap_output/):
  per_slide_grid_lineage.png      — color = dominant cell-type lineage (10 groups; only present color)
  per_slide_grid_epithelial.png   — color = Epithelial subtype (17 cell types; non-Epi grey)
  per_slide_grid_stromal.png      — color = Stromal subtype (16 cell types; non-Stro grey)

Subtype 색칠 정책: 각 spot 의 dominant cell type (argmax) 의 lineage 가
대상 카테고리면 해당 sub-type 색, 아니면 옅은 grey. 즉 *해당 lineage 의
spot 분포 + 그 안 sub-type 다양성* 을 본다.

Usage:
    .venv/bin/python lung_pilot/umap_subtype_grid.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap

ROOT = Path(__file__).resolve().parent.parent
INFER_DIR = ROOT / "lung_pilot" / "inference_output"
DINO_DIR = ROOT / "lung_pilot" / "dino_output"
OUT_DIR = ROOT / "lung_pilot" / "umap_output"
EMB_DIR = OUT_DIR / "embeddings"
GROUPS_CSV = ROOT / "inference" / "analysis_spatial" / "cell_type_groups.csv"

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

UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)
GREY = (0.85, 0.85, 0.85, 0.4)


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


def short_slide(s):
    return s.replace("TCGA-05-", "").replace("-01A-01-", "-")


def get_embedding(slide, rep_name, X):
    cache_path = EMB_DIR / f"{slide}_{rep_name}_umap2d.npy"
    if cache_path.exists():
        return np.load(cache_path)
    print(f"  fitting UMAP for {slide} × {rep_name} ({X.shape[1]}-d)...")
    z = umap.UMAP(**UMAP_KW).fit_transform(X)
    np.save(cache_path, z)
    return z


def plot_grid(emb_cache, payloads, color_fn, title, out_path,
              legend_handles, legend_title, legend_ncol=1):
    n_row, n_col = len(SLIDES), len(REPS)
    fig, axes = plt.subplots(n_row, n_col, figsize=(5.5 * n_col, 5.0 * n_row))
    for i, s in enumerate(SLIDES):
        colors = color_fn(payloads[s])
        for j, (rep_name, _, dim) in enumerate(REPS):
            ax = axes[i, j]
            z = emb_cache[(s, rep_name)]
            ax.scatter(z[:, 0], z[:, 1], c=colors, s=4, alpha=0.65, linewidths=0)
            if i == 0:
                ax.set_title(f"{rep_name}  ({dim}-d)", fontsize=11)
            if j == 0:
                ax.set_ylabel(
                    f"{short_slide(s)}\n({len(payloads[s]['df'])} spots)",
                    fontsize=10,
                )
            ax.set_xticks([]); ax.set_yticks([])
    fig.legend(handles=legend_handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=8.5,
               title=legend_title, ncol=legend_ncol)
    fig.suptitle(title, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 0.83 if legend_ncol == 1 else 0.78, 0.97])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    EMB_DIR.mkdir(parents=True, exist_ok=True)
    grp = pd.read_csv(GROUPS_CSV)

    payloads = {s: load_slide(s) for s in SLIDES}
    ct_cols = payloads[SLIDES[0]]["ct_cols"]
    lineage_map = dict(zip(grp["cell_type"], grp["group"]))

    # 1. UMAP embedding cache (fit if missing)
    print("caching UMAP embeddings (re-uses .npy if present)...")
    emb_cache = {}
    for s in SLIDES:
        for rep_name, fn, _ in REPS:
            X = fn(payloads[s])
            emb_cache[(s, rep_name)] = get_embedding(s, rep_name, X)
    print("  cache ready.")

    # 2. Figure A — dominant lineage (10 groups in legend)
    all_lineages = sorted(set(grp["group"]))
    lin_cmap = plt.get_cmap("tab10", max(10, len(all_lineages)))
    lin_color_map = {l: lin_cmap(i % lin_cmap.N) for i, l in enumerate(all_lineages)}

    lineage_arr = np.array([lineage_map[c] for c in ct_cols])

    def lineage_color_fn(payload):
        dom_idx = payload["pred"].argmax(axis=1)
        dom_lineage = lineage_arr[dom_idx]
        return np.array([lin_color_map[l] for l in dom_lineage])

    legend_lin = [
        plt.Line2D([], [], marker="o", linestyle="", color=lin_color_map[l],
                   markersize=8, label=l)
        for l in all_lineages
    ]
    plot_grid(
        emb_cache, payloads, lineage_color_fn,
        "Per-slide UMAP grid (3 slides × 4 representations) — color = dominant cell-type lineage (10 groups; only spots with dominant ∈ group get colored)",
        OUT_DIR / "per_slide_grid_lineage.png",
        legend_lin, "lineage",
    )
    print(f"saved: {OUT_DIR / 'per_slide_grid_lineage.png'}")

    # 3. Figure B — Epithelial subtype (airway 14 + alveolar 3 = 17)
    epi_types = sorted(grp.loc[grp["group"].str.startswith("Epithelial"), "cell_type"].tolist())
    epi_cmap = plt.get_cmap("tab20", max(20, len(epi_types)))
    epi_color_map = {t: epi_cmap(i % epi_cmap.N) for i, t in enumerate(epi_types)}

    def epi_color_fn(payload):
        dom_idx = payload["pred"].argmax(axis=1)
        dom_lineage = lineage_arr[dom_idx]
        dom_type = np.array(ct_cols)[dom_idx]
        out = np.zeros((len(dom_idx), 4), dtype=float)
        for i, (t, l) in enumerate(zip(dom_type, dom_lineage)):
            out[i] = epi_color_map[t] if l.startswith("Epithelial") else GREY
        return out

    legend_epi = [
        plt.Line2D([], [], marker="o", linestyle="", color=epi_color_map[t],
                   markersize=8, label=t)
        for t in epi_types
    ]
    legend_epi.append(
        plt.Line2D([], [], marker="o", linestyle="", color=GREY,
                   markersize=8, label="non-Epithelial dominant")
    )
    plot_grid(
        emb_cache, payloads, epi_color_fn,
        "Per-slide UMAP grid — color = dominant Epithelial cell type (17 sub-types; airway 14 + alveolar 3). non-Epi spots in grey.",
        OUT_DIR / "per_slide_grid_epithelial.png",
        legend_epi, "Epithelial cell type",
    )
    print(f"saved: {OUT_DIR / 'per_slide_grid_epithelial.png'}")

    # 4. Figure C — Stromal subtype (fibroblast 6 + muscle 6 + other 4 = 16)
    stro_types = sorted(grp.loc[grp["group"].str.startswith("Stromal"), "cell_type"].tolist())
    stro_cmap = plt.get_cmap("tab20", max(20, len(stro_types)))
    stro_color_map = {t: stro_cmap(i % stro_cmap.N) for i, t in enumerate(stro_types)}

    def stro_color_fn(payload):
        dom_idx = payload["pred"].argmax(axis=1)
        dom_lineage = lineage_arr[dom_idx]
        dom_type = np.array(ct_cols)[dom_idx]
        out = np.zeros((len(dom_idx), 4), dtype=float)
        for i, (t, l) in enumerate(zip(dom_type, dom_lineage)):
            out[i] = stro_color_map[t] if l.startswith("Stromal") else GREY
        return out

    legend_stro = [
        plt.Line2D([], [], marker="o", linestyle="", color=stro_color_map[t],
                   markersize=8, label=t)
        for t in stro_types
    ]
    legend_stro.append(
        plt.Line2D([], [], marker="o", linestyle="", color=GREY,
                   markersize=8, label="non-Stromal dominant")
    )
    plot_grid(
        emb_cache, payloads, stro_color_fn,
        "Per-slide UMAP grid — color = dominant Stromal cell type (16 sub-types; fibroblast 6 + muscle 6 + other 4). non-Stro spots in grey.",
        OUT_DIR / "per_slide_grid_stromal.png",
        legend_stro, "Stromal cell type",
    )
    print(f"saved: {OUT_DIR / 'per_slide_grid_stromal.png'}")

    print("done.")


if __name__ == "__main__":
    main()
