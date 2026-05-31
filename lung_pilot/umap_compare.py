"""UMAP 비교 — Hist2Cell 3 rep + DINOv2 ViT-B/14, 3 슬라이드.

산출 (lung_pilot/umap_output/):
  per_slide_<slide>.png      — 한 슬라이드의 4 rep UMAP, color = dominant cell-type lineage
  cross_slide_combined.png   — 3 슬라이드 합쳐 4 rep UMAP, color = slide (batch effect 확인)

Representations (4):
  prediction_log1p — predictions.npy [N,80] 의 log1p (cell2location abundance, row_sum 1.3~62.8)
  features_fused   — features_fused.npy [N,256] (Hist2Cell 의 fused_head 직전, visual+GAT+TF)
  features_resnet  — features_resnet.npy [N,512] (Hist2Cell 의 ResNet18 backbone, graph 미포함)
  features_dinov2  — features_dinov2.npy [N,768] (DINOv2 ViT-B/14 CLS, 외부 self-supervised)

Usage:
    .venv/bin/python lung_pilot/umap_compare.py
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
GROUPS_CSV = ROOT / "inference" / "analysis_spatial" / "cell_type_groups.csv"

SLIDES = [
    "TCGA-05-4245-01A-01-BS1",
    "TCGA-05-4245-01A-01-TS1",
    "TCGA-05-4390-01A-01-BS1",
]

REPS = [
    ("prediction_log1p", lambda d: np.log1p(d["pred"]),  80),
    ("features_fused",   lambda d: d["ff"],              256),
    ("features_resnet",  lambda d: d["fr"],              512),
    ("features_dinov2",  lambda d: d["dn"],              768),
]

UMAP_KW = dict(n_neighbors=15, min_dist=0.1, metric="euclidean", random_state=42)


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


def load_lineage_map(ct_cols):
    grp = pd.read_csv(GROUPS_CSV)
    m = dict(zip(grp["cell_type"], grp["group"]))
    return np.array([m.get(c, "Unknown") for c in ct_cols])


def categorical_colors(labels):
    uniq = sorted(set(labels))
    cmap = plt.get_cmap("tab20", max(20, len(uniq)))
    return {u: cmap(i % cmap.N) for i, u in enumerate(uniq)}, uniq


def plot_per_slide(slide, payload, lineage_arr, out_path):
    dom_idx = payload["pred"].argmax(axis=1)
    dom_lineage = lineage_arr[dom_idx]
    color_map, uniq = categorical_colors(dom_lineage)
    point_colors = [color_map[l] for l in dom_lineage]

    fig, axes = plt.subplots(1, len(REPS), figsize=(5.5 * len(REPS), 6.5))
    for ax, (rep_name, fn, dim) in zip(axes, REPS):
        X = fn(payload)
        z = umap.UMAP(**UMAP_KW).fit_transform(X)
        ax.scatter(z[:, 0], z[:, 1], c=point_colors, s=4, alpha=0.6, linewidths=0)
        ax.set_title(f"{rep_name}  ({dim}-d)")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.set_xticks([]); ax.set_yticks([])

    handles = [plt.Line2D([], [], marker="o", linestyle="",
                          color=color_map[u], markersize=8, label=u) for u in uniq]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5),
               fontsize=9, title="dominant lineage")
    fig.suptitle(
        f"{slide}  —  {len(payload['df'])} spots, color = dominant cell-type lineage (argmax → group)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 0.96])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_cross_slide(payloads, slides, out_path):
    slide_colors = dict(zip(slides, plt.get_cmap("tab10").colors[: len(slides)]))
    fig, axes = plt.subplots(1, len(REPS), figsize=(5.5 * len(REPS), 6.5))
    for ax, (rep_name, fn, dim) in zip(axes, REPS):
        Xs, slide_lbl = [], []
        for s in slides:
            X = fn(payloads[s])
            Xs.append(X)
            slide_lbl.extend([s] * len(X))
        X_all = np.vstack(Xs)
        slide_lbl = np.array(slide_lbl)
        z = umap.UMAP(**UMAP_KW).fit_transform(X_all)
        for s in slides:
            mask = slide_lbl == s
            ax.scatter(z[mask, 0], z[mask, 1], c=[slide_colors[s]],
                       s=3, alpha=0.45, linewidths=0,
                       label=s if ax is axes[0] else None)
        ax.set_title(f"{rep_name}  ({dim}-d, combined)")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.set_xticks([]); ax.set_yticks([])
    handles = [plt.Line2D([], [], marker="o", linestyle="",
                          color=slide_colors[s], markersize=8, label=s) for s in slides]
    fig.legend(handles=handles, loc="center right",
               bbox_to_anchor=(1.0, 0.5), fontsize=9, title="slide")
    total = sum(len(payloads[s]["df"]) for s in slides)
    fig.suptitle(
        f"Cross-slide combined UMAP  —  3 slides, {total} spots total, color = slide   "
        f"(slides separate = batch effect; slides mixed = tissue-generic representation)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 0.96])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    global INFER_DIR, DINO_DIR, OUT_DIR
    import argparse
    ap = argparse.ArgumentParser(description="UMAP 4-rep 비교 (경로 미지정 시 224 기본)")
    ap.add_argument("--infer-dir", default=str(INFER_DIR))
    ap.add_argument("--dino-dir", default=str(DINO_DIR))
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    a = ap.parse_args()
    INFER_DIR, DINO_DIR, OUT_DIR = Path(a.infer_dir), Path(a.dino_dir), Path(a.out_dir)
    print(f"INFER={INFER_DIR}\nDINO={DINO_DIR}\nOUT={OUT_DIR}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payloads = {s: load_slide(s) for s in SLIDES}
    ct_cols = payloads[SLIDES[0]]["ct_cols"]
    lineage_arr = load_lineage_map(ct_cols)

    for s in SLIDES:
        out = OUT_DIR / f"per_slide_{s}.png"
        print(f"[per-slide] {s} ({len(payloads[s]['df'])} spots) → {out}")
        plot_per_slide(s, payloads[s], lineage_arr, out)

    out = OUT_DIR / "cross_slide_combined.png"
    print(f"[cross-slide] → {out}")
    plot_cross_slide(payloads, SLIDES, out)

    print("done.")


if __name__ == "__main__":
    main()
