"""Joint CCA — 두 슬라이드를 한 모델로 학습.

지금까지의 모든 per-slide CCA → 부호 정렬 → 합본 접근은 두 슬라이드의
axis 의미가 달라지면 합본 통계가 ad hoc 정렬에 기댔다.  본 빌더는 그
구조를 *대체* 하는 본 분석.

파이프라인:
 1. 두 슬라이드의 raw H (80 cell type), P (4168 공통 gene) 로드
 2. 슬라이드별 z-score (= batch correction) — 평균/스케일 차이 제거
 3. 두 행렬 vstack → H_joint (94, 80), P_joint (94, 4168)
 4. modality 별 PCA10 + 단일 CCA (3 component)
 5. permutation null — proteomics ROI 순서 shuffle, 1000회
 6. axis × {risk, section group, slide} 정렬 검정
 7. loadings 추출 — 두 슬라이드 공통 모듈

산출:
 joint_cca_scores.csv   (94 ROI × {slide, section, group, risk,
                                    H_c1, P_c1, H_c2, P_c2, H_c3, P_c3})
 joint_cca_summary.csv  (canonical r + permutation null + slide·risk r)
 joint_cca_loadings.csv (axis × cell_type/gene × loading)
 joint_cca_showcase.png (4 panel)
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from _proof_ver2_lib import (
    SlideConfig, build_roi_signatures, load_proteomics_matrix,
    align_modalities, N_PCS, N_CCA_COMP, N_PERM, RANDOM_SEED,
)

SECTION_GROUP = {
    "a":"Tumor-h","b":"Tumor-l","c":"Tcell-h","d":"Tcell-l","t":"Tumor-ctrl",
    "e":"Tumor-h","f":"Tumor-l","g":"Tcell-h","h":"Tcell-l","v":"Tumor-ctrl",
}
GROUP_COLOR = {
    "Tumor-h":"#d62728","Tumor-l":"#1f77b4","Tumor-ctrl":"#7f7f7f",
    "Tcell-h":"#2ca02c","Tcell-l":"#9467bd",
}
GROUP_ORDER = ["Tumor-h","Tumor-l","Tumor-ctrl","Tcell-h","Tcell-l"]
SLIDE_MARKER = {"slide1":"o","slide2":"^"}
TUMOR_GROUPS = {"Tumor-h","Tumor-l","Tumor-ctrl"}


def slide_cfg(slide_id: str) -> SlideConfig:
    if slide_id == "slide1":
        return SlideConfig(
            name="slide1", pred_csv=Path("/home/sjhong/hist2cell/inference/slide1_085_12_v2/predictions.csv"),
            roi_pkl=ROOT/"1_085_12"/"1_085_12_ROI_groups.pkl",
            npy=ROOT/"1_085_12"/"meteo_1_085_12_coords.npy",
            section_label={}, section_color={}, sample_section_prefixes="abcdt", out_dir=HERE,
        )
    return SlideConfig(
        name="slide2", pred_csv=Path("/home/sjhong/hist2cell/inference/slide2_152_19_v2/predictions.csv"),
        roi_pkl=ROOT/"1_152_19"/"1_152_19_ROI_groups.pkl",
        npy=ROOT/"1_152_19"/"meteo_1_152_19_coords.npy",
        section_label={}, section_color={}, sample_section_prefixes="efghv", out_dir=HERE,
    )


def load_risk(slide_id: str) -> dict:
    folder = "1_085_12" if slide_id == "slide1" else "1_152_19"
    with open(ROOT/folder/f"{folder}_ROI_groups_risk_scores.pkl", "rb") as f:
        return pickle.load(f)


def build_joint_matrices():
    """두 슬라이드의 raw H, P 를 로드하고 공통 gene 만 남긴 후 슬라이드별
    z-score 정규화 한 다음 vertical stack."""
    out = {}
    raw_H = []
    raw_P = []
    meta_rows = []
    for sid in ["slide1", "slide2"]:
        cfg = slide_cfg(sid)
        sig_df, cell_cols = build_roi_signatures(cfg)
        log2_f, slide_cols = load_proteomics_matrix(cfg)
        common, H, P, sig_aligned, genes = align_modalities(
            sig_df, log2_f, slide_cols, cell_cols)
        out[sid] = {"common": common, "H": H, "P": P, "genes": genes,
                    "cell_cols": cell_cols, "sections": [t[0] for t in common]}
    # gene 교집합
    g1 = set(out["slide1"]["genes"]); g2 = set(out["slide2"]["genes"])
    common_genes = sorted(g1 & g2)
    print(f"[*] gene intersection = {len(common_genes)} "
          f"(slide1 raw {len(g1)}, slide2 raw {len(g2)})")
    # P 행렬을 공통 gene 으로 잘라내기
    for sid in ["slide1", "slide2"]:
        gene_idx = [out[sid]["genes"].index(g) for g in common_genes]
        out[sid]["P_common"] = out[sid]["P"][:, gene_idx]
    # 슬라이드별 z-score (batch correction)
    H1z = StandardScaler().fit_transform(out["slide1"]["H"])
    H2z = StandardScaler().fit_transform(out["slide2"]["H"])
    P1z = StandardScaler().fit_transform(out["slide1"]["P_common"])
    P2z = StandardScaler().fit_transform(out["slide2"]["P_common"])
    H_joint = np.vstack([H1z, H2z])
    P_joint = np.vstack([P1z, P2z])
    slides = (["slide1"] * len(out["slide1"]["common"]) +
              ["slide2"] * len(out["slide2"]["common"]))
    tubes  = out["slide1"]["common"] + out["slide2"]["common"]
    sections = out["slide1"]["sections"] + out["slide2"]["sections"]
    cell_cols = out["slide1"]["cell_cols"]
    return H_joint, P_joint, slides, tubes, sections, cell_cols, common_genes


def joint_cca(H, P, n_pcs=N_PCS, n_components=N_CCA_COMP, seed=RANDOM_SEED):
    n_pcs = min(n_pcs, H.shape[0]-1, H.shape[1], P.shape[1])
    pca_h = PCA(n_components=n_pcs, random_state=seed)
    pca_p = PCA(n_components=n_pcs, random_state=seed)
    H_pc = pca_h.fit_transform(H)
    P_pc = pca_p.fit_transform(P)
    cca = CCA(n_components=n_components, max_iter=1000)
    cca.fit(H_pc, P_pc)
    Hc, Pc = cca.transform(H_pc, P_pc)
    rs = [float(pearsonr(Hc[:, i], Pc[:, i])[0]) for i in range(n_components)]
    h_load = pca_h.components_.T @ cca.x_loadings_
    p_load = pca_p.components_.T @ cca.y_loadings_
    return dict(rs=rs, Hc=Hc, Pc=Pc, h_load=h_load, p_load=p_load,
                pca_h_var=pca_h.explained_variance_ratio_,
                pca_p_var=pca_p.explained_variance_ratio_,
                n_pcs=n_pcs)


def perm_null_top_r(H, P, n_perm=N_PERM, n_pcs=N_PCS, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    n = H.shape[0]
    null = np.zeros(n_perm)
    n_pcs = min(n_pcs, n-1, H.shape[1], P.shape[1])
    for k in range(n_perm):
        perm = rng.permutation(n)
        pca_h = PCA(n_components=n_pcs, random_state=seed).fit_transform(H)
        pca_p = PCA(n_components=n_pcs, random_state=seed).fit_transform(P[perm])
        cca = CCA(n_components=1, max_iter=500).fit(pca_h, pca_p)
        Hc1, Pc1 = cca.transform(pca_h, pca_p)
        null[k] = pearsonr(Hc1[:, 0], Pc1[:, 0])[0]
    return null


def perm_p_pearson(x, y, observed_r, n_perm=1000, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm)
    for k in range(n_perm):
        ys = rng.permutation(y)
        null[k] = pearsonr(x, ys)[0]
    return float(np.mean(np.abs(null) >= abs(observed_r)))


def main():
    print("[*] building joint matrices …")
    H, P, slides, tubes, sections, cell_cols, genes = build_joint_matrices()
    print(f"    joint H = {H.shape}, joint P = {P.shape}, n_ROI={H.shape[0]}")

    print("[*] joint CCA …")
    res = joint_cca(H, P)
    print(f"    canonical r (3 axes) = {[f'{r:+.3f}' for r in res['rs']]}")
    print(f"    n_pcs = {res['n_pcs']}, PCA H var3 = {res['pca_h_var'][:3].sum():.2%}, "
          f"PCA P var3 = {res['pca_p_var'][:3].sum():.2%}")

    print("[*] permutation null (1000) …")
    null = perm_null_top_r(H, P)
    obs_top = res["rs"][0]
    p_top = float(np.mean(np.abs(null) >= abs(obs_top)))
    null_mean = float(np.mean(null)); null_95_lo, null_95_hi = np.percentile(null, [2.5, 97.5])
    print(f"    observed top r = {obs_top:+.3f},  null mean = {null_mean:+.3f},  "
          f"95% [{null_95_lo:+.3f}, {null_95_hi:+.3f}],  p = {p_top}")

    # risk scores
    risk1 = load_risk("slide1"); risk2 = load_risk("slide2")
    risk = np.array([risk1[t] if s == "slide1" else risk2[t]
                     for t, s in zip(tubes, slides)])

    # axis vs risk + slide effect ───────────────────────────
    groups = [SECTION_GROUP[s] for s in sections]
    df = pd.DataFrame({
        "slide": slides, "tube_id": tubes, "section": sections,
        "group": groups, "risk": risk,
        "H_c1": res["Hc"][:, 0], "P_c1": res["Pc"][:, 0],
        "H_c2": res["Hc"][:, 1], "P_c2": res["Pc"][:, 1],
        "H_c3": res["Hc"][:, 2], "P_c3": res["Pc"][:, 2],
    })
    df.to_csv(HERE / "joint_cca_scores.csv", index=False)

    # risk correlation per axis × modality × {all, Tumor-only, slide1, slide2}
    rows = []
    is_tumor = df["group"].isin(TUMOR_GROUPS)
    for axis in [1, 2, 3]:
        for mod_label, col in [("Hist2Cell", f"H_c{axis}"),
                                ("Proteomics", f"P_c{axis}")]:
            for subset_label, mask in [
                ("all_ROI",     pd.Series([True] * len(df))),
                ("Tumor_only",  is_tumor),
                ("slide1_only", df["slide"] == "slide1"),
                ("slide2_only", df["slide"] == "slide2"),
                ("slide1_Tumor", (df["slide"] == "slide1") & is_tumor),
                ("slide2_Tumor", (df["slide"] == "slide2") & is_tumor),
            ]:
                sub = df[mask]
                if len(sub) < 5: continue
                x = sub[col].values; y = sub["risk"].values
                r_p, _ = pearsonr(x, y)
                r_s, _ = spearmanr(x, y)
                p_perm = perm_p_pearson(x, y, r_p)
                rows.append({"axis": axis, "modality": mod_label,
                             "subset": subset_label, "n": len(sub),
                             "pearson_r": r_p, "p_perm": p_perm,
                             "spearman_rho": r_s,
                             "r_squared_pct": 100*r_p**2})
    corr = pd.DataFrame(rows)
    corr.to_csv(HERE / "joint_cca_risk_correlations.csv", index=False)
    print("[*] joint CCA risk correlations (top 10 by |r|):")
    print(corr.assign(abs_r=corr["pearson_r"].abs())
                .sort_values("abs_r", ascending=False)
                .head(10)[["axis","modality","subset","n","pearson_r",
                            "r_squared_pct","p_perm","spearman_rho"]]
                .to_string(index=False))

    # slide effect on each axis — t-test 대신 mean per slide
    slide_eff = []
    for axis in [1, 2, 3]:
        for col in [f"H_c{axis}", f"P_c{axis}"]:
            m1 = df[df["slide"]=="slide1"][col].mean()
            m2 = df[df["slide"]=="slide2"][col].mean()
            slide_eff.append({"axis": axis, "modality_col": col,
                              "slide1_mean": m1, "slide2_mean": m2,
                              "abs_diff": abs(m1-m2)})
    slide_eff_df = pd.DataFrame(slide_eff)
    print("\n[*] slide-effect on each axis (joint CCA score means):")
    print(slide_eff_df.to_string(index=False))

    summary = {
        "n_ROI_joint": H.shape[0],
        "n_genes_common": P.shape[1],
        "canonical_r_axis1": res["rs"][0],
        "canonical_r_axis2": res["rs"][1],
        "canonical_r_axis3": res["rs"][2],
        "null_mean_top": null_mean,
        "null_95_lo_top": null_95_lo,
        "null_95_hi_top": null_95_hi,
        "perm_p_top": p_top,
    }
    pd.DataFrame([summary]).to_csv(HERE / "joint_cca_summary.csv", index=False)

    # loadings
    h_load_df = pd.DataFrame(res["h_load"], index=cell_cols,
                              columns=[f"axis{i+1}" for i in range(3)])
    h_load_df.insert(0, "feature_type", "cell_type")
    p_load_df = pd.DataFrame(res["p_load"], index=genes,
                              columns=[f"axis{i+1}" for i in range(3)])
    p_load_df.insert(0, "feature_type", "gene")
    pd.concat([h_load_df, p_load_df]).to_csv(HERE / "joint_cca_loadings.csv")

    # ── figure ───────────────────────────────────────────────────────
    print("[*] building 4-panel joint CCA showcase …")
    fig = plt.figure(figsize=(20, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.30,
                           left=0.055, right=0.985, top=0.93, bottom=0.07)
    axA = fig.add_subplot(gs[0, 0])  # axis 1 paired (H vs P)
    axB = fig.add_subplot(gs[0, 1])  # axis 1 vs risk
    axC = fig.add_subplot(gs[0, 2])  # axis 2 vs risk
    axD = fig.add_subplot(gs[1, 0])  # slide effect — H_c1, P_c1 per slide
    axE = fig.add_subplot(gs[1, 1])  # risk vs axis 2 with slide split
    axF = fig.add_subplot(gs[1, 2])  # axis 1/2/3 vs risk grid heatmap

    # A — axis 1 paired
    for sid in ["slide1","slide2"]:
        for g in GROUP_ORDER:
            sub = df[(df.slide==sid) & (df.group==g)]
            if sub.empty: continue
            axA.scatter(sub["H_c1"], sub["P_c1"],
                         marker=SLIDE_MARKER[sid], c=GROUP_COLOR[g],
                         s=72 if sid=="slide1" else 64,
                         edgecolors="black", linewidths=0.6, alpha=0.85,
                         label=f"{g} ({sid})")
    lo = min(df["H_c1"].min(), df["P_c1"].min())-0.3
    hi = max(df["H_c1"].max(), df["P_c1"].max())+0.3
    axA.plot([lo,hi],[lo,hi],"k--",lw=1,alpha=0.5)
    axA.set_xlim(lo,hi); axA.set_ylim(lo,hi)
    axA.set_xlabel("Hist2Cell joint axis-1 score")
    axA.set_ylabel("Proteomics joint axis-1 score")
    axA.set_title(f"A) Joint CCA axis 1 paired  (canonical r = {res['rs'][0]:+.3f})\n"
                   f"null mean {null_mean:+.3f}, 95% [{null_95_lo:+.3f}, {null_95_hi:+.3f}], "
                   f"p={p_top:.3f}",
                   fontsize=10)
    axA.grid(alpha=0.2)
    axA.legend(loc="upper left", fontsize=6, frameon=False, ncol=2)

    # B — axis 1 score (Hist2Cell) vs risk
    rb_h_all = corr[(corr.axis==1)&(corr.modality=="Hist2Cell")&(corr.subset=="all_ROI")].iloc[0]
    rb_h_tu  = corr[(corr.axis==1)&(corr.modality=="Hist2Cell")&(corr.subset=="Tumor_only")].iloc[0]
    for sid in ["slide1","slide2"]:
        for g in GROUP_ORDER:
            sub = df[(df.slide==sid) & (df.group==g)]
            if sub.empty: continue
            axB.scatter(sub["risk"], sub["H_c1"],
                         marker=SLIDE_MARKER[sid], c=GROUP_COLOR[g],
                         s=72 if sid=="slide1" else 64,
                         edgecolors="black", linewidths=0.6, alpha=0.85)
    coef = np.polyfit(df["risk"], df["H_c1"], 1)
    xs = np.linspace(df["risk"].min(), df["risk"].max(), 100)
    axB.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
    axB.set_xlabel("ROI risk score")
    axB.set_ylabel("Hist2Cell joint axis-1 score")
    axB.set_title(f"B) joint axis-1 (H2C) vs risk\n"
                   f"all r={rb_h_all.pearson_r:+.3f} (r²={rb_h_all.r_squared_pct:.0f}%, p={rb_h_all.p_perm:.3f})  |  "
                   f"Tumor r={rb_h_tu.pearson_r:+.3f} (r²={rb_h_tu.r_squared_pct:.0f}%, p={rb_h_tu.p_perm:.3f})",
                   fontsize=10)
    axB.grid(alpha=0.2)

    # C — axis 2 score (Hist2Cell) vs risk
    rc_h_all = corr[(corr.axis==2)&(corr.modality=="Hist2Cell")&(corr.subset=="all_ROI")].iloc[0]
    rc_h_tu  = corr[(corr.axis==2)&(corr.modality=="Hist2Cell")&(corr.subset=="Tumor_only")].iloc[0]
    for sid in ["slide1","slide2"]:
        for g in GROUP_ORDER:
            sub = df[(df.slide==sid) & (df.group==g)]
            if sub.empty: continue
            axC.scatter(sub["risk"], sub["H_c2"],
                         marker=SLIDE_MARKER[sid], c=GROUP_COLOR[g],
                         s=72 if sid=="slide1" else 64,
                         edgecolors="black", linewidths=0.6, alpha=0.85)
    coef = np.polyfit(df["risk"], df["H_c2"], 1)
    xs = np.linspace(df["risk"].min(), df["risk"].max(), 100)
    axC.plot(xs, np.polyval(coef, xs), "k--", lw=1.2, alpha=0.6)
    axC.set_xlabel("ROI risk score")
    axC.set_ylabel("Hist2Cell joint axis-2 score")
    axC.set_title(f"C) joint axis-2 (H2C) vs risk\n"
                   f"all r={rc_h_all.pearson_r:+.3f} (r²={rc_h_all.r_squared_pct:.0f}%, p={rc_h_all.p_perm:.3f})  |  "
                   f"Tumor r={rc_h_tu.pearson_r:+.3f} (r²={rc_h_tu.r_squared_pct:.0f}%, p={rc_h_tu.p_perm:.3f})",
                   fontsize=10)
    axC.grid(alpha=0.2)

    # D — slide-effect boxplot per axis
    box_data = []; box_labels = []; box_colors=[]
    for axis in [1,2,3]:
        for col in [f"H_c{axis}", f"P_c{axis}"]:
            for sid, c in [("slide1","#444444"), ("slide2","#aaaaaa")]:
                box_data.append(df[df.slide==sid][col].values)
                short = col.replace("_c","_a")
                box_labels.append(f"{short}\n{sid}")
                box_colors.append(c)
    bp = axD.boxplot(box_data, labels=box_labels, patch_artist=True, widths=0.55, showfliers=False)
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c); patch.set_alpha(0.6); patch.set_edgecolor("black")
    for med in bp["medians"]: med.set_color("black"); med.set_linewidth(1.4)
    axD.axhline(0, c="k", lw=0.5, alpha=0.3)
    axD.tick_params(axis="x", rotation=45, labelsize=7)
    axD.set_ylabel("joint axis score")
    axD.set_title("D) Slide effect  -  axis score distribution per slide\n"
                   "(large slide gap = axis captures batch identity, not biology)",
                   fontsize=10)
    axD.grid(axis="y", alpha=0.2)

    # E — axis 2 vs risk with per-slide regression lines
    for sid, m_clr in [("slide1", "#1f77b4"), ("slide2", "#d62728")]:
        sub = df[(df.slide==sid)]
        axE.scatter(sub["risk"], sub["H_c2"], marker=SLIDE_MARKER[sid],
                     c=[GROUP_COLOR[g] for g in sub["group"]],
                     s=72 if sid=="slide1" else 64,
                     edgecolors="black", linewidths=0.6, alpha=0.85)
        coef = np.polyfit(sub["risk"], sub["H_c2"], 1)
        xs = np.linspace(sub["risk"].min(), sub["risk"].max(), 100)
        r_p, _ = pearsonr(sub["risk"], sub["H_c2"])
        axE.plot(xs, np.polyval(coef, xs), c=m_clr, lw=1.5, alpha=0.7,
                  label=f"{sid} r={r_p:+.3f} (r²={100*r_p**2:.0f}%)")
    axE.set_xlabel("ROI risk score")
    axE.set_ylabel("Hist2Cell joint axis-2 score")
    axE.set_title("E) joint axis-2 vs risk  -  per-slide regression\n"
                   "(do both slides align same direction on the *same* axis?)",
                   fontsize=10)
    axE.legend(loc="upper left", fontsize=8, frameon=False)
    axE.grid(alpha=0.2)

    # F — risk correlation heatmap, joint CCA edition
    pivot = corr[corr["subset"].isin(["all_ROI", "Tumor_only", "slide1_only",
                                       "slide2_only", "slide1_Tumor", "slide2_Tumor"])]\
              .pivot_table(index=["axis","modality"], columns="subset",
                           values="pearson_r")
    cols_order = ["all_ROI", "Tumor_only", "slide1_only", "slide1_Tumor",
                  "slide2_only", "slide2_Tumor"]
    pivot = pivot[cols_order]
    data = pivot.values
    im = axF.imshow(data, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    axF.set_xticks(range(pivot.shape[1]))
    axF.set_xticklabels(pivot.columns.tolist(), fontsize=7, rotation=30, ha="right")
    axF.set_yticks(range(pivot.shape[0]))
    axF.set_yticklabels([f"axis {a} {m}" for a,m in pivot.index.tolist()], fontsize=8)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i,j]
            txt_c = "white" if abs(val)>0.35 else "black"
            axF.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8, color=txt_c)
    fig.colorbar(im, ax=axF, shrink=0.85, label="Pearson r")
    axF.set_title("F) joint CCA  axis × modality × subset  vs  risk Pearson r",
                   fontsize=10)

    fig.suptitle(
        "Joint CCA (n=94 ROI, 4168 common genes, slide-z-score)  -  "
        "single model learned across both slides, no post-hoc sign alignment",
        fontsize=12, fontweight="bold", y=0.985)
    fig.savefig(HERE / "joint_cca_showcase.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[✓] joint_cca_showcase.png")


if __name__ == "__main__":
    main()
