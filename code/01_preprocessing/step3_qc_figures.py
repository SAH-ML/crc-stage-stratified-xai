#!/usr/bin/env python3
"""
Step 3 (figures) — Publication-quality QC visualisations at 600 DPI.

Generates:
  Figure 2. Cohort QC overview:
    (A) Sample composition per cohort (tumour/normal/stage)
    (B) Expression distribution boxplots (library-size / normalisation sanity)
    (C) PCA of the TCGA training cohort, coloured by tumour vs normal
    (D) PCA of the TCGA training cohort, coloured by AJCC stage
    (E) t-SNE of the TCGA training cohort, coloured by tumour vs normal
    (F) Cross-cohort PCA (common genes), coloured by cohort — illustrates the
        platform/batch structure that motivates keeping platforms separate.

All panels are rendered at 600 DPI in PNG and JPG with large, legible fonts,
clear axis labels, ticks, and legends, following the figure standards adopted
for the authors' prior work.

License: MIT   |   Random seed: 42
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

SEED = 42
np.random.seed(SEED)

PROC = Path("data/processed")
FIG = Path("figures")
FIG.mkdir(parents=True, exist_ok=True)

# ---- Global publication style (large, legible) ----
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "font.family": "DejaVu Sans",
    "font.size": 13,
    "axes.titlesize": 15,
    "axes.titleweight": "bold",
    "axes.labelsize": 13.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "legend.fontsize": 11,
    "legend.title_fontsize": 11.5,
    "axes.linewidth": 1.1,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Colour-blind-friendly palette
C_TUMOR = "#D55E00"
C_NORMAL = "#0072B2"
STAGE_COLORS = {"I": "#56B4E9", "II": "#009E73", "III": "#E69F00", "IV": "#CC79A7"}
COHORT_COLORS = {"TCGA-COAD": "#000000", "GSE50760": "#D55E00",
                 "GSE39582": "#009E73", "GSE17536": "#56B4E9"}


def load_all():
    expr = pd.read_csv(PROC / "TCGA_expression_qc.csv", index_col=0)
    labels = pd.read_csv(PROC / "TCGA_labels.csv", index_col=0)
    common = pd.read_csv(PROC / "common_genes.csv")["gene"].tolist()
    val = {}
    for name in ["GSE50760", "GSE39582", "GSE17536"]:
        val[name] = pd.read_csv(PROC / f"{name}_expression_common.csv", index_col=0)
    tcga_common = pd.read_csv(PROC / "TCGA_expression_common.csv", index_col=0)
    return expr, labels, common, val, tcga_common


def log_norm(df):
    """log2(x+1) for consistent scaling in plots (TCGA already log-ish; safe)."""
    return np.log2(df.clip(lower=0) + 1)


def main():
    expr, labels, common, val, tcga_common = load_all()

    # Align labels to expression columns
    labels = labels.loc[expr.columns]
    tumor_mask = labels["tumor_vs_normal"] == "tumor"

    # ---- Prepare matrices for dimensionality reduction (samples x genes) ----
    # Use top-2000 most variable genes for PCA/t-SNE stability & speed
    v = expr.var(axis=1).sort_values(ascending=False)
    top_genes = v.head(2000).index
    X = expr.loc[top_genes].T.values
    Xs = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=SEED)
    pcs = pca.fit_transform(Xs)
    pc_var = pca.explained_variance_ratio_ * 100

    tsne = TSNE(n_components=2, random_state=SEED, perplexity=30,
                init="pca", learning_rate="auto")
    tsne_xy = tsne.fit_transform(Xs)

    # ---- Cross-cohort PCA on common genes (illustrates platform structure) ----
    parts = [("TCGA-COAD", log_norm(tcga_common))]
    for name in ["GSE50760", "GSE39582", "GSE17536"]:
        parts.append((name, log_norm(val[name])))
    # build combined sample x gene matrix on common genes
    combined = []
    cohort_labels = []
    for name, mat in parts:
        m = mat.loc[common]               # genes x samples
        combined.append(m.T)              # samples x genes
        cohort_labels += [name] * m.shape[1]
    combined = pd.concat(combined, axis=0)
    cohort_labels = np.array(cohort_labels)
    Xc = StandardScaler().fit_transform(combined.values)
    pca_c = PCA(n_components=2, random_state=SEED)
    pcs_c = pca_c.fit_transform(Xc)
    pc_var_c = pca_c.explained_variance_ratio_ * 100

    # ======================= FIGURE =======================
    fig = plt.figure(figsize=(16, 19))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.34, wspace=0.24,
                           height_ratios=[1, 1, 1])

    # ---- (A) Sample composition ----
    axA = fig.add_subplot(gs[0, 0])
    comp = {
        "TCGA\nTumour": (tumor_mask).sum(),
        "TCGA\nNormal": (~tumor_mask).sum(),
        "GSE50760": val["GSE50760"].shape[1],
        "GSE39582": val["GSE39582"].shape[1],
        "GSE17536": val["GSE17536"].shape[1],
    }
    bars = axA.bar(list(comp.keys()), list(comp.values()),
                   color=[C_TUMOR, C_NORMAL, "#D55E00", "#009E73", "#56B4E9"],
                   edgecolor="black", linewidth=1.1)
    for b, v_ in zip(bars, comp.values()):
        axA.text(b.get_x() + b.get_width()/2, b.get_height() + 4, str(v_),
                 ha="center", va="bottom", fontsize=11.5, fontweight="bold")
    axA.set_ylabel("Number of samples")
    axA.set_title("(A) Sample composition across cohorts", loc="left")
    axA.set_ylim(0, max(comp.values()) * 1.15)
    axA.grid(axis="y", alpha=0.3)

    # ---- (B) Expression distribution boxplots (TCGA, subset of samples) ----
    axB = fig.add_subplot(gs[0, 1])
    sub = expr.iloc[:, :25]  # first 25 samples for legibility
    sub_log = log_norm(sub)
    bp = axB.boxplot([sub_log[c].values for c in sub_log.columns],
                     showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", linewidth=1.3))
    for patch, c in zip(bp["boxes"], sub_log.columns):
        is_t = labels.loc[c, "tumor_vs_normal"] == "tumor"
        patch.set_facecolor(C_TUMOR if is_t else C_NORMAL)
        patch.set_alpha(0.7)
    axB.set_ylabel("log\u2082(expression + 1)")
    axB.set_xlabel("TCGA samples (first 25 shown)")
    axB.set_title("(B) Expression distribution (normalisation check)", loc="left")
    axB.set_xticks([])
    axB.grid(axis="y", alpha=0.3)
    # legend
    from matplotlib.patches import Patch
    axB.legend(handles=[Patch(facecolor=C_TUMOR, edgecolor="black", label="Tumour"),
                        Patch(facecolor=C_NORMAL, edgecolor="black", label="Normal")],
               loc="upper right", frameon=True)

    # ---- (C) PCA coloured by tumour vs normal ----
    axC = fig.add_subplot(gs[1, 0])
    for lab, col in [("tumor", C_TUMOR), ("normal", C_NORMAL)]:
        m = labels["tumor_vs_normal"].values == lab
        axC.scatter(pcs[m, 0], pcs[m, 1], s=42, c=col, alpha=0.78,
                    edgecolors="black", linewidths=0.4,
                    label=lab.capitalize())
    axC.set_xlabel(f"PC1 ({pc_var[0]:.1f}% variance)")
    axC.set_ylabel(f"PC2 ({pc_var[1]:.1f}% variance)")
    axC.set_title("(C) PCA of TCGA cohort \u2014 tumour vs normal", loc="left")
    axC.legend(title="Tissue", loc="best", frameon=True)
    axC.grid(alpha=0.25)

    # ---- (D) PCA coloured by AJCC stage (tumours only) ----
    axD = fig.add_subplot(gs[1, 1])
    for stg in ["I", "II", "III", "IV"]:
        m = (labels["stage"].values == stg) & tumor_mask.values
        if m.sum() > 0:
            axD.scatter(pcs[m, 0], pcs[m, 1], s=46, c=STAGE_COLORS[stg],
                        alpha=0.82, edgecolors="black", linewidths=0.4,
                        label=f"Stage {stg} (n={m.sum()})")
    axD.set_xlabel(f"PC1 ({pc_var[0]:.1f}% variance)")
    axD.set_ylabel(f"PC2 ({pc_var[1]:.1f}% variance)")
    axD.set_title("(D) PCA of TCGA tumours \u2014 by AJCC stage", loc="left")
    axD.legend(title="AJCC stage", loc="best", frameon=True)
    axD.grid(alpha=0.25)

    # ---- (E) t-SNE coloured by tumour vs normal ----
    axE = fig.add_subplot(gs[2, 0])
    for lab, col in [("tumor", C_TUMOR), ("normal", C_NORMAL)]:
        m = labels["tumor_vs_normal"].values == lab
        axE.scatter(tsne_xy[m, 0], tsne_xy[m, 1], s=42, c=col, alpha=0.78,
                    edgecolors="black", linewidths=0.4, label=lab.capitalize())
    axE.set_xlabel("t-SNE dimension 1")
    axE.set_ylabel("t-SNE dimension 2")
    axE.set_title("(E) t-SNE of TCGA cohort \u2014 tumour vs normal", loc="left")
    axE.legend(title="Tissue", loc="best", frameon=True)
    axE.grid(alpha=0.25)

    # ---- (F) Cross-cohort PCA on common genes ----
    axF = fig.add_subplot(gs[2, 1])
    for name in ["TCGA-COAD", "GSE50760", "GSE39582", "GSE17536"]:
        m = cohort_labels == name
        axF.scatter(pcs_c[m, 0], pcs_c[m, 1], s=30,
                    c=COHORT_COLORS[name], alpha=0.6,
                    edgecolors="none", label=f"{name} (n={m.sum()})")
    axF.set_xlabel(f"PC1 ({pc_var_c[0]:.1f}% variance)")
    axF.set_ylabel(f"PC2 ({pc_var_c[1]:.1f}% variance)")
    axF.set_title("(F) Cross-cohort PCA (common genes) \u2014 platform structure",
                  loc="left")
    axF.legend(title="Cohort", loc="best", frameon=True, markerscale=1.5)
    axF.grid(alpha=0.25)

    fig.savefig(FIG / "Figure2_QC_overview.png", dpi=600, bbox_inches="tight",
                facecolor="white")
    fig.savefig(FIG / "Figure2_QC_overview.jpg", dpi=600, bbox_inches="tight",
                facecolor="white", pil_kwargs={"quality": 95})
    print("Saved Figure2_QC_overview.png / .jpg at 600 DPI")

    # report explained variance to results
    with open("results/step3_qc_stats.txt", "w") as f:
        f.write(f"TCGA PCA PC1 variance: {pc_var[0]:.2f}%\n")
        f.write(f"TCGA PCA PC2 variance: {pc_var[1]:.2f}%\n")
        f.write(f"Cross-cohort PCA PC1 variance: {pc_var_c[0]:.2f}%\n")
        f.write(f"Cross-cohort PCA PC2 variance: {pc_var_c[1]:.2f}%\n")
        f.write(f"Top variable genes used: {len(top_genes)}\n")
    print("Wrote results/step3_qc_stats.txt")


if __name__ == "__main__":
    main()
