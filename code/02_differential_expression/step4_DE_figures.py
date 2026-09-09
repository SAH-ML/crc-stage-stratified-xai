#!/usr/bin/env python3
"""
Step 4 (figures) — Differential expression visualisations at 600 DPI.

Generates (caption-free; captions belong in the manuscript text):
  Figure 3. Differential expression overview:
    (A) Volcano plot, global tumour vs normal.
    (B) Number of significant DEGs per AJCC stage (up/down).
    (C) Upset-style overlap of DEG sets across the four stages.
    (D) Volcano plot for Stage IV vs normal (representative per-stage panel).

All panels rendered at 600 DPI in PNG + JPG, large legible fonts.

License: MIT   |   Random seed: 42
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from adjustText import adjust_text

RES = Path("results/DE_python")
FIG = Path("figures")
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 15, "axes.titleweight": "bold",
    "axes.labelsize": 13.5, "axes.labelweight": "bold",
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "legend.fontsize": 11, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

C_UP = "#D55E00"
C_DOWN = "#0072B2"
C_NS = "#BBBBBB"
STAGE_COLORS = {"I": "#56B4E9", "II": "#009E73", "III": "#E69F00", "IV": "#CC79A7"}
LOGFC, FDR = 1.0, 0.05


def volcano(ax, df, title, label_top=8):
    fdr = df["FDR_per_contrast"].fillna(1.0).values
    lfc = df["log2FC"].values
    neglog = -np.log10(np.clip(fdr, 1e-300, 1))
    sig_up = (fdr < FDR) & (lfc >= LOGFC)
    sig_dn = (fdr < FDR) & (lfc <= -LOGFC)
    ns = ~(sig_up | sig_dn)
    ax.scatter(lfc[ns], neglog[ns], s=8, c=C_NS, alpha=0.4, edgecolors="none")
    ax.scatter(lfc[sig_up], neglog[sig_up], s=14, c=C_UP, alpha=0.7,
               edgecolors="none", label=f"Up ({sig_up.sum()})")
    ax.scatter(lfc[sig_dn], neglog[sig_dn], s=14, c=C_DOWN, alpha=0.7,
               edgecolors="none", label=f"Down ({sig_dn.sum()})")
    ax.axvline(LOGFC, ls="--", c="grey", lw=1)
    ax.axvline(-LOGFC, ls="--", c="grey", lw=1)
    ax.axhline(-np.log10(FDR), ls="--", c="grey", lw=1)
    # label a few top genes, with automatic de-overlap
    top = df.assign(neglog=neglog).nlargest(label_top, "neglog")
    texts = []
    for _, row in top.iterrows():
        x = row["log2FC"]
        y = -np.log10(max(row["FDR_per_contrast"], 1e-300))
        texts.append(ax.text(x, y, row["gene"], fontsize=9, fontweight="bold"))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
                expand_points=(1.4, 1.6), expand_text=(1.3, 1.5),
                force_text=(0.4, 0.6), only_move={"points": "y", "text": "xy"})
    ax.set_xlabel("log\u2082 fold-change (tumour vs normal)")
    ax.set_ylabel("\u2212log\u2081\u2080(FDR)")
    ax.set_title(title, loc="left")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(alpha=0.2)


def main():
    g = pd.read_csv(RES / "DE_global_tumor_vs_normal.csv")
    stage_dfs = {s: pd.read_csv(RES / f"DE_stage_{s}_vs_normal.csv")
                 for s in ["I", "II", "III", "IV"]}

    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.22)

    # (A) global volcano
    axA = fig.add_subplot(gs[0, 0])
    volcano(axA, g, "(A) Global: tumour vs normal")

    # (B) DEG counts per stage
    axB = fig.add_subplot(gs[0, 1])
    stages = ["I", "II", "III", "IV"]
    ups, downs = [], []
    for s in stages:
        d = stage_dfs[s]
        sig = d["significant_per_contrast"]
        ups.append(int((sig & (d["log2FC"] > 0)).sum()))
        downs.append(int((sig & (d["log2FC"] < 0)).sum()))
    x = np.arange(len(stages))
    w = 0.4
    b1 = axB.bar(x - w/2, ups, w, color=C_UP, edgecolor="black",
                 linewidth=1, label="Up-regulated")
    b2 = axB.bar(x + w/2, downs, w, color=C_DOWN, edgecolor="black",
                 linewidth=1, label="Down-regulated")
    for bars in (b1, b2):
        for b in bars:
            axB.text(b.get_x()+b.get_width()/2, b.get_height()+30,
                     str(int(b.get_height())), ha="center", va="bottom",
                     fontsize=10, fontweight="bold")
    axB.set_xticks(x)
    axB.set_xticklabels([f"Stage {s}" for s in stages])
    axB.set_ylabel("Number of significant DEGs")
    axB.set_title("(B) Significant DEGs per AJCC stage", loc="left")
    axB.legend(frameon=True)
    axB.grid(axis="y", alpha=0.3)
    axB.set_ylim(0, max(max(ups), max(downs))*1.18)

    # (C) DEG set overlap across stages (bar of shared/unique)
    axC = fig.add_subplot(gs[1, 0])
    sig_sets = {}
    for s in stages:
        d = stage_dfs[s]
        sig_sets[s] = set(d.loc[d["significant_per_contrast"], "gene"])
    # shared across all four
    shared_all = set.intersection(*sig_sets.values())
    # unique to each stage
    unique_counts = {}
    for s in stages:
        others = set.union(*[sig_sets[o] for o in stages if o != s])
        unique_counts[s] = len(sig_sets[s] - others)
    cats = ["Shared\n(all 4 stages)"] + [f"Unique to\nStage {s}" for s in stages]
    vals = [len(shared_all)] + [unique_counts[s] for s in stages]
    cols = ["#555555"] + [STAGE_COLORS[s] for s in stages]
    bars = axC.bar(cats, vals, color=cols, edgecolor="black", linewidth=1)
    for b in bars:
        axC.text(b.get_x()+b.get_width()/2, b.get_height()+10,
                 str(int(b.get_height())), ha="center", va="bottom",
                 fontsize=10.5, fontweight="bold")
    axC.set_ylabel("Number of DEGs")
    axC.set_title("(C) Stage-shared vs stage-unique DEGs", loc="left")
    axC.grid(axis="y", alpha=0.3)
    axC.set_ylim(0, max(vals)*1.18)

    # (D) representative per-stage volcano (Stage IV)
    axD = fig.add_subplot(gs[1, 1])
    volcano(axD, stage_dfs["IV"], "(D) Stage IV vs normal")

    fig.savefig(FIG / "Figure3_DE_overview.png", dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure3_DE_overview.jpg", dpi=600,
                bbox_inches="tight", facecolor="white",
                pil_kwargs={"quality": 95})
    print("Saved Figure3_DE_overview.png / .jpg at 600 DPI")

    # write a small stats file
    with open("results/step4_de_stats.txt", "w") as f:
        f.write(f"Global DEGs: {int(g['significant'].sum())}\n")
        for s in stages:
            d = stage_dfs[s]
            f.write(f"Stage {s} DEGs (per-contrast): "
                    f"{int(d['significant_per_contrast'].sum())}\n")
        f.write(f"DEGs shared across all 4 stages: {len(shared_all)}\n")
        for s in stages:
            f.write(f"DEGs unique to Stage {s}: {unique_counts[s]}\n")
    print("Wrote results/step4_de_stats.txt")


if __name__ == "__main__":
    main()
