#!/usr/bin/env python3
"""
Step 7 (figures) - Stage-stratified SHAP visualisation. 600 DPI, caption-free.

Panels:
  (A) Heatmap: mean|SHAP| per gene (rows) x {global, Stage I-IV} (cols).
  (B) Rank-shift (bump) plot: gene importance rank across stages.
  (C) Stage-pair concordance matrix (Spearman rho).
  (D) Per-stage SHAP ranking stability (bootstrap mean rho with threshold line).

License: MIT | Random seed: 42
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

RES = Path("results/step7"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 15, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.labelweight": "bold",
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "legend.fontsize": 10.5, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42})

STAGES = ["I", "II", "III", "IV"]
SC = {"I": "#56B4E9", "II": "#009E73", "III": "#E69F00", "IV": "#CC79A7"}


def main():
    imp = pd.read_csv(RES / "shap_importance_by_stage.csv")
    genes = imp["gene"].tolist()
    concord = json.load(open(RES / "concordance.json"))
    stab = json.load(open(RES / "stage_stability.json"))

    # order genes by global importance
    imp = imp.sort_values("global", ascending=False).reset_index(drop=True)
    genes = imp["gene"].tolist()

    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.28)

    # (A) heatmap of mean|SHAP|
    axA = fig.add_subplot(gs[0, 0])
    cols = ["global"] + [f"stage_{s}" for s in STAGES]
    M = imp[cols].values
    # normalise each column to [0,1] for visual comparability of *ranking*
    Mn = M / M.max(axis=0, keepdims=True)
    im = axA.imshow(Mn, aspect="auto", cmap="viridis")
    axA.set_xticks(range(len(cols)))
    axA.set_xticklabels(["Global", "I", "II", "III", "IV"])
    axA.set_yticks(range(len(genes))); axA.set_yticklabels(genes)
    axA.set_xlabel("Stage"); axA.set_title("(A) Mean |SHAP| by stage (col-normalised)", loc="left")
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04)
    cb.set_label("Relative importance", fontsize=10)

    # (B) rank-shift bump plot
    axB = fig.add_subplot(gs[0, 1])
    ranks = {}
    for s in STAGES:
        col = imp[f"stage_{s}"].values
        order = np.argsort(col)[::-1]
        r = np.empty(len(genes), dtype=int)
        for rank, gi in enumerate(order):
            r[gi] = rank + 1
        ranks[s] = r
    xpos = range(len(STAGES))
    # highlight genes whose rank changes the most
    rank_matrix = np.array([ranks[s] for s in STAGES]).T  # genes x stages
    rank_range = rank_matrix.max(axis=1) - rank_matrix.min(axis=1)
    highlight = set(np.argsort(rank_range)[::-1][:4])
    for gi, g in enumerate(genes):
        y = [ranks[s][gi] for s in STAGES]
        if gi in highlight:
            axB.plot(xpos, y, "-o", lw=2.4, markersize=7, label=g, zorder=3)
        else:
            axB.plot(xpos, y, "-", lw=1, color="#CCCCCC", zorder=1)
    axB.set_xticks(list(xpos)); axB.set_xticklabels([f"Stage {s}" for s in STAGES])
    axB.set_ylabel("Importance rank (1 = most important)")
    axB.invert_yaxis()
    axB.set_title("(B) Importance rank shift across stages", loc="left")
    axB.legend(title="Largest shifts", loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=True)
    axB.grid(alpha=0.25)

    # (C) concordance matrix
    axC = fig.add_subplot(gs[1, 0])
    allstages = STAGES
    Cm = np.eye(len(allstages))
    for i, a in enumerate(allstages):
        for j, b in enumerate(allstages):
            if i == j:
                continue
            key = f"{a}_vs_{b}" if f"{a}_vs_{b}" in concord else f"{b}_vs_{a}"
            if key in concord:
                Cm[i, j] = concord[key]["rho"]
    im2 = axC.imshow(Cm, cmap="RdYlGn", vmin=0.5, vmax=1.0)
    axC.set_xticks(range(len(allstages))); axC.set_xticklabels(allstages)
    axC.set_yticks(range(len(allstages))); axC.set_yticklabels(allstages)
    for i in range(len(allstages)):
        for j in range(len(allstages)):
            axC.text(j, i, f"{Cm[i,j]:.2f}", ha="center", va="center",
                     fontsize=11, fontweight="bold",
                     color="black")
    axC.set_title("(C) Stage-pair SHAP concordance (Spearman \u03c1)", loc="left")
    axC.set_xlabel("AJCC stage"); axC.set_ylabel("AJCC stage")
    fig.colorbar(im2, ax=axC, fraction=0.046, pad=0.04)

    # (D) stability bars
    axD = fig.add_subplot(gs[1, 1])
    rhos = [stab[s]["mean_rho"] for s in STAGES]
    sds = [stab[s]["sd_rho"] for s in STAGES]
    bars = axD.bar([f"Stage {s}" for s in STAGES], rhos, yerr=sds, capsize=5,
                   color=[SC[s] for s in STAGES], edgecolor="black", linewidth=1.1)
    axD.axhline(0.7, ls="--", color="red", lw=1.5, label="Stability threshold (\u03c1=0.7)")
    for b, r in zip(bars, rhos):
        axD.text(b.get_x()+b.get_width()/2, r+0.02, f"{r:.2f}",
                 ha="center", va="bottom", fontsize=11, fontweight="bold")
    axD.set_ylabel("Bootstrap mean Spearman \u03c1")
    axD.set_title("(D) Per-stage SHAP ranking stability", loc="left")
    axD.set_ylim(0, 1.05); axD.legend(frameon=True); axD.grid(axis="y", alpha=0.3)

    fig.savefig(FIG / "Figure5_stage_shap.png", dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure5_stage_shap.jpg", dpi=600,
                bbox_inches="tight", facecolor="white", pil_kwargs={"quality": 95})
    print("Saved Figure5_stage_shap.png / .jpg at 600 DPI")


if __name__ == "__main__":
    main()
