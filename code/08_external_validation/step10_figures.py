#!/usr/bin/env python3
"""
Step 10 (figures) - Three-tier external validation. 600 DPI, caption-free.

Panels:
  (A) Tier 1 classification replication: AUC/AUPRC on GSE50760 with CIs.
  (B) Tier 2 survival validation: C-index per cohort with CIs; chance line.
  (C) Summary forest plot across all tiers (point + 95% CI).
  (D) Honest contrast: classification (strong) vs survival (weak) for the panel.

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

RES = Path("results/step10"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 14.5, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.labelweight": "bold",
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "legend.fontsize": 10.5, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42})


def main():
    r = json.load(open(RES/"external_validation_summary.json"))
    t1 = r["tier1_gse50760"]; g39 = r["tier2_gse39582"]; g17 = r["tier2_gse17536"]

    fig = plt.figure(figsize=(16,13))
    gs = gridspec.GridSpec(2,2,figure=fig,hspace=0.32,wspace=0.26)

    # (A) Tier 1 classification
    axA = fig.add_subplot(gs[0,0])
    metrics = ["AUC","AUPRC"]
    vals = [t1["auc"], t1["auprc"]]
    los = [t1["auc"]-t1["auc_ci"][0], t1["auprc"]-t1["auprc_ci"][0]]
    his = [t1["auc_ci"][1]-t1["auc"], t1["auprc_ci"][1]-t1["auprc"]]
    bars = axA.bar(metrics, vals, yerr=[los,his], capsize=7,
                   color=["#0072B2","#009E73"], edgecolor="black",
                   linewidth=1.1, error_kw={"elinewidth":1.3})
    for b,v in zip(bars,vals):
        axA.text(b.get_x()+b.get_width()/2, v-0.06, f"{v:.3f}",
                 ha="center", va="top", fontsize=12, fontweight="bold", color="white")
    axA.axhline(0.5, ls="--", color="red", lw=1.4, label="Chance")
    axA.set_ylabel("Score (95% CI)")
    axA.set_title(f"(A) Tier 1: classification replication\nGSE50760 (n={t1['n']}, RNA-seq)", loc="left")
    axA.set_ylim(0,1.05); axA.legend(frameon=True, loc="lower right"); axA.grid(axis="y",alpha=0.3)

    # (B) Tier 2 survival
    axB = fig.add_subplot(gs[0,1])
    cohorts = ["GSE39582","GSE17536"]
    cis = [g39["cindex"], g17["cindex"]]
    clos = [g39["cindex"]-g39["cindex_ci"][0], g17["cindex"]-g17["cindex_ci"][0]]
    chis = [g39["cindex_ci"][1]-g39["cindex"], g17["cindex_ci"][1]-g17["cindex"]]
    epvs = [g39["epv"], g17["epv"]]
    bars = axB.bar(cohorts, cis, yerr=[clos,chis], capsize=7,
                   color=["#E69F00","#CC79A7"], edgecolor="black",
                   linewidth=1.1, error_kw={"elinewidth":1.3})
    axB.axhline(0.5, ls="--", color="red", lw=1.6, label="Chance (C=0.5)")
    for b,v,ep in zip(bars,cis,epvs):
        axB.text(b.get_x()+b.get_width()/2, v+0.015, f"{v:.3f}\nEPV={ep:.1f}",
                 ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    axB.set_ylabel("Harrell's C-index (95% CI)")
    axB.set_title("(B) Tier 2: survival validation (signature transfer)", loc="left")
    axB.set_ylim(0.30,0.75); axB.legend(frameon=True); axB.grid(axis="y",alpha=0.3)

    # (C) forest plot all tiers
    axC = fig.add_subplot(gs[1,0])
    labels = ["GSE50760\n(AUC, classification)","GSE39582\n(C-index, survival)",
              "GSE17536\n(C-index, survival)"]
    points = [t1["auc"], g39["cindex"], g17["cindex"]]
    lo = [t1["auc_ci"][0], g39["cindex_ci"][0], g17["cindex_ci"][0]]
    hi = [t1["auc_ci"][1], g39["cindex_ci"][1], g17["cindex_ci"][1]]
    cols = ["#0072B2","#E69F00","#CC79A7"]
    ypos = np.arange(len(labels))[::-1]
    for y,p,l,h,c in zip(ypos,points,lo,hi,cols):
        axC.plot([l,h],[y,y],"-",color=c,lw=2.5)
        axC.plot(p,y,"o",color=c,markersize=12,markeredgecolor="black",markeredgewidth=0.7)
    axC.axvline(0.5, ls="--", color="red", lw=1.5, label="Chance / null")
    axC.set_yticks(ypos); axC.set_yticklabels(labels)
    axC.set_xlabel("Metric value (point + 95% CI)")
    axC.set_title("(C) External validation forest plot", loc="left")
    axC.set_xlim(0.35,1.05); axC.legend(frameon=True, loc="lower right"); axC.grid(axis="x",alpha=0.3)

    # (D) honest contrast
    axD = fig.add_subplot(gs[1,1])
    cats = ["Classification\n(tumour vs normal)","Survival\n(prognosis)"]
    # representative: tier-1 AUC vs best tier-2 c-index
    vals = [t1["auc"], max(g39["cindex"], g17["cindex"])]
    colors = ["#009E73","#D55E00"]
    bars = axD.bar(cats, vals, color=colors, edgecolor="black", linewidth=1.1)
    axD.axhline(0.5, ls="--", color="red", lw=1.5, label="Chance")
    for b,v in zip(bars,vals):
        axD.text(b.get_x()+b.get_width()/2, v+0.015, f"{v:.3f}",
                 ha="center", va="bottom", fontsize=12, fontweight="bold")
    axD.set_ylabel("External performance")
    axD.set_title("(D) Panel is a strong classifier, weak prognosticator", loc="left")
    axD.set_ylim(0,1.05); axD.legend(frameon=True); axD.grid(axis="y",alpha=0.3)

    fig.savefig(FIG/"Figure8_external_validation_tiers.png", dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG/"Figure8_external_validation_tiers.jpg", dpi=600,
                bbox_inches="tight", facecolor="white", pil_kwargs={"quality":95})
    print("Saved Figure8_external_validation_tiers.png / .jpg at 600 DPI")


if __name__ == "__main__":
    main()
