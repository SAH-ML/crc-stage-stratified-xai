#!/usr/bin/env python3
"""
Step 7b (figures) - External validation visualisation (GSE39582). 600 DPI, caption-free.

Panels:
  (A) Cross-cohort per-stage concordance (TCGA vs GSE39582), with significance.
  (B) ESM1 importance-rank trajectory in both cohorts (the replicated signal).
  (C) Side-by-side stage concordance matrices (TCGA | GSE39582).
  (D) Per-stage importance scatter, TCGA vs GSE39582 (Stage III, best-replicating).

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
from scipy.stats import spearmanr

RES = Path("results/step7"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 14.5, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.labelweight": "bold",
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "legend.fontsize": 10.5, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42})

STAGES = ["I","II","III","IV"]
SC = {"I":"#56B4E9","II":"#009E73","III":"#E69F00","IV":"#CC79A7"}
PANEL = ['CDH3','OTOP2','DHRS7C','AADACL2','KRT80','ETV4','ESM1','TMEFF2','LGI1','KHDRBS2']


def main():
    ext = json.load(open(RES/"external_validation_gse39582.json"))
    tcga = pd.read_csv(RES/"shap_importance_by_stage.csv").set_index("gene").loc[PANEL]
    gse = pd.read_csv(RES/"gse39582_shap_importance_by_stage.csv").set_index("gene").loc[PANEL]
    cross = ext["cross_cohort_perstage_concordance"]

    fig = plt.figure(figsize=(16,13))
    gs = gridspec.GridSpec(2,2,figure=fig,hspace=0.32,wspace=0.26)

    # (A) cross-cohort concordance bars
    axA = fig.add_subplot(gs[0,0])
    rhos=[cross[s]["rho"] for s in STAGES]; ps=[cross[s]["p"] for s in STAGES]
    bars=axA.bar([f"Stage {s}" for s in STAGES], rhos,
                 color=[SC[s] for s in STAGES], edgecolor="black", linewidth=1.1)
    for b,r,p in zip(bars,rhos,ps):
        star = "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.05 else "n.s."))
        axA.text(b.get_x()+b.get_width()/2, r+0.02, f"{r:.2f}\n{star}",
                 ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    axA.axhline(0.7, ls="--", color="red", lw=1.4, label="Strong agreement (\u03c1=0.7)")
    axA.set_ylabel("Cross-cohort Spearman \u03c1\n(TCGA vs GSE39582)")
    axA.set_title("(A) Per-stage importance replication", loc="left")
    axA.set_ylim(0,1.08); axA.legend(frameon=True, loc="lower right"); axA.grid(axis="y",alpha=0.3)

    # (B) ESM1 trajectory both cohorts
    axB = fig.add_subplot(gs[0,1])
    esm1 = ext["key_gene_trajectories"]["ESM1"]
    x=range(len(STAGES))
    axB.plot(x,[esm1["tcga"][s] for s in STAGES],"-o",lw=2.6,markersize=9,
             color="#0072B2",label="TCGA (discovery)")
    axB.plot(x,[esm1["gse39582"][s] for s in STAGES],"-s",lw=2.6,markersize=9,
             color="#D55E00",label="GSE39582 (validation)")
    axB.set_xticks(list(x)); axB.set_xticklabels([f"Stage {s}" for s in STAGES])
    axB.set_ylabel("ESM1 importance rank\n(1 = most important)")
    axB.invert_yaxis()
    axB.set_title("(B) ESM1 rises to late-stage prominence in BOTH cohorts", loc="left")
    axB.legend(frameon=True); axB.grid(alpha=0.25)

    # (C) concordance matrices side by side (use insets)
    axC = fig.add_subplot(gs[1,0])
    # build TCGA and GSE matrices
    def cmat(impdf, prefix):
        M=np.eye(4)
        for i,a in enumerate(STAGES):
            for j,b in enumerate(STAGES):
                if i<j:
                    r=spearmanr(impdf[f"{prefix}{a}"], impdf[f"{prefix}{b}"])[0]
                    M[i,j]=r; M[j,i]=r
        return M
    Mt=cmat(tcga,"stage_"); Mg=cmat(gse,"stage_")
    # plot TCGA in left half, GSE in right via two imshows is messy; show difference
    Mdiff = Mg - Mt
    im=axC.imshow(Mdiff, cmap="RdBu_r", vmin=-0.4, vmax=0.4)
    axC.set_xticks(range(4)); axC.set_xticklabels(STAGES)
    axC.set_yticks(range(4)); axC.set_yticklabels(STAGES)
    for i in range(4):
        for j in range(4):
            axC.text(j,i,f"{Mdiff[i,j]:+.2f}",ha="center",va="center",
                     fontsize=10,fontweight="bold")
    axC.set_title("(C) Concordance difference (GSE39582 \u2212 TCGA)", loc="left")
    axC.set_xlabel("AJCC stage"); axC.set_ylabel("AJCC stage")
    cb=fig.colorbar(im,ax=axC,fraction=0.046,pad=0.04); cb.set_label("\u0394\u03c1",fontsize=10)

    # (D) Stage III scatter (best replicating)
    axD = fig.add_subplot(gs[1,1])
    s="III"
    axD.scatter(tcga[f"stage_{s}"], gse[f"stage_{s}"], s=90,
                color="#E69F00", edgecolors="black", linewidths=0.6, zorder=3)
    for g in PANEL:
        axD.annotate(g,(tcga.loc[g,f"stage_{s}"], gse.loc[g,f"stage_{s}"]),
                     fontsize=8.5, xytext=(4,3), textcoords="offset points")
    axD.set_xlabel("TCGA mean|SHAP| (Stage III)")
    axD.set_ylabel("GSE39582 mean|contribution| (Stage III)")
    axD.set_title(f"(D) Stage III importance agreement (\u03c1={cross['III']['rho']:.2f})", loc="left")
    axD.grid(alpha=0.25)

    fig.savefig(FIG/"Figure6_external_validation.png",dpi=600,
                bbox_inches="tight",facecolor="white")
    fig.savefig(FIG/"Figure6_external_validation.jpg",dpi=600,
                bbox_inches="tight",facecolor="white",pil_kwargs={"quality":95})
    print("Saved Figure6_external_validation.png / .jpg at 600 DPI")


if __name__ == "__main__":
    main()
