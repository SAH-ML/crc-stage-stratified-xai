#!/usr/bin/env python3
"""
Step 9 (figures) - Survival benchmark visualisation. 600 DPI, caption-free.

Panels:
  (A) C-index per model with 95% bootstrap CI; chance line at 0.5.
  (B) Per-fold C-index distribution (strip + mean).
  (C) Kaplan-Meier by risk group (best model, median split) with log-rank p.
  (D) EPV diagnostic: events per stage vs the EPV>=10 threshold (honesty panel).

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
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

RES = Path("results/step9"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 14.5, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.labelweight": "bold",
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
    "legend.fontsize": 10.5, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42})

MC = {"CoxPH": "#0072B2", "RSF": "#009E73", "DeepSurv": "#D55E00"}


def main():
    res = pd.read_csv(RES/"survival_benchmark_overall.csv")
    km = pd.read_csv(RES/"risk_scores_overall.csv")
    epv = pd.read_csv(RES/"epv_table.csv")
    summary = json.load(open(RES/"survival_summary.json"))

    fig = plt.figure(figsize=(16,13))
    gs = gridspec.GridSpec(2,2,figure=fig,hspace=0.30,wspace=0.24)

    # (A) C-index with CI
    axA = fig.add_subplot(gs[0,0])
    order = res.sort_values("cindex", ascending=False)
    xs = np.arange(len(order))
    yerr = np.array([order["cindex"]-order["cindex_lo"],
                     order["cindex_hi"]-order["cindex"]])
    bars = axA.bar(xs, order["cindex"], yerr=yerr, capsize=6,
                   color=[MC[m] for m in order["model"]],
                   edgecolor="black", linewidth=1.1, error_kw={"elinewidth":1.2})
    axA.axhline(0.5, ls="--", color="red", lw=1.6, label="Chance (C=0.5)")
    for b,v in zip(bars, order["cindex"]):
        axA.text(b.get_x()+b.get_width()/2, v+0.015, f"{v:.3f}",
                 ha="center", va="bottom", fontsize=11, fontweight="bold")
    axA.set_xticks(xs); axA.set_xticklabels(order["model"])
    axA.set_ylabel("Harrell's C-index (95% CI)")
    axA.set_title("(A) Survival discrimination by model", loc="left")
    axA.set_ylim(0.30, 0.75); axA.legend(frameon=True); axA.grid(axis="y", alpha=0.3)

    # (B) per-fold strip
    axB = fig.add_subplot(gs[0,1])
    for i,m in enumerate(order["model"]):
        mean = order[order["model"]==m]["fold_cindex_mean"].values[0]
        sd = order[order["model"]==m]["fold_cindex_sd"].values[0]
        # reconstruct approx fold spread visually from mean+/-sd
        axB.errorbar(i, mean, yerr=sd, fmt="o", markersize=11,
                     color=MC[m], capsize=6, elinewidth=1.5,
                     markeredgecolor="black", markeredgewidth=0.6)
    axB.axhline(0.5, ls="--", color="red", lw=1.4)
    axB.set_xticks(range(len(order))); axB.set_xticklabels(order["model"])
    axB.set_ylabel("Per-fold C-index (mean \u00b1 SD)")
    axB.set_title("(B) Cross-validation fold consistency", loc="left")
    axB.set_ylim(0.30, 0.75); axB.grid(axis="y", alpha=0.3)

    # (C) Kaplan-Meier by risk group
    axC = fig.add_subplot(gs[1,0])
    kmf = KaplanMeierFitter()
    years = km["time"]/365.25
    for grp,col in [("low","#0072B2"),("high","#D55E00")]:
        mask = km["risk_group"]==grp
        kmf.fit(years[mask], km["event"][mask], label=f"{grp} risk (n={mask.sum()})")
        kmf.plot_survival_function(ax=axC, color=col, ci_show=True, linewidth=2.2)
    lr = logrank_test(years[km["risk_group"]=="low"], years[km["risk_group"]=="high"],
                      km["event"][km["risk_group"]=="low"], km["event"][km["risk_group"]=="high"])
    axC.set_xlabel("Time (years)"); axC.set_ylabel("Overall survival probability")
    axC.set_title(f"(C) KM by risk group, {summary['best_model_overall']} "
                  f"(log-rank p={lr.p_value:.3f})", loc="left")
    axC.legend(frameon=True, loc="lower left"); axC.grid(alpha=0.25)
    axC.set_ylim(0,1.02)

    # (D) EPV honesty panel
    axD = fig.add_subplot(gs[1,1])
    perstage = epv[epv["subgroup"].str.startswith("Stage")]
    bars = axD.bar(perstage["subgroup"], perstage["events"],
                   color="#999999", edgecolor="black", linewidth=1.1)
    axD.axhline(100, ls="--", color="red", lw=1.6,
                label="Events needed for EPV\u226510 (=100)")
    for b,ev in zip(bars, perstage["events"]):
        axD.text(b.get_x()+b.get_width()/2, ev+1.5, str(int(ev)),
                 ha="center", va="bottom", fontsize=11, fontweight="bold")
    axD.set_ylabel("Number of death events")
    axD.set_title("(D) EPV diagnostic: all stages below threshold \u2192 exploratory",
                  loc="left")
    axD.legend(frameon=True); axD.grid(axis="y", alpha=0.3)
    axD.set_ylim(0, 110)

    fig.savefig(FIG/"Figure7_survival_benchmark.png", dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG/"Figure7_survival_benchmark.jpg", dpi=600,
                bbox_inches="tight", facecolor="white", pil_kwargs={"quality":95})
    print("Saved Figure7_survival_benchmark.png / .jpg at 600 DPI")


if __name__ == "__main__":
    main()
