#!/usr/bin/env python3
"""
Step 6 (figures) - Multi-metric model benchmark, imbalance diagnosis, DeLong
selection, and feature stability. 600 DPI, caption-free.

Panels:
  (A) Multi-metric comparison (AUC, AUPRC, MCC, balanced acc) per model.
  (B) Per-fold AUC distribution (strip + mean).
  (C) ROC curves (pooled out-of-fold predictions).
  (D) Feature-selection frequency across folds (stability).

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
from sklearn.metrics import roc_curve, auc

RES = Path("results/step6"); FIG = Path("figures"); FIG.mkdir(exist_ok=True)
plt.rcParams.update({
    "savefig.dpi": 600, "font.family": "DejaVu Sans", "font.size": 13,
    "axes.titlesize": 15, "axes.titleweight": "bold",
    "axes.labelsize": 13.5, "axes.labelweight": "bold",
    "xtick.labelsize": 11, "ytick.labelsize": 11.5,
    "legend.fontsize": 10.5, "axes.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42})

MC = {"ElasticNet": "#0072B2", "RandomForest": "#009E73", "SVM": "#E69F00",
      "XGBoost": "#D55E00", "LightGBM": "#CC79A7", "MLP": "#56B4E9"}


def main():
    summary = pd.read_csv(RES / "model_benchmark_summary.csv")
    perfold = pd.read_csv(RES / "per_fold_metrics.csv")
    roc_data = json.load(open(RES / "roc_data.json"))
    freq = pd.read_csv(RES / "feature_selection_frequency.csv", index_col=0)
    unw = summary[summary["variant"] == "unweighted"].set_index("model")
    models = list(MC.keys())

    # bootstrap CIs (for error bars) if available
    ci_path = RES / "bootstrap_ci.csv"
    ci = pd.read_csv(ci_path) if ci_path.exists() else None

    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.22)

    # (A) grouped multi-metric bars with 95% bootstrap CI error bars
    axA = fig.add_subplot(gs[0, 0])
    metrics = ["auc_mean", "auprc_mean", "mcc_mean", "balanced_acc_mean"]
    ci_metric = {"auc_mean": "auc", "auprc_mean": "auprc",
                 "mcc_mean": "mcc", "balanced_acc_mean": "balanced_acc"}
    mlabels = ["AUC", "AUPRC", "MCC", "Bal.Acc"]
    x = np.arange(len(models)); w = 0.2
    mcolors = ["#0072B2", "#009E73", "#D55E00", "#E69F00"]
    for i, (met, ml, col) in enumerate(zip(metrics, mlabels, mcolors)):
        vals = [unw.loc[m, met] for m in models]
        errs = None
        if ci is not None:
            lo, hi = [], []
            for m in models:
                row = ci[(ci["model"] == m) & (ci["metric"] == ci_metric[met])]
                if len(row):
                    pt = row["point"].values[0]
                    lo.append(pt - row["ci_low"].values[0])
                    hi.append(row["ci_high"].values[0] - pt)
                else:
                    lo.append(0); hi.append(0)
            errs = np.array([lo, hi])
        axA.bar(x + (i - 1.5) * w, vals, w, label=ml, color=col,
                edgecolor="black", linewidth=0.7,
                yerr=errs, capsize=2, error_kw={"elinewidth": 0.8})
    axA.set_xticks(x); axA.set_xticklabels(models, rotation=20)
    axA.set_ylabel("Score (95% bootstrap CI)")
    axA.set_title("(A) Multi-metric comparison (imbalance-robust, 95% CI)", loc="left")
    axA.set_ylim(0.55, 1.05); axA.legend(ncol=4, loc="lower center", frameon=True)
    axA.grid(axis="y", alpha=0.3)

    # (B) per-fold AUC strip (unweighted)
    axB = fig.add_subplot(gs[0, 1])
    pf = perfold[perfold["variant"] == "unweighted"]
    for i, m in enumerate(models):
        vals = pf.loc[pf["model"] == m, "auc"].values
        jit = np.random.RandomState(42).normal(0, 0.04, len(vals))
        axB.scatter(np.full(len(vals), i) + jit, vals, s=55, color=MC[m],
                    edgecolors="black", linewidths=0.5, alpha=0.85, zorder=3)
        axB.hlines(vals.mean(), i - 0.25, i + 0.25, color="black", lw=2, zorder=4)
    axB.set_xticks(range(len(models))); axB.set_xticklabels(models, rotation=20)
    axB.set_ylabel("Per-fold AUC")
    axB.set_title("(B) Per-fold AUC distribution (bar = mean)", loc="left")
    axB.set_ylim(0.85, 1.02); axB.grid(axis="y", alpha=0.3)

    # (C) ROC
    axC = fig.add_subplot(gs[1, 0])
    for m in models:
        ys, ss = roc_data[m]["y"], roc_data[m]["p"]
        fpr, tpr, _ = roc_curve(ys, ss)
        axC.plot(fpr, tpr, lw=2, color=MC[m], label=f"{m} (AUC={auc(fpr,tpr):.3f})")
    axC.plot([0, 1], [0, 1], "--", color="grey", lw=1)
    axC.set_xlabel("False positive rate"); axC.set_ylabel("True positive rate")
    axC.set_title("(C) ROC curves (pooled out-of-fold)", loc="left")
    axC.legend(loc="lower right", frameon=True); axC.grid(alpha=0.25)

    # (D) feature stability
    axD = fig.add_subplot(gs[1, 1])
    fs = freq.sort_values("folds_selected", ascending=True)
    cols = ["#009E73" if v == fs["folds_selected"].max() else
            ("#56B4E9" if v >= 3 else "#BBBBBB") for v in fs["folds_selected"]]
    axD.barh(fs.index, fs["folds_selected"], color=cols,
             edgecolor="black", linewidth=0.9)
    axD.set_xlabel("Folds selected in (out of 5)")
    axD.set_title("(D) Feature-selection stability across folds", loc="left")
    axD.set_xlim(0, 5.5); axD.axvline(3, ls="--", color="black", lw=1, alpha=0.6)
    axD.grid(axis="x", alpha=0.3)

    fig.savefig(FIG / "Figure4_model_benchmark.png", dpi=600,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure4_model_benchmark.jpg", dpi=600,
                bbox_inches="tight", facecolor="white", pil_kwargs={"quality": 95})
    print("Saved Figure4_model_benchmark.png / .jpg at 600 DPI")


if __name__ == "__main__":
    main()
