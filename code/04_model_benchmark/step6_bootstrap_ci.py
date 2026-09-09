#!/usr/bin/env python3
"""
Step 6 (bootstrap CIs) - 95% confidence intervals for all benchmark metrics.

Reads the pooled out-of-fold predictions saved by step6_model_benchmark.py and
computes bias-corrected percentile bootstrap 95% confidence intervals for AUC,
AUPRC, MCC, balanced accuracy, sensitivity, and specificity, for every model.

This implements the uncertainty quantification committed in Methods 2.7.2
(1,000-iteration bootstrap, percentile 95% CI) and establishes the bootstrap-CI
machinery reused in later steps (SHAP stability, survival, external validation).

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             matthews_corrcoef, balanced_accuracy_score,
                             recall_score, confusion_matrix)

SEED = 42
rng = np.random.RandomState(SEED)
RES = Path("results/step6")
N_BOOT = 1000


def metric_set(y, p, thr=0.5):
    pred = (np.asarray(p) >= thr).astype(int)
    y = np.asarray(y)
    out = {}
    # AUC / AUPRC require both classes present in the resample
    if len(np.unique(y)) < 2:
        return None
    out["auc"] = roc_auc_score(y, p)
    out["auprc"] = average_precision_score(y, p)
    out["mcc"] = matthews_corrcoef(y, pred)
    out["balanced_acc"] = balanced_accuracy_score(y, pred)
    out["sensitivity"] = recall_score(y, pred, pos_label=1, zero_division=0)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    out["specificity"] = tn / (tn + fp) if (tn + fp) else np.nan
    return out


def bootstrap_ci(y, p, n_boot=N_BOOT):
    y = np.asarray(y); p = np.asarray(p)
    n = len(y)
    point = metric_set(y, p)
    boot = {k: [] for k in point}
    tries = 0
    collected = 0
    while collected < n_boot and tries < n_boot * 5:
        tries += 1
        idx = rng.randint(0, n, n)             # resample with replacement
        ms = metric_set(y[idx], p[idx])
        if ms is None:                         # skip resamples missing a class
            continue
        for k, v in ms.items():
            boot[k].append(v)
        collected += 1
    rows = {}
    for k in point:
        arr = np.array(boot[k])
        lo, hi = np.percentile(arr, [2.5, 97.5])
        rows[k] = (point[k], lo, hi)
    return rows, collected


def main():
    roc_data = json.load(open(RES / "roc_data.json"))
    models = list(roc_data.keys())

    records = []
    for m in models:
        y, p = roc_data[m]["y"], roc_data[m]["p"]
        cis, ncollected = bootstrap_ci(y, p)
        for metric, (pt, lo, hi) in cis.items():
            records.append({
                "model": m, "metric": metric,
                "point": round(pt, 4),
                "ci_low": round(lo, 4), "ci_high": round(hi, 4),
                "ci_string": f"{pt:.3f} [{lo:.3f}-{hi:.3f}]",
                "n_bootstrap": ncollected,
            })
        print(f"{m}: AUC {cis['auc'][0]:.3f} "
              f"[{cis['auc'][1]:.3f}-{cis['auc'][2]:.3f}] "
              f"(n_boot={ncollected})", flush=True)

    df = pd.DataFrame(records)
    df.to_csv(RES / "bootstrap_ci.csv", index=False)

    # wide table: one row per model, AUC/AUPRC/MCC CI strings
    wide = df.pivot(index="model", columns="metric", values="ci_string")
    wide.to_csv(RES / "bootstrap_ci_wide.csv")
    print("\n=== 95% bootstrap CIs (point [low-high]) ===")
    print(wide[["auc", "auprc", "mcc", "balanced_acc"]].to_string())
    print(f"\nSaved bootstrap_ci.csv and bootstrap_ci_wide.csv to {RES}/")


if __name__ == "__main__":
    main()
