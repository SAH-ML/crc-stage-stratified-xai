#!/usr/bin/env python3
"""
Step 9b - Stage-wise bootstrap confidence intervals for the SHAP concordance
statistics (revision analysis -> Supplementary Table S8).

Addresses peer-review comments (Reviewer 1 Comment 1; Reviewer 2 Comment 2) asking
for confidence intervals / sensitivity analysis showing the stage-specific SHAP
rankings are robust and not driven by a few influential samples.

Method: the per-gene stage importance vectors produced by Step 7
(results/05_stage_shap/shap_importance_by_stage.csv) are resampled with replacement
over the ten panel genes (1,000 iterations, seed 42); for each resample the Spearman
rank correlation between a stage pair (or stage vs global) is recomputed, yielding a
95% bootstrap interval around each published concordance point estimate. Point
estimates and permutation p-values are those already reported by Step 7
(results/05_stage_shap/concordance.json), so this analysis is fully consistent with
the main-text values.

Inputs : results/05_stage_shap/shap_importance_by_stage.csv,
         results/05_stage_shap/concordance.json
Outputs: results/09_revision_analyses/tableS8_stage_bootstrap_ci.csv

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 42
rng = np.random.RandomState(SEED)
N_BOOT = 1000

RS7 = Path("results/05_stage_shap")
RES = Path("results/09_revision_analyses")
RES.mkdir(parents=True, exist_ok=True)

COL = {"global": "global", "I": "stage_I", "II": "stage_II",
       "III": "stage_III", "IV": "stage_IV"}
ORDER = [("global", "I"), ("global", "II"), ("global", "III"), ("global", "IV"),
         ("I", "II"), ("I", "III"), ("I", "IV"),
         ("II", "III"), ("II", "IV"), ("III", "IV")]


def boot_ci(a_vec, b_vec):
    n = len(a_vec)
    rhos = []
    for _ in range(N_BOOT):
        idx = rng.randint(0, n, n)
        r, _ = spearmanr(a_vec[idx], b_vec[idx])
        rhos.append(r)
    rhos = np.array(rhos)
    return np.nanpercentile(rhos, 2.5), np.nanpercentile(rhos, 97.5)


def main():
    imp = pd.read_csv(RS7 / "shap_importance_by_stage.csv")
    conc = json.load(open(RS7 / "concordance.json"))
    rows = []
    for a, b in ORDER:
        key = f"global_vs_{b}" if a == "global" else f"{a}_vs_{b}"
        rho = conc[key]["rho"]
        p = conc[key]["p"]
        lo, hi = boot_ci(imp[COL[a]].values, imp[COL[b]].values)
        la = "Global" if a == "global" else f"Stage {a}"
        rows.append({"Comparison": f"{la} vs Stage {b}",
                     "Spearman rho": f"{rho:.3f}",
                     "95% CI": f"[{lo:.3f}, {hi:.3f}]",
                     "p-value": f"{p:.3f}"})
        print(f"{la} vs Stage {b}: rho={rho:.3f} [{lo:.3f}, {hi:.3f}] p={p:.3f}")
    pd.DataFrame(rows).to_csv(RES / "tableS8_stage_bootstrap_ci.csv", index=False)
    print(f"Saved {RES/'tableS8_stage_bootstrap_ci.csv'}")


if __name__ == "__main__":
    main()
