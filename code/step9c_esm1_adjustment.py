#!/usr/bin/env python3
"""
Step 9c - Robustness of the ESM1 stage association to adjustment (revision analysis
-> Supplementary Table S11).

Addresses peer-review comments (Reviewer 1 Comment 7; Reviewer 2 Comment 5) asking
whether ESM1 retains its stage-associated signal after adjustment for other
covariates, and to distinguish predictive attribution from biological causality.

Method (tractable with bulk public data):
  (1) Unadjusted Spearman correlation of ESM1 expression with ordinal AJCC stage.
  (2) Partial correlation of ESM1 vs stage, controlling for the other nine panel
      genes (residualise both ESM1 and stage on the other genes by OLS, then
      correlate residuals). Tests for an independent stage-associated component.
  (3) Mann-Whitney comparison of ESM1 expression in late (III+IV) vs early (I+II).

Full adjustment for tumour-microenvironment composition (stromal/immune fractions)
would require cell-type deconvolution data not available for this bulk cohort and
is identified in the manuscript as a next step.

Inputs : data/TCGA_expression_qc.csv.gz, data/TCGA_labels.csv,
         results/04_model_benchmark/final_feature_panel.csv
Outputs: results/09_revision_analyses/esm1_adjustment.json

License: MIT   |   Random seed: 42
"""
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import spearmanr, mannwhitneyu

warnings.filterwarnings("ignore")
SEED = 42
np.random.seed(SEED)

DATA = Path("data")
RESB = Path("results/04_model_benchmark")
RES = Path("results/09_revision_analyses")
RES.mkdir(parents=True, exist_ok=True)

STAGE_MAP = {"I": 1, "II": 2, "III": 3, "IV": 4}


def residualise(y, Z):
    b, _, _, _ = lstsq(Z, y, rcond=None)
    return y - Z @ b


def main():
    expr = pd.read_csv(DATA / "TCGA_expression_qc.csv.gz", index_col=0)
    labels = pd.read_csv(DATA / "TCGA_labels.csv", index_col=0)
    panel = pd.read_csv(RESB / "final_feature_panel.csv")["gene"].tolist()

    X = expr.loc[panel].T
    meta = labels.loc[X.index]
    mask = (meta["tumor_vs_normal"] == "tumor") & (meta["stage"].isin(STAGE_MAP))
    Xt = X[mask]
    st = meta.loc[Xt.index, "stage"].map(STAGE_MAP).values.astype(float)

    # (1) unadjusted
    rho_u, p_u = spearmanr(Xt["ESM1"].values, st)

    # (2) partial, controlling the other nine panel genes
    others = [g for g in panel if g != "ESM1"]
    Z = np.column_stack([np.ones(len(Xt)), Xt[others].values])
    e_res = residualise(Xt["ESM1"].values.astype(float), Z)
    s_res = residualise(st, Z)
    rho_a, p_a = spearmanr(e_res, s_res)

    # (3) late vs early
    early = Xt.loc[meta.loc[Xt.index, "stage"].isin(["I", "II"]), "ESM1"]
    late = Xt.loc[meta.loc[Xt.index, "stage"].isin(["III", "IV"]), "ESM1"]
    _, p_el = mannwhitneyu(late, early, alternative="greater")

    out = {
        "esm1_vs_stage_unadjusted": {"rho": round(float(rho_u), 3), "p": float(p_u)},
        "esm1_vs_stage_adjusted_partial": {
            "rho": round(float(rho_a), 3), "p": float(p_a),
            "note": "partial correlation controlling for the other nine panel genes"},
        "esm1_late_vs_early": {
            "median_early": round(float(early.median()), 3),
            "median_late": round(float(late.median()), 3),
            "mwu_p": float(p_el)},
    }
    json.dump(out, open(RES / "esm1_adjustment.json", "w"), indent=2)
    print(json.dumps(out, indent=2))
    print(f"Saved {RES/'esm1_adjustment.json'}")


if __name__ == "__main__":
    main()
