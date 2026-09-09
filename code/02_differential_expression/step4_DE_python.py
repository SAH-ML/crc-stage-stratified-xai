#!/usr/bin/env python3
"""
Step 4 (Python) — Differential Expression Analysis: development engine.

This is the PYTHON limma-equivalent used to drive pipeline development. It
implements the core statistical logic of limma — per-gene linear modelling with
empirical-Bayes moderated variance — adapted for the already-log2-normalised
TCGA-COAD expression matrix from UCSC Xena.

IMPORTANT — relationship to the canonical method:
  The OFFICIAL differentially expressed gene (DEG) lists reported in the paper
  are produced by the canonical R limma-voom script (step4_DE_limma_voom.R).
  This Python implementation is the development engine and the reproducibility
  cross-check; concordance between the two is quantified by
  step4_DE_concordance.py and reported in the Supplementary Material. Downstream
  steps in the FINAL paper are re-run on the canonical limma-voom gene list.

Statistical approach (limma-style, for log2 data):
  * Fit, per gene, an ordinary-least-squares linear model of expression on the
    contrast of interest.
  * Compute the ordinary t-statistic and residual variance per gene.
  * Apply empirical-Bayes variance moderation (Smyth 2004): shrink each gene's
    residual variance toward a pooled prior estimated across all genes, giving a
    moderated t-statistic with augmented degrees of freedom. This is the same
    shrinkage principle limma uses, reproduced here in Python.
  * Convert moderated t to p-values; apply Benjamini-Hochberg FDR.

Two FDR regimes are reported, per the Methods (multiple-testing safeguard):
  (i) per-contrast FDR (within each comparison), and
  (ii) pooled FDR across all stage-stratified gene-by-contrast tests.

Analyses performed:
  1. Global: tumour vs normal (all tumours pooled).
  2. Per-stage: Stage_k tumours vs normal, for k in {I, II, III, IV}.

License: MIT   |   Random seed: 42
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy import special

SEED = 42
np.random.seed(SEED)

PROC = Path("data/processed")
RES = Path("results/DE_python")
RES.mkdir(parents=True, exist_ok=True)

LOGFC_THRESHOLD = 1.0     # |log2 fold-change| >= 1
FDR_THRESHOLD = 0.05      # FDR < 0.05


# ---------------------------------------------------------------------------
# Empirical-Bayes moderation (limma / Smyth 2004 core)
# ---------------------------------------------------------------------------

def squeeze_var(s2, df):
    """
    Empirical-Bayes variance moderation (limma's squeezeVar).
    Estimates a prior (d0, s0^2) from the distribution of per-gene residual
    variances via the method of moments on log-variances, then returns the
    posterior (moderated) variances and the prior degrees of freedom.

    Parameters
    ----------
    s2 : per-gene residual variances (array)
    df : residual degrees of freedom (scalar, equal across genes here)

    Returns
    -------
    s2_post : moderated variances
    df_prior : prior degrees of freedom (d0)
    df_total : df + df_prior (moderated t denominator df)
    """
    s2 = np.asarray(s2, dtype=float)
    s2 = np.where(s2 <= 0, np.nan, s2)
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)  # adjust to log-chisq mean
    emean = np.nanmean(e)
    evar = np.nanvar(e, ddof=1)

    # variance of log residual variance under chisq model = trigamma(df/2)
    evar_pred = special.polygamma(1, df / 2)
    if evar <= evar_pred:
        # near-infinite prior df -> heavy shrinkage to common variance
        df_prior = np.inf
        s2_prior = np.exp(emean + special.digamma(df / 2) - np.log(df / 2))
        s2_post = np.full_like(s2, s2_prior)
        return s2_post, df_prior, np.inf

    # solve trigamma(df_prior/2) = evar - trigamma(df/2)
    target = evar - evar_pred
    # invert trigamma numerically
    def trigamma_inverse(x):
        # Newton iteration (limma's trigammaInverse)
        y = 0.5 + 1.0 / x
        for _ in range(50):
            tg = special.polygamma(1, y)
            d = special.polygamma(2, y)
            delta = tg * (1 - tg / x) / d
            y += delta
            if np.all(np.abs(delta) < 1e-8):
                break
        return y
    df_prior = 2 * trigamma_inverse(target)
    s2_prior = np.exp(emean + special.digamma(df_prior / 2) - np.log(df_prior / 2))

    s2_post = (df_prior * s2_prior + df * s2) / (df_prior + df)
    df_total = df + df_prior
    return s2_post, df_prior, df_total


def moderated_ttest(expr, group_mask):
    """
    Two-group moderated t-test (limma-style) on log2 expression.

    expr : genes x samples DataFrame (log2 scale)
    group_mask : boolean array over samples; True = group1 (e.g. tumour),
                 False = group2 (e.g. normal)

    Returns a DataFrame with log2FC, moderated t, p-value.
    """
    X = expr.values
    g1 = X[:, group_mask]
    g2 = X[:, ~group_mask]
    n1, n2 = g1.shape[1], g2.shape[1]

    mean1 = g1.mean(axis=1)
    mean2 = g2.mean(axis=1)
    logfc = mean1 - mean2  # already log2 scale -> difference is log2FC

    # pooled residual variance
    var1 = g1.var(axis=1, ddof=1)
    var2 = g2.var(axis=1, ddof=1)
    df = n1 + n2 - 2
    sp2 = ((n1 - 1) * var1 + (n2 - 1) * var2) / df
    se_factor = np.sqrt(1.0 / n1 + 1.0 / n2)

    # empirical-Bayes moderation of variance
    sp2_mod, df_prior, df_total = squeeze_var(sp2, df)

    se_mod = np.sqrt(sp2_mod) * se_factor
    with np.errstate(divide="ignore", invalid="ignore"):
        t_mod = logfc / se_mod
    # moderated p-value uses augmented df
    if np.isinf(df_total):
        p = 2 * stats.norm.sf(np.abs(t_mod))
    else:
        p = 2 * stats.t.sf(np.abs(t_mod), df_total)

    out = pd.DataFrame({
        "gene": expr.index,
        "log2FC": logfc,
        "moderated_t": t_mod,
        "p_value": p,
    }).set_index("gene")
    return out


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR, robust to NaN p-values (genes with undefined
    test statistics, e.g. zero variance in a subgroup, are assigned FDR=NaN and
    excluded from the ranking denominator)."""
    p = np.asarray(pvals, dtype=float)
    fdr = np.full(len(p), np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return fdr
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr_valid = np.empty(n)
    fdr_valid[order] = np.clip(ranked, 0, 1)
    fdr[valid] = fdr_valid
    return fdr


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    expr = pd.read_csv(PROC / "TCGA_expression_qc.csv", index_col=0)
    labels = pd.read_csv(PROC / "TCGA_labels.csv", index_col=0).loc[expr.columns]

    tumor = labels["tumor_vs_normal"].values == "tumor"
    normal = labels["tumor_vs_normal"].values == "normal"
    stages = labels["stage"].values

    all_results = {}

    # ---- 1. Global: tumour vs normal ----
    print("=== Global DE: tumour vs normal ===")
    res = moderated_ttest(expr, tumor)
    res["FDR_per_contrast"] = bh_fdr(res["p_value"].values)
    res["significant"] = ((res["FDR_per_contrast"] < FDR_THRESHOLD) &
                          (res["log2FC"].abs() >= LOGFC_THRESHOLD))
    res = res.sort_values("p_value")
    res.to_csv(RES / "DE_global_tumor_vs_normal.csv")
    n_sig = res["significant"].sum()
    print(f"  Significant DEGs (FDR<{FDR_THRESHOLD}, |log2FC|>={LOGFC_THRESHOLD}): {n_sig}")
    print(f"    Up: {(res['significant'] & (res['log2FC']>0)).sum()}, "
          f"Down: {(res['significant'] & (res['log2FC']<0)).sum()}")
    all_results["global"] = res

    # ---- 2. Per-stage: Stage_k vs normal ----
    pooled_p = []          # collect all p-values for pooled FDR
    pooled_index = []
    for stg in ["I", "II", "III", "IV"]:
        stage_mask = (stages == stg) & tumor
        n_stage = stage_mask.sum()
        print(f"\n=== Per-stage DE: Stage {stg} (n={n_stage}) vs normal ===")
        # build subset: this stage's tumours + all normals
        keep = stage_mask | normal
        sub_expr = expr.loc[:, keep]
        sub_group = stage_mask[keep]   # True = this-stage tumour
        res_s = moderated_ttest(sub_expr, sub_group)
        res_s["FDR_per_contrast"] = bh_fdr(res_s["p_value"].values)
        res_s["significant_per_contrast"] = (
            (res_s["FDR_per_contrast"] < FDR_THRESHOLD) &
            (res_s["log2FC"].abs() >= LOGFC_THRESHOLD))
        res_s = res_s.sort_values("p_value")
        all_results[f"stage_{stg}"] = res_s
        # collect for pooled FDR
        pooled_p.append(res_s["p_value"].values)
        pooled_index += [(stg, g) for g in res_s.index]
        n_sig_s = res_s["significant_per_contrast"].sum()
        print(f"  Significant (per-contrast FDR): {n_sig_s}")

    # ---- pooled FDR across ALL stage-by-gene tests ----
    print("\n=== Pooled multiple-testing correction across all stage tests ===")
    pooled_p_all = np.concatenate(pooled_p)
    pooled_fdr_all = bh_fdr(pooled_p_all)
    # assign pooled FDR back to each stage result
    offset = 0
    for stg in ["I", "II", "III", "IV"]:
        res_s = all_results[f"stage_{stg}"]
        n = len(res_s)
        # pooled_index order matches concatenation order; map by position
        res_s = res_s.copy()
        # recover this stage's slice (same order as appended)
        # NOTE: res_s was sorted; re-align by reindexing on original order
        # Simplest: recompute pooled FDR per stage using the global pooled threshold
        offset += n
    # Re-do cleanly: build one long frame then split
    long_rows = []
    for stg in ["I", "II", "III", "IV"]:
        r = all_results[f"stage_{stg}"][["log2FC", "p_value"]].copy()
        r["stage"] = stg
        r["gene"] = r.index
        long_rows.append(r)
    long_df = pd.concat(long_rows, ignore_index=True)
    long_df["FDR_pooled"] = bh_fdr(long_df["p_value"].values)
    long_df["significant_pooled"] = (
        (long_df["FDR_pooled"] < FDR_THRESHOLD) &
        (long_df["log2FC"].abs() >= LOGFC_THRESHOLD))
    # write per-stage files WITH pooled FDR merged in
    for stg in ["I", "II", "III", "IV"]:
        sub = long_df[long_df["stage"] == stg].set_index("gene")
        res_s = all_results[f"stage_{stg}"].join(sub[["FDR_pooled", "significant_pooled"]])
        res_s.to_csv(RES / f"DE_stage_{stg}_vs_normal.csv")
        print(f"  Stage {stg}: per-contrast sig={res_s['significant_per_contrast'].sum()}, "
              f"pooled sig={res_s['significant_pooled'].sum()}")

    # ---- summary table ----
    summary = []
    g = all_results["global"]
    summary.append({"contrast": "Global tumour-vs-normal",
                    "n_significant": int(g["significant"].sum()),
                    "n_up": int((g["significant"] & (g["log2FC"]>0)).sum()),
                    "n_down": int((g["significant"] & (g["log2FC"]<0)).sum())})
    for stg in ["I", "II", "III", "IV"]:
        sub = long_df[long_df["stage"] == stg]
        per = all_results[f"stage_{stg}"]["significant_per_contrast"].sum()
        summary.append({"contrast": f"Stage {stg} vs normal",
                        "n_significant": int(per),
                        "n_up": int((per>0) and ((all_results[f'stage_{stg}']['significant_per_contrast']) & (all_results[f'stage_{stg}']['log2FC']>0)).sum()),
                        "n_down": int(((all_results[f'stage_{stg}']['significant_per_contrast']) & (all_results[f'stage_{stg}']['log2FC']<0)).sum())})
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RES / "DE_summary.csv", index=False)
    print("\n=== Step 4 (Python) summary ===")
    print(summary_df.to_string(index=False))
    print(f"\nResults written to {RES}/")


if __name__ == "__main__":
    main()
