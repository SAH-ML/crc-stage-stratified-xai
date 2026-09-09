#!/usr/bin/env python3
"""
Step 3 — Preprocessing, label preparation, gene harmonisation, and QC.

This script:
  1. Loads the TCGA-COAD training cohort (expression + clinical + survival).
  2. Prepares tumour-vs-normal and AJCC-stage labels (sub-stages collapsed).
  3. Loads the three validation cohorts (gene-level matrices from Step 3a /
     Step 2) and their metadata.
  4. Performs gene-symbol harmonisation across all cohorts (common-gene index).
  5. Applies basic QC filtering (near-zero-variance and high-missingness genes)
     to the TRAINING cohort only.
  6. Saves processed matrices and label tables for downstream steps.

IMPORTANT (leakage prevention): no batch correction, feature selection, or
normalisation that depends on the outcome is performed here. ComBat and feature
selection are deferred to *within* cross-validation folds in Step 6, exactly as
specified in the Methods. This script performs only outcome-independent cleaning
and structuring.

Part of the reproducible pipeline for:
  "Stage-Stratified Explainable AI Reveals Shifting Transcriptomic Drivers
   of Colorectal Cancer Progression"

Author: [Author list]
License: MIT
Random seed: 42 (fixed throughout the project)
"""

from pathlib import Path
import pandas as pd
import numpy as np

SEED = 42
np.random.seed(SEED)

RAW = Path("data/raw")
PROC = Path("data/processed")
RESULTS = Path("results")
PROC.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. TCGA-COAD training cohort
# ---------------------------------------------------------------------------

def collapse_stage(raw):
    """Collapse AJCC sub-stages to parent stage (I, II, III, IV)."""
    if pd.isna(raw):
        return None
    s = str(raw).upper().replace("STAGE", "").strip()
    if "DISCREPANCY" in s:
        return None
    for roman in ["IV", "III", "II", "I"]:
        if s.startswith(roman):
            return roman
    return None


def load_tcga():
    print("=== TCGA-COAD (training cohort) ===")
    expr = pd.read_csv(RAW / "expression.tsv", sep="\t", index_col=0)
    print(f"  Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

    clin = pd.read_csv(RAW / "clinical.tsv", sep="\t", low_memory=False)
    clin = clin.set_index("sampleID")

    # Sample type: tumour vs normal (TCGA barcodes: -01 tumour, -11 normal)
    clin["stage_collapsed"] = clin["pathologic_stage"].apply(collapse_stage)

    # Align clinical to expression samples
    common = [s for s in expr.columns if s in clin.index]
    expr = expr[common]
    clin = clin.loc[common]

    # Label tumour vs normal from sample_type
    is_tumor = clin["sample_type"] == "Primary Tumor"
    is_normal = clin["sample_type"] == "Solid Tissue Normal"

    labels = pd.DataFrame(index=common)
    labels["sample_type"] = clin["sample_type"]
    labels["tumor_vs_normal"] = np.where(is_tumor, "tumor",
                                  np.where(is_normal, "normal", "other"))
    labels["stage"] = clin["stage_collapsed"]

    # Restrict to tumour + normal (drop metastatic/recurrent for the main design)
    keep = labels["tumor_vs_normal"].isin(["tumor", "normal"])
    expr = expr.loc[:, keep]
    labels = labels.loc[keep]

    print(f"  After tumour/normal restriction: {expr.shape[1]} samples")
    print(f"    Tumour: {(labels['tumor_vs_normal']=='tumor').sum()}, "
          f"Normal: {(labels['tumor_vs_normal']=='normal').sum()}")
    print(f"  Stage distribution (tumours): ")
    print(labels.loc[labels['tumor_vs_normal']=='tumor','stage']
          .value_counts(dropna=False).to_string())

    # Survival
    surv = pd.read_csv(RAW / "survival.tsv", sep="\t", index_col="sample")
    return expr, labels, surv


# ---------------------------------------------------------------------------
# 2. QC filtering (training cohort only; outcome-independent)
# ---------------------------------------------------------------------------

def qc_filter(expr, missing_frac=0.10, var_quantile=0.0):
    """
    Remove genes with >missing_frac missing values and zero/near-zero variance.
    Outcome-independent, so safe to apply outside CV folds.
    """
    n0 = expr.shape[0]
    # drop high-missingness genes
    miss = expr.isna().mean(axis=1)
    expr = expr.loc[miss <= missing_frac]
    # fill any residual missing with row (gene) median
    expr = expr.apply(lambda r: r.fillna(r.median()), axis=1)
    # drop zero-variance genes
    variances = expr.var(axis=1)
    expr = expr.loc[variances > 0]
    print(f"  QC: {n0} -> {expr.shape[0]} genes "
          f"(removed {n0 - expr.shape[0]} high-missing/zero-variance)")
    return expr


# ---------------------------------------------------------------------------
# 3. Validation cohorts
# ---------------------------------------------------------------------------

def load_validation():
    cohorts = {}

    # GSE50760 (RNA-seq, Tier 1) — gene symbols already
    g = pd.read_csv(RAW / "GSE50760_expression_matrix.csv", index_col=0)
    m = pd.read_csv(RAW / "GSE50760_metadata.csv")
    cohorts["GSE50760"] = {"expr": g, "meta": m}
    print(f"\n=== GSE50760 (Tier 1, RNA-seq) ===")
    print(f"  Expression: {g.shape[0]} genes x {g.shape[1]} samples")

    # GSE39582 (microarray, Tier 2) — from Step 3a
    g = pd.read_csv(PROC / "GSE39582_expression_gene.csv", index_col=0)
    m = pd.read_csv(RAW / "GSE39582_metadata.csv")
    cohorts["GSE39582"] = {"expr": g, "meta": m}
    print(f"=== GSE39582 (Tier 2, microarray) ===")
    print(f"  Expression: {g.shape[0]} genes x {g.shape[1]} samples")

    # GSE17536 (microarray, Tier 2)
    g = pd.read_csv(PROC / "GSE17536_expression_gene.csv", index_col=0)
    m = pd.read_csv(RAW / "GSE17536_metadata.csv")
    cohorts["GSE17536"] = {"expr": g, "meta": m}
    print(f"=== GSE17536 (Tier 2, microarray) ===")
    print(f"  Expression: {g.shape[0]} genes x {g.shape[1]} samples")

    return cohorts


# ---------------------------------------------------------------------------
# 4. Gene harmonisation
# ---------------------------------------------------------------------------

def harmonise_genes(tcga_expr, cohorts):
    """Find the common gene-symbol set across all cohorts."""
    gene_sets = [set(tcga_expr.index)]
    for name, c in cohorts.items():
        gene_sets.append(set(c["expr"].index))
    common = set.intersection(*gene_sets)
    print(f"\n=== Gene harmonisation ===")
    print(f"  TCGA genes: {len(gene_sets[0])}")
    for name, c in cohorts.items():
        print(f"  {name} genes: {c['expr'].shape[0]}")
    print(f"  Common genes across ALL cohorts: {len(common)}")
    return sorted(common)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    tcga_expr, tcga_labels, tcga_surv = load_tcga()
    print("\n=== QC filtering (TCGA training cohort) ===")
    tcga_expr = qc_filter(tcga_expr)

    cohorts = load_validation()
    common_genes = harmonise_genes(tcga_expr, cohorts)

    # Save processed training matrix (full gene set; common set saved separately)
    tcga_expr.to_csv(PROC / "TCGA_expression_qc.csv")
    tcga_labels.to_csv(PROC / "TCGA_labels.csv")
    tcga_surv.to_csv(PROC / "TCGA_survival.csv")

    # Save common-gene restricted matrices for cross-cohort work
    pd.Series(common_genes, name="gene").to_csv(
        PROC / "common_genes.csv", index=False)
    tcga_expr.loc[common_genes].to_csv(PROC / "TCGA_expression_common.csv")
    for name, c in cohorts.items():
        c["expr"].loc[common_genes].to_csv(
            PROC / f"{name}_expression_common.csv")

    # Summary table
    summary = pd.DataFrame({
        "cohort": ["TCGA-COAD", "GSE50760", "GSE39582", "GSE17536"],
        "role": ["Training", "Tier1 (RNA-seq)", "Tier2 (microarray)",
                 "Tier2 (microarray)"],
        "n_samples": [tcga_expr.shape[1], cohorts["GSE50760"]["expr"].shape[1],
                      cohorts["GSE39582"]["expr"].shape[1],
                      cohorts["GSE17536"]["expr"].shape[1]],
        "n_genes_full": [tcga_expr.shape[0],
                         cohorts["GSE50760"]["expr"].shape[0],
                         cohorts["GSE39582"]["expr"].shape[0],
                         cohorts["GSE17536"]["expr"].shape[0]],
        "n_genes_common": [len(common_genes)] * 4,
    })
    summary.to_csv(RESULTS / "step3_cohort_summary.csv", index=False)
    print("\n=== Step 3 summary ===")
    print(summary.to_string(index=False))
    print(f"\nProcessed files written to {PROC}/")


if __name__ == "__main__":
    main()
