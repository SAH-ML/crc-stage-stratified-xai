# Step 3 — Preprocessing, Harmonisation & Quality Control

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

Step 3 takes the raw data downloaded in Steps 1–2 and produces clean, structured,
gene-harmonised matrices plus quality-control figures. **No outcome-dependent
transformation (batch correction, feature selection, supervised normalisation) is
performed here** — those are deferred to *within* cross-validation folds in Step 6
to prevent information leakage, exactly as specified in the Methods.

## Scripts (run in this order)

| Script | What it does |
|--------|--------------|
| `step3a_extract_microarray.py` | Parses GSE39582 and GSE17536 `family.soft.gz` files into gene-level expression matrices; maps Affymetrix GPL570 probes to gene symbols (maxMean collapse). Memory-efficient streaming parser. |
| `step3_preprocess.py` | Loads TCGA-COAD (expression + clinical + survival), prepares tumour/normal and AJCC-stage labels (sub-stages collapsed I–IV), QC-filters the training cohort (drops high-missingness/zero-variance genes), harmonises gene symbols across all four cohorts, and writes processed matrices. |
| `step3_qc_figures.py` | Generates **Figure 2** (cohort QC overview) at 600 DPI in PNG + JPG. |

## Inputs (from Steps 1–2, placed in `data/raw/`)

- `expression.tsv`, `clinical.tsv`, `survival.tsv` — TCGA-COAD (UCSC Xena)
- `GSE50760_expression_matrix.csv`, `GSE50760_metadata.csv` — Tier-1 RNA-seq
- `GSE39582_family.soft.gz`, `GSE39582_metadata.csv` — Tier-2 microarray
- `GSE17536_family.soft.gz`, `GSE17536_metadata.csv` — Tier-2 microarray

## Outputs (`data/processed/`)

- `TCGA_expression_qc.csv` — QC-filtered training expression (20,050 genes × 327 samples)
- `TCGA_labels.csv` — tumour/normal + AJCC stage per sample
- `TCGA_survival.csv` — overall-survival annotation
- `common_genes.csv` — 16,602 genes shared by all four cohorts
- `{COHORT}_expression_common.csv` — each cohort restricted to common genes

## Key results

| Cohort | Role | Samples | Genes (full) | Genes (common) |
|--------|------|---------|--------------|----------------|
| TCGA-COAD | Training | 327 (286 tumour / 41 normal) | 20,050 | 16,602 |
| GSE50760 | Tier 1 (RNA-seq) | 54 | 23,505 | 16,602 |
| GSE39582 | Tier 2 (microarray) | 585 | 22,880 | 16,602 |
| GSE17536 | Tier 2 (microarray) | 177 | 22,880 | 16,602 |

TCGA tumour stage distribution: **I = 45, II = 110, III = 80, IV = 39** (+12 unstaged).

**16,602 common genes** are retained across all cohorts — a large shared feature
space that confirms the decision to keep RNA-seq training separate from microarray
validation (rather than merging platforms, which would have collapsed the common
gene set dramatically).

## Figure 2 (key QC findings)

- **(C, E)** Tumour vs normal separate cleanly → strong learnable classification signal.
- **(D)** AJCC stages are intermixed on global PCA → stage differences are subtle and
  not captured by crude variance, **motivating the stage-stratified SHAP approach**.
- **(F)** Cohorts separate completely by platform on cross-cohort PCA → visual
  justification for not pooling platforms into training and for using the
  signature-transfer protocol in cross-platform validation.

## Reproducibility

- Random seed fixed at **42** throughout.
- Run from the project root: `python scripts/step3a_extract_microarray.py && python scripts/step3_preprocess.py && python scripts/step3_qc_figures.py`
- Python 3.10+; dependencies: `pandas numpy scikit-learn matplotlib seaborn`.
