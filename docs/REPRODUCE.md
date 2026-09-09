# Step-by-Step Reproduction Guide

This guide reproduces every result in the paper. Each step is independent and writes to
`results/`; the expected outputs are already committed, so you can verify your run.

**Prerequisites:** Python 3.10+ (`pip install -r requirements.txt`), and R 4.x with
Bioconductor (limma, edgeR) for the canonical DE step. Random seed 42 is fixed throughout.

---

## Step 0 — Get the data

Follow `docs/DATA_SOURCES.md`. Minimum to reproduce from scratch:
1. Download the three TCGA-COAD files from UCSC Xena (expression, clinical, survival).
2. Download GSE50760, GSE39582, GSE17536 series matrices from GEO.

A pre-processed `data/TCGA_labels.csv` is provided; the expression matrix is large and
is rebuilt by Step 1 (or download `TCGA_expression_qc.csv.gz` from the release assets).

## Step 1 — Preprocessing  (`code/01_preprocessing`)
```bash
python code/01_preprocessing/step3_preprocess.py
python code/01_preprocessing/step3a_extract_microarray.py   # for GEO validation cohorts
python code/01_preprocessing/step3_qc_figures.py            # Figure 2
```
Produces the QC'd expression matrix + labels. See its README for details.

## Step 2 — Differential expression  (`code/02_differential_expression`)
```bash
# Python engine (fast, for the pipeline):
python code/02_differential_expression/step4_DE_python.py
# Canonical R limma (the numbers reported in the paper):
Rscript code/02_differential_expression/step4_DE_limma_voom.R
python code/02_differential_expression/step4_DE_figures.py  # Figure 3
```
> If you are new to R, run the R script on Posit Cloud (https://posit.cloud) — upload
> the expression CSV + labels, open the .R file, click **Source**. ~10-min Bioconductor install.

## Step 2b — Concordance check  (`code/02b_limma_concordance`)
```bash
python code/02b_limma_concordance/step4_DE_concordance_final.py
```
Confirms Python vs R limma agree (global Jaccard 0.995, log2FC Spearman 1.000). Output in
`results/02b_limma_concordance/`.

## Step 3 — Network & enrichment  (`code/03_network_enrichment`)
The PPI network and enrichment are run on web tools (STRING, DAVID, KEGG, Reactome,
ClueGO — see DATA_SOURCES.md). Their exported outputs are in `results/03_network_enrichment/`.
Then:
```bash
python code/03_network_enrichment/step5_01_triangulated_enrichment.py
python code/03_network_enrichment/step5_02_two_level_analysis.py
python code/03_network_enrichment/step5_03_saudi_crossref.py
python code/03_network_enrichment/step5_figures.py          # Figure 9
```
Produces the 216 robust pathways, the two-level panel localisation, and the Saudi
cross-reference.

## Step 4 — Model benchmark  (`code/04_model_benchmark`)
```bash
python code/04_model_benchmark/step6_model_benchmark.py     # six algorithms
python code/04_model_benchmark/step6_bootstrap_ci.py        # 1000-iter CIs + imbalance analysis
python code/04_model_benchmark/step6_figures.py             # Figure 4
```
Selects the Elastic Net; produces Table 1 (class-imbalance evidence).

## Step 5 — Stage-stratified SHAP  (`code/05_stage_shap`)
```bash
python code/05_stage_shap/step7_stage_shap.py               # per-stage SHAP + concordance
python code/05_stage_shap/step7b_external_stage_validation.py  # GSE39582 replication
python code/05_stage_shap/step7c_transfer_diagnostic.py
python code/05_stage_shap/step7_figures.py                  # Figure 5
python code/05_stage_shap/step7b_figures.py                 # Figure 6
```
Produces the stage-importance gradient and its external replication (ESM1 late-stage).

## Step 6 — Biological corroboration  (`code/06_biological_corroboration`)
See its README — literature + GEPIA2 corroboration of the ten panel genes
(GEPIA2 run at http://gepia2.cancer-pku.cn/). Outputs in `results/06_biological_corroboration/`.

## Step 7 — Survival benchmark  (`code/07_survival`)
```bash
python code/07_survival/step9_survival_benchmark.py         # Cox / DeepSurv / RSF
python code/07_survival/step9_figures.py                    # Figure 7
```
Produces the honest null prognostic result (C-index ≈ 0.5).

## Step 8 — External validation  (`code/08_external_validation`)
```bash
python code/08_external_validation/step10_external_validation.py   # three-tier
python code/08_external_validation/step10_figures.py               # Figure 8
```

---

## Verifying your run

Each step prints a summary and writes CSV/JSON to `results/`. Compare against the
committed files — numbers should match to the precision reported in the paper. If a
step depends on a web-tool export (Steps 3, 6), the exports are already provided so the
downstream Python is fully reproducible.

## Figures

All 9 main + 3 supplementary figures (PNG + JPG, 600 dpi) are in `figures/`. The
figure-generation scripts above regenerate them from `results/`.
