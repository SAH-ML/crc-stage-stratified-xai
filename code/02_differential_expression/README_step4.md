# Step 4 — Differential Expression Analysis

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Two-implementation design (important)

Differential expression is implemented **twice**, with clearly defined roles:

| Script | Language | Role |
|--------|----------|------|
| `step4_DE_limma_voom.R` | R / Bioconductor | **CANONICAL / OFFICIAL.** Produces the DEG numbers reported in the paper, matching the Methods (limma-voom). |
| `step4_DE_python.py` | Python | **Development engine + reproducibility cross-check.** Faithful limma-style implementation (per-gene linear model + empirical-Bayes variance moderation) used to drive the pipeline end-to-end without R. |
| `step4_DE_concordance.py` | Python | **Trust anchor.** Compares the two implementations (set overlap, Jaccard, Spearman rank correlation) and writes `DE_concordance_report.csv` for the Supplementary Material. |

**Why two?** The development environment is Python-only, so the Python engine lets
the full pipeline run and the science be evaluated immediately. The canonical R
limma-voom is the method named in the manuscript; before submission it is run
(Phase 2), the two are shown to agree, and **all downstream steps are re-run on the
canonical limma-voom gene list** so every reported number traces to one source.

## Data-scale note (documented honestly)

UCSC Xena TCGA-COAD expression is `log2(norm_count + 1)`, *not* raw counts.
Classic voom expects raw counts; for an already-normalised log2 matrix the
recommended limma path is **limma-trend** (`eBayes(..., trend = TRUE)`), which the
R script uses. The Python engine mirrors this by modelling the log2 values
directly with empirical-Bayes moderation. The manuscript Methods will state this
exactly so the description matches the implementation.

## What it computes

- **Global** DE: tumour vs normal (all tumours pooled).
- **Per-stage** DE: Stage I/II/III/IV tumours vs normal (the candidate feature
  pools for stage-stratified SHAP).
- **Two FDR regimes** (multiple-testing safeguard from the Methods):
  per-contrast FDR and pooled FDR across all stage-by-gene tests.

## Results (Python development engine)

| Contrast | Significant DEGs | Up | Down |
|----------|------------------|----|----|
| Global tumour-vs-normal | 4,807 | 1,569 | 3,238 |
| Stage I vs normal | 4,996 | 1,511 | 3,485 |
| Stage II vs normal | 4,857 | 1,542 | 3,315 |
| Stage III vs normal | 4,661 | 1,608 | 3,053 |
| Stage IV vs normal | 4,971 | 1,809 | 3,162 |

DEGs shared across all four stages: **3,926**. Stage-unique DEGs peak at the
extremes — Stage I (292) and Stage IV (368) — an early, descriptive hint of the
stage-dependent signal the SHAP analysis will formally test.

Top global DEGs are well-established CRC genes (CDH3, KRT80, ETV4, ESM1 up;
OTOP2, BEST4, TMEFF2 down), confirming biological validity.

## Outputs

- `results/DE_python/DE_global_tumor_vs_normal.csv`
- `results/DE_python/DE_stage_{I,II,III,IV}_vs_normal.csv` (with per-contrast and pooled FDR)
- `results/DE_python/DE_summary.csv`
- `figures/Figure3_DE_overview.{png,jpg}` (600 DPI, caption-free)
- *(Phase 2)* `results/DE_limma_voom/*.csv` and `results/DE_concordance_report.csv`

## Reproducibility

- Random seed fixed at **42**.
- Python: `python scripts/step4_DE_python.py && python scripts/step4_DE_figures.py`
- R (Phase 2): run `scripts/step4_DE_limma_voom.R` in RStudio/Posit Cloud or
  Colab-R; Claude provides step-by-step guidance.
- Concordance (Phase 2): `python scripts/step4_DE_concordance.py`
