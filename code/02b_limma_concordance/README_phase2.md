# Phase 2 — Canonical R limma DE + Concordance with Python Engine

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

The paper's Methods specify differential expression (DE) by **limma**. This phase
runs the canonical limma analysis in R (Bioconductor) and verifies that the Python
development engine — used throughout the pipeline for tractability — agrees with
it. The paper reports the **canonical R limma numbers**; the Python engine is
retained as a validated reproducibility cross-check.

## How the R analysis was run

`step4_DE_limma_voom.R` (provided) was run on Posit Cloud (RStudio, R 4.6.0,
Bioconductor limma + edgeR). Inputs: `TCGA_expression_qc.csv` (20,530 genes ×
329 samples, log2 scale), `TCGA_labels.csv` (tumour/normal + AJCC stage).
Because UCSC-Xena TCGA expression is already log2-normalised (not raw counts),
limma's `lmFit` + `eBayes(trend=TRUE)` (limma-trend) path is used, which is the
recommended limma approach for log-expression matrices. Thresholds: FDR<0.05 and
|log2FC|>=1. Seed 42.

## Canonical R-limma result (reported in the paper)

Global tumour-vs-normal DEGs: **4,810** (1,567 up, 3,243 down).
Per stage: I=4,979; II=4,865; III=4,682; IV=4,971 (see `canonical_R_limma/DE_summary.csv`).

## Concordance: R limma vs Python engine

`step4_DE_concordance_final.py` quantifies agreement
(`DE_concordance_summary.csv`):

| Contrast | R DEGs | Python DEGs | Jaccard | log2FC Spearman | Direction agree |
|----------|--------|-------------|---------|-----------------|-----------------|
| Global | 4,810 | 4,807 | **0.995** | **1.000** | 100% |
| Stage I | 4,979 | 4,996 | 0.963 | 0.999 | 100% |
| Stage II | 4,865 | 4,857 | 0.989 | 1.000 | 100% |
| Stage III | 4,682 | 4,661 | 0.986 | 1.000 | 100% |
| Stage IV | 4,971 | 4,971 | **1.000** | 1.000 | 100% |

The two independent implementations agree **near-perfectly**: log2FC rank
correlation 1.000, direction agreement 100%, significant-DEG-set Jaccard
0.96–1.00. The few discordant genes lie at the FDR/log2FC significance boundary
and are immaterial.

## Panel-gene check (the critical one)

All 10 SHAP-panel genes are **significant in canonical R limma with matching
direction and near-identical fold-changes** (`DE_concordance_panel_genes.csv`):

| Gene | R log2FC | Py log2FC | R FDR | dir match |
|------|----------|-----------|-------|-----------|
| CDH3 | +6.548 | +6.544 | 8.4e-128 | ✓ |
| OTOP2 | -9.480 | -9.475 | 8.4e-128 | ✓ |
| KRT80 | +6.710 | +6.713 | 2.6e-109 | ✓ |
| ETV4 | +5.475 | +5.476 | 1.3e-101 | ✓ |
| ESM1 | +5.564 | +5.563 | 2.2e-99 | ✓ |
| DHRS7C | -2.757 | -2.756 | 1.2e-115 | ✓ |
| AADACL2 | -2.761 | -2.760 | 5.3e-114 | ✓ |
| LGI1 | -4.567 | -4.563 | 6.7e-90 | ✓ |
| TMEFF2 | -3.344 | -3.354 | 2.2e-82 | ✓ |
| KHDRBS2 | -2.474 | -2.472 | 1.2e-77 | ✓ |

**Implication:** the entire downstream pipeline (10-gene panel, model, SHAP,
stage-stratified findings, external validation) holds under canonical limma. No
reported result changes; the panel is not an artefact of the development engine.

## Outputs
- `canonical_R_limma/DE_*.csv` — the official limma DEG tables (paper source)
- `DE_concordance_summary.csv` — per-contrast agreement metrics
- `DE_concordance_panel_genes.csv` — panel-gene comparison
- `DE_concordance_verdict.json` — machine-readable verdict

## Reproducible code
```bash
# 1. (in R / Posit Cloud) produce canonical limma DEGs
Rscript step4_DE_limma_voom.R

# 2. (Python) concordance check vs the development engine
python step4_DE_concordance_final.py --rdir <R_outputs> --pydir <python_outputs>
```
Dependencies: R (limma, edgeR); Python (pandas, numpy, scipy). Deterministic; seed 42.
