# Stage-Stratified Explainable AI for Colorectal Cancer Transcriptomics

Reproducible code and data for the manuscript:

> **Stage-Stratified Explainable Machine Learning Identifies a Reproducible
> Transcriptomic Classifier with Stage-Dependent Feature Importance in Colorectal
> Cancer: A Multi-Cohort Study with Saudi Population Contextualisation**
>
> Asif Hassan Syed¹, Sultan Alhayyani² — King Abdulaziz University
> ¹ Dept. of Computer Science, FCIT Rabigh; ² Dept. of Chemistry, FSA Rabigh
> Correspondence: shassan1@kau.edu.sa · ORCID 0000-0002-7288-3098

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 1. What this repository contains

A complete, end-to-end reproduction of every analysis in the paper: from raw public
transcriptomic data, through differential expression, network/enrichment analysis, a
six-algorithm classification benchmark, stage-stratified SHAP explainability, survival
benchmarking, and three-tier external validation, to the final figures and tables.

```
.
├── code/        # all analysis scripts, organised by pipeline step
├── data/        # processed inputs + pointers to raw public sources
├── results/     # exact numerical outputs (CSV/JSON) for every step
├── figures/     # all manuscript figures (PNG + JPG, 600 dpi)
├── docs/        # data-source details, environment, step-by-step guide
├── requirements.txt
└── LICENSE      # MIT
```

## 2. Key result (one-line summary)

A ten-gene panel (CDH3, OTOP2, DHRS7C, AADACL2, KRT80, ETV4, ESM1, TMEFF2, LGI1,
KHDRBS2) classifies tumour vs normal with external AUC ≈ 0.97; its SHAP feature-
importance structure shifts measurably and **reproducibly** across AJCC stages
(ESM1 rising to late-stage prominence in two cohorts); the panel is a strong
classifier but **not** a prognostic model (survival C-index ≈ 0.5).

## 3. Data sources (all public — see docs/DATA_SOURCES.md for full detail)

| Cohort | Role | Platform | Access |
|--------|------|----------|--------|
| **TCGA-COAD** | Training / discovery | RNA-seq | UCSC Xena — https://xenabrowser.net/ (hub: https://tcga.xenahubs.net/) |
| **GSE50760** | External classification | RNA-seq | GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760 |
| **GSE39582** | External stage + survival | Microarray | GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582 |
| **GSE17536** | External survival | Microarray | GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17536 |

No raw controlled-access data are used. The Saudi molecular data are used **only as a
published-literature pathway comparison** (Alsolme et al. 2023, *Diagnostics*; Alatwi
et al. 2025, *Front. Oncol.*) — never as training or validation input.

> **Note on ICGC:** the protocol originally planned ICGC COAD-READ for survival
> validation, but the ICGC data portal was retired (June 2024) and now requires DACO
> approval. It was replaced with the fully open-access GEO cohorts above.

## 4. Quick start

```bash
git clone https://github.com/SAH-ML/crc-stage-stratified-xai.git
cd crc-stage-stratified-xai
python -m venv venv && source venv/bin/activate     # (Windows: venv\Scripts\activate)
pip install -r requirements.txt
```

Then follow **docs/REPRODUCE.md** step by step. Each `code/NN_*/` folder has its own
README and is runnable independently; `results/` already contains the expected outputs
so you can verify your run matches.

## 5. Pipeline overview

| Step | Folder | What it does |
|------|--------|--------------|
| 1 | `code/01_preprocessing` | Download/QC TCGA-COAD; build expression matrix + labels |
| 2 | `code/02_differential_expression` | DE (Python engine) + canonical **R limma** |
| 2b | `code/02b_limma_concordance` | Verify Python vs R limma agree (Jaccard 0.995) |
| 3 | `code/03_network_enrichment` | STRING/MCODE network; 4-platform triangulated enrichment; two-level analysis; Saudi cross-ref |
| 4 | `code/04_model_benchmark` | Six-algorithm benchmark; bootstrap CIs; class-imbalance analysis |
| 5 | `code/05_stage_shap` | Stage-stratified SHAP; cross-stage concordance; external replication |
| 6 | `code/06_biological_corroboration` | Literature + GEPIA2 corroboration of the panel |
| 7 | `code/07_survival` | Cox / DeepSurv / RSF survival benchmark (the honest null) |
| 8 | `code/08_external_validation` | Three-tier external validation |

## 6. Reproducibility guarantees

- **Random seed 42** throughout.
- **Leakage-free**: all data-dependent steps fit on training folds only (see paper §2.3.4).
- Every step's `results/` folder holds the exact CSV/JSON the paper reports, so you can
  diff your run against ours.
- Canonical DE numbers come from **R limma** (Bioconductor); the Python engine is a
  validated cross-check (concordance in `results/02b_limma_concordance/`).

## 7. Citation

If you use this code or data, please cite the paper (full citation to be added on
publication) and this repository.

## 8. License

MIT — see [LICENSE](LICENSE). You are free to use, modify, and distribute with attribution.
