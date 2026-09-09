# Step 9 — Survival Benchmark on the 10-Gene SHAP Panel (TCGA-COAD)

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

Benchmark three survival architectures (Methods 2.7.1) — Cox Proportional Hazards,
Random Survival Forest (RSF), and DeepSurv — on the 10-gene SHAP panel, using
TCGA-COAD overall survival, to assess whether the panel carries prognostic signal.

## Events-per-variable (EPV) honesty — pre-registered (Methods 2.7.1, Table 6 Issue 3)

The 10-gene panel needs EPV ≥ 10 (≥ 100 events) for a confirmatory model.
Actual TCGA OS events:

| Subgroup | n | events | EPV | Designation |
|----------|---|--------|-----|-------------|
| Stage I | 45 | 3 | 0.3 | exploratory |
| Stage II | 110 | 22 | 2.2 | exploratory |
| Stage III | 80 | 21 | 2.1 | exploratory |
| Stage IV | 39 | 19 | 1.9 | exploratory |
| **Overall (pooled)** | 283 | 69 | **6.9** | primary (EPV caveat) |

**Every per-stage subgroup falls below the threshold**, so all per-stage survival
results are labelled **exploratory**, never confirmatory. The overall pooled model
(EPV 6.9) is also below the conventional threshold and is reported with an explicit
caveat. This was computed and written out first (`epv_table.csv`); nothing is hidden.

## Result — the panel does NOT carry a strong survival signal (honest)

| Model | C-index [95% CI] | 3-yr AUC | 5-yr AUC |
|-------|------------------|----------|----------|
| Cox PH | 0.525 [0.450–0.606] | 0.589 | 0.482 |
| DeepSurv | 0.506 [0.426–0.586] | 0.520 | 0.504 |
| RSF | 0.504 [0.425–0.584] | 0.502 | 0.465 |

**A C-index of 0.5 is chance-level. All three models hover at 0.50–0.52, and
every 95% CI includes 0.5** — none demonstrates statistically significant survival
discrimination. The Kaplan-Meier risk-group separation is not significant
(log-rank p = 0.376).

## Interpretation (honest)

This result is biologically coherent: the 10-gene panel was selected to
distinguish **tumour from normal** (Step 6), which is a different biological
question from **predicting survival time**. Genes that separate cancer from
healthy tissue are not necessarily prognostic. The panel was optimised for
classification, not prognosis, and the survival benchmark confirms it does not
transfer to the survival task. This is reported plainly rather than over-interpreted.

**Implication for the paper's framing:** survival is not a strength of this panel
and cannot serve as a load-bearing pillar of the paper. The contribution rests on
the stage-stratified SHAP methodology and the externally-replicated importance
structure (Step 7), not on survival prediction.

## Outputs

- `epv_table.csv` — events-per-variable per subgroup (the honesty artefact)
- `survival_benchmark_overall.csv` — C-index + CIs + time-dependent AUC per model
- `survival_perstage_exploratory.csv` — per-stage C-index (exploratory, labelled)
- `risk_scores_overall.csv` — out-of-fold risk scores + median-split risk groups
- `survival_summary.json` — summary + interpretation
- `figures/Figure7_survival_benchmark.{png,jpg}` — 600 DPI, caption-free

## Reproducibility

- Random seed fixed at **42** throughout.
- Run: `python scripts/step9_survival_benchmark.py && python scripts/step9_figures.py`
- Dependencies: `lifelines scikit-survival pycox torch torchtuples scikit-learn pandas numpy matplotlib`
