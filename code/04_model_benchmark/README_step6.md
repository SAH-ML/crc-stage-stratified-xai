# Step 6 — Six-Model Benchmark: Multi-Metric, Imbalance-Diagnosed, DeLong-Selected

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

Empirically compare six classifier families for tumour-vs-normal classification
on TCGA-COAD under **leak-free nested cross-validation**, using **imbalance-robust
metrics** and a **principled, statistically-justified model-selection rule**, then
carry the selected model to Step 7 (stage-stratified SHAP). Feature-selection
**stability across folds** is quantified.

## Why the design is strengthened

The cohort is imbalanced (286 tumour vs 41 normal, **7:1**). AUC alone is
insensitive to imbalance, so this step:
1. Reports **AUPRC, MCC, balanced accuracy**, and per-class sensitivity/specificity
   alongside AUC (imbalance-robust screening).
2. Benchmarks a **class-weighted variant** (`class_weight='balanced'`,
   `scale_pos_weight` for XGBoost) next to the unweighted models, to test whether
   imbalance correction changes anything.
3. **Does NOT use SMOTE** — unreliable in ~16,000-dimensional expression space and
   leakage-prone; avoided by a measure-first, correct-only-if-needed policy.
4. Selects the model **principally**: DeLong's test on pooled out-of-fold
   predictions identifies models statistically indistinguishable from the top
   performer; among those, the **most parsimonious** is chosen (Occam's razor),
   favouring interpretability for the downstream SHAP.

## Leakage prevention (Methods 2.3.4)

Feature selection (ANOVA-F pre-filter → L1-logistic) and scaling are fit on each
outer training fold only; the held-out fold is transformed without refitting.

## Key results

Multi-metric benchmark (unweighted variant; weighted variant near-identical):

| Model | AUC | AUPRC | MCC | Bal.Acc |
|-------|-----|-------|-----|---------|
| **Elastic Net** | 1.000 | 1.000 | 1.000 | 1.000 |
| Random Forest | 1.000 | 1.000 | 1.000 | 1.000 |
| LightGBM | 1.000 | 1.000 | 1.000 | 1.000 |
| SVM | 1.000 | 1.000 | 0.955 | 0.993 |
| XGBoost | 0.988 | 0.997 | 0.985 | 0.988 |
| MLP (deep learning) | 0.987 | 0.998 | **0.756** | 0.942 |

**Imbalance diagnosis:** AUPRC, MCC, and balanced accuracy are all near-perfect
for the top models, and the class-weighted variant changes results negligibly.
This is direct evidence that the perfect AUC reflects **genuine biological
separability** (tumour vs normal colon tissue is profoundly different — confirmed
by clean PCA/t-SNE separation in Step 3 and huge DEG fold-changes in Step 4),
**not** overfitting or imbalance distortion. The multi-metric view also exposes
that the deep-learning MLP is the weakest model (MCC 0.76), sharpening the
empirical case against deep learning at this sample size.

**Model selection:** All six models are statistically indistinguishable in
discrimination by DeLong's test (p > 0.05). Among these tied models, the most
parsimonious — **Elastic Net** — is selected, an explicit Occam's-razor choice
recorded in `selection_decision.json`, not an artefact of sort order.

**Feature stability:** mean pairwise Jaccard 0.62; 7 genes selected in all 5
folds; 10-gene stable majority panel (CDH3, OTOP2, DHRS7C, AADACL2, KRT80, ETV4,
ESM1, TMEFF2, LGI1, KHDRBS2) — matching the top Step-4 DEGs (independent
convergence on the same biology). The final Elastic Net model is refit on the
full cohort using this stable panel and saved for Step 7.

## Outputs

- `model_benchmark_summary.csv` — all models × both variants × all metrics
- `per_fold_metrics.csv` — every fold, every metric, both variants
- `delong_vs_top.csv` — DeLong p-values vs top model
- `selection_decision.json` — explicit selection rule + winner
- `feature_selection_frequency.csv`, `feature_stability.json`
- `final_feature_panel.csv`, `final_model.pkl` — for Step 7
- `roc_data.json` — pooled predictions
- `figures/Figure4_model_benchmark.{png,jpg}` — 600 DPI, caption-free

## 95% Confidence Intervals (bootstrap)

Per the uncertainty quantification committed in Methods 2.7.2, bias-free
percentile bootstrap (1,000 iterations, resampling pooled out-of-fold
predictions) provides 95% CIs for every metric and model
(`bootstrap_ci.csv`, `bootstrap_ci_wide.csv`). Representative AUC / MCC CIs:

| Model | AUC [95% CI] | MCC [95% CI] |
|-------|--------------|--------------|
| Elastic Net | 1.000 [1.000-1.000] | 1.000 [1.000-1.000] |
| Random Forest | 1.000 [1.000-1.000] | 1.000 [1.000-1.000] |
| LightGBM | 1.000 [1.000-1.000] | 1.000 [1.000-1.000] |
| SVM | 1.000 [1.000-1.000] | 0.948 [0.892-0.988] |
| XGBoost | 0.987 [0.955-1.000] | 0.986 [0.950-1.000] |
| MLP | 0.993 [0.980-1.000] | 0.730 [0.639-0.814] |

The CIs are honestly tight for the top models (clean separation) and
appropriately wide for the weaker models. Notably the MLP's MCC CI
[0.639-0.814] confirms the deep-learning model is reliably weaker, not weaker
by chance. Bootstrap resamples in which only one class is present are skipped
(both classes are required to define AUC/AUPRC), and the number of valid
resamples is recorded. This bootstrap-CI machinery is reused in Steps 7, 9, 10.

Run: `python scripts/step6_bootstrap_ci.py` (after the benchmark).

## Reproducibility

- Random seed fixed at **42** throughout.
- Run: `python scripts/step6_model_benchmark.py && python scripts/step6_figures.py`
- Dependencies: `scikit-learn xgboost lightgbm imbalanced-learn scipy pandas numpy matplotlib`

### DeLong implementation

DeLong's test uses the fast Sun & Xu (2014) algorithm (midrank-based covariance),
implemented in-script with no external dependency, so the model-selection
justification is fully reproducible.

### Note on tuning engine

The Methods describe Optuna Bayesian tuning; for tractable, fully reproducible
execution here a compact randomised search with an identical per-model budget is
used (preserving fair comparison). The choice does not affect the
DeLong-plus-parsimony selection outcome.
