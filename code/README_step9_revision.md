# Step 9 — Revision Analyses (Peer-Review Response)

Confirmatory analyses added during peer review of the manuscript
*"Stage-Stratified Explainable Machine Learning Identifies a Reproducible
Transcriptomic Classifier with Stage-Dependent Feature Importance in Colorectal
Cancer: A Multi-Cohort Study with Saudi Population Contextualisation"*
(Biomedicines, manuscript ID biomedicines-4590563).

These scripts reproduce the additional Supplementary results requested by the
reviewers. All use random seed 42 and the same ten-gene panel and leakage-free
conventions as the main pipeline.

## Scripts

| Script | Reviewer comment | Produces |
|--------|------------------|----------|
| `step9a_permutation_repeated_cv.py` | R1.6, R2.4 | `permutation_repeated_cv.json` → Supplementary Table S10 |
| `step9b_stage_bootstrap_ci.py` | R1.1, R2.2 | `tableS8_stage_bootstrap_ci.csv` → Supplementary Table S8 |
| `step9c_esm1_adjustment.py` | R1.7, R2.5 | `esm1_adjustment.json` → Supplementary Table S11 |
| `step9d_permutation_figure.py` | R1.6 | `FigureS4_permutation_null.png` (600 dpi) → Supplementary Figure S4 |

## How to run

From the repository root (so the relative `data/` and `results/` paths resolve):

```bash
python code/09_revision_analyses/step9a_permutation_repeated_cv.py
python code/09_revision_analyses/step9b_stage_bootstrap_ci.py
python code/09_revision_analyses/step9c_esm1_adjustment.py
python code/09_revision_analyses/step9d_permutation_figure.py
```

Outputs are written to `results/09_revision_analyses/`.

## Key results

- **Permutation test:** observed classifier AUC = 1.000; permuted-label null mean
  = 0.50 (maximum 0.69); permutation p < 0.001. The result is robust to the
  Elastic Net hyperparameters. Repeated nested CV (10 seeds): mean AUC 1.000,
  SD 0.000. → the near-perfect internal AUC reflects the genuine, large
  tumour-versus-normal transcriptomic distance, not optimism or overfitting.
- **Stage-wise bootstrap CIs:** point estimates match the main text; the
  monotonic gradient is preserved (Stage I vs IV weakest/widest; Stage III vs IV
  strong and tightly bounded).
- **ESM1 adjustment:** the ESM1–stage association retains significance after
  partial-correlation adjustment for the other nine panel genes (partial
  rho = 0.120, p = 0.045); effect size is modest, consistent with the cautious
  interpretation in the manuscript.

## Dependencies

Same as the main pipeline (`requirements.txt`): numpy, pandas, scipy,
scikit-learn, shap, matplotlib. Seed 42 throughout.

License: MIT.
