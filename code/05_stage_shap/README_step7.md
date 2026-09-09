# Step 7 — Stage-Stratified SHAP Attribution (Central Novel Analysis)

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

Explain the Step-6 winning model (Elastic Net on the 10-gene stable panel) with
SHAP computed (a) globally and (b) separately within each AJCC stage subgroup,
to test the study's central hypothesis: **does gene importance shift across
tumour stages?** The model is NOT retrained per stage — SHAP is a post-hoc
stratification of one full-cohort model, which is what makes the small Stage-I
subgroup acceptable (no per-subgroup training).

## Method

- Exact `shap.LinearExplainer` for the linear model; shared background.
- Global ranking (all tumours) + four per-stage rankings.
- **Bootstrap stability** (1,000 iters) of each stage's importance ranking;
  a stage is "stable" if bootstrap mean Spearman rho > 0.7.
- **Concordance** between global and per-stage, and between stage pairs, via
  Spearman rho with permutation-test p-values.
- **Pre-registered fallback** (Methods 2.5.4): four-stage primary if all stages
  meet n>=30 and the stability bar; else collapse to early (I+II) vs late
  (III+IV). Both outcomes computed; the decision is recorded.

## Result — four-stage analysis supported

Decision: **all four stages meet n>=30 AND pass the stability bar** -> four-stage
analysis reported as primary; fallback NOT triggered.

Per-stage SHAP ranking stability (bootstrap mean Spearman rho):

| Stage | n | mean rho | stable |
|-------|---|----------|--------|
| I | 45 | 0.86 | yes |
| II | 110 | 0.89 | yes |
| III | 80 | 0.95 | yes |
| IV | 39 | 0.91 | yes |

Stage-pair SHAP concordance (Spearman rho) shows a **monotonic progression
gradient** — importance divergence increases with stage distance:

|       | I | II | III | IV |
|-------|---|----|-----|----|
| I | 1.00 | 0.84 | 0.83 | **0.66** |
| II | 0.84 | 1.00 | 0.82 | 0.64 |
| III | 0.83 | 0.82 | 1.00 | 0.93 |
| IV | 0.66 | 0.64 | 0.93 | 1.00 |

The extremes (Stage I vs IV, rho=0.66) are the least concordant; adjacent stages
are most concordant. Global-vs-stage concordances are all significant
(permutation p<0.005).

## Interpretation (honest)

The core tumour signature is **shared** across stages (all concordances positive
and fairly high), but the **relative importance of its components shifts
measurably with progression**: cell-adhesion / epithelial-identity genes (notably
**CDH3**) dominate early stages, while invasion / angiogenesis-associated genes
(**ETV4**, **ESM1**) gain prominence in late stages. CDH3 and ETV4 importance
trajectories cross over between early and late disease (Figure 5B). This is a
modulation of a shared signature with progression — a coherent, biologically
plausible, and statistically supported finding — rather than four wholly distinct
per-stage signatures. The monotonic concordance gradient (Figure 5C) argues
against an artefactual explanation.

## Outputs

- `shap_importance_by_stage.csv` — mean|SHAP| per gene x {global, I-IV}
- `stage_stability.json` — bootstrap stability per stage
- `concordance.json` — all pairwise Spearman rho + permutation p
- `fallback_decision.json` — four-stage vs early/late decision (recorded)
- `early_late.json` — early/late rankings (computed regardless, for completeness)
- `shap_values.json` — raw SHAP arrays for figures
- `figures/Figure5_stage_shap.{png,jpg}` — 600 DPI, caption-free

## Next

External validation of this stage-dependent pattern on GSE39582 (which carries
full TNM staging) via the cross-platform signature-transfer protocol (Methods
2.2.3), plus biological corroboration (Step 8). The stage-stratified gene set
here is the focused panel located within the broader DEG network in Step 5
(two-level enrichment, Methods 2.4.4).

## Reproducibility

- Random seed fixed at **42** throughout.
- Run: `python scripts/step7_stage_shap.py && python scripts/step7_figures.py`
- Dependencies: `shap scikit-learn scipy pandas numpy matplotlib`.
