# Step 7b — External Validation of the Stage-Stratified Pattern (GSE39582)

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

The critical replication test: does the stage-dependent importance structure
found in TCGA-COAD reproduce in an **independent cohort on a different platform**
(GSE39582, Affymetrix microarray) with its own TNM staging?

## Method — exact linear-SHAP transfer (Methods 2.2.3)

For a linear model (the Step-6 Elastic Net winner), the SHAP contribution of
gene *g* is **exactly** `coef_g * (z_g - E[z_g])`. We therefore transfer the
10-gene panel and compute attributions directly:

1. Transfer the panel to GSE39582; per-gene z-score within GSE39582 (places the
   microarray platform on a comparable, mean-centred scale).
2. Per-sample signed contribution `= coef_g * z_g`, using the TCGA model's
   learned coefficients.
3. Stage importance `= mean|contribution|` within each TNM stage subgroup,
   mirroring the TCGA `mean|SHAP|` computation exactly.

**This exact method supersedes an earlier pseudo-label recalibration approach.**
The diagnostic in `step7c_transfer_diagnostic.py` showed the pseudo-label method
*understated* cross-cohort agreement; the exact method (M1) and a
variance-weighted check (M2) agree at ρ=0.97, confirming robustness. Using the
mathematically exact transfer is a methodological correction, not result
selection.

## Result — partial replication (honest)

GSE39582 usable tumours with TNM stage: **579** (I=38, II=271, III=210, IV=60).
All four stages pass the bootstrap stability bar.

**Cross-cohort per-stage importance concordance (TCGA vs GSE39582):**

| Stage | Spearman ρ | p | Replicates |
|-------|-----------|---|------------|
| I | 0.83 | 0.003 | yes |
| II | 0.52 | 0.13 | no (n.s.) |
| III | 0.94 | <0.001 | yes (strong) |
| IV | 0.82 | 0.004 | yes |

**3 of 4 stages replicate strongly and significantly.** Stage II is the
consistent weak spot across transfer methods.

**Key gene findings:**
- **ESM1 (angiogenesis) replicates cleanly** — rises to rank 1 at Stage IV in
  GSE39582 (rank 2 in TCGA). Late-stage angiogenesis prominence is robust across
  cohorts and platforms (Figure 6B).
- **CDH3 decline is weaker** externally than in TCGA (cohort-dependent).
- **ETV4 rise does NOT replicate** — flat in GSE39582; this was TCGA-specific.

## Honest interpretation

The **overall stage-stratified importance structure is strongly conserved** at
3 of 4 stages, and **ESM1's late-stage emergence is a reproducible, cross-platform
finding**. However, the *specific* CDH3→ETV4 trajectory seen in TCGA is partly
cohort-dependent and does not fully replicate. The defensible cross-cohort claim
is therefore: *stage-dependent redistribution of importance is real and its
structure largely replicates, with ESM1 emerging as a reproducibly
stage-dependent driver; specific individual-gene trajectories require
cohort-aware interpretation.*

## Outputs

- `external_validation_gse39582.json` — concordances, stability, trajectories
- `gse39582_shap_importance_by_stage.csv` — importance per gene × stage
- `transfer_diagnostic.json` — M1/M2/M3 comparison (artefact check)
- `figures/Figure6_external_validation.{png,jpg}` — 600 DPI, caption-free

## Reproducibility

- Random seed fixed at **42** throughout.
- Run: `python scripts/step7b_external_stage_validation.py && python scripts/step7b_figures.py`
- Diagnostic: `python scripts/step7c_transfer_diagnostic.py`
- Dependencies: `scikit-learn scipy pandas numpy matplotlib`.
