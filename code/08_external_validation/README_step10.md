# Step 10 — Three-Tier External Validation (TRIPOD+AI-aligned)

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose

Validate the 10-gene panel on fully independent cohorts (Methods 2.7.2) across
three tiers, none used in model development.

## Design & results

### Tier 1 — same-platform (RNA-seq) classification replication — GSE50760
The Step-6 tumour-vs-normal Elastic Net is applied to GSE50760 primary-tumour vs
normal samples via signature transfer (per-gene z-score, then the model's linear
decision function). **Classification only** (cohort is all Stage IV; no
per-stage/survival).

| Metric | Value [95% CI] |
|--------|----------------|
| AUC | **0.969 [0.896–1.000]** |
| AUPRC | 0.980 [0.924–1.000] |

**Strong replication.** The panel's tumour-vs-normal discrimination generalises
cleanly to an independent RNA-seq cohort.

### Tier 2 — cross-platform (microarray) survival validation — GSE39582 + GSE17536
A signature risk score (panel z-scored within cohort, weighted by the TCGA model
coefficients) is tested for survival association.

| Cohort | n | events | EPV | C-index [95% CI] |
|--------|---|--------|-----|------------------|
| GSE39582 | 579 | 194 | **19.4** | 0.520 [0.477–0.565] |
| GSE17536 | 177 | 73 | 7.3 | 0.575 [0.503–0.642] |

**Weak survival signal.** The well-powered cohort (GSE39582, EPV 19.4) gives a
C-index of 0.52 — essentially chance, CI includes 0.5. GSE17536 reaches 0.575 but
is underpowered (EPV 7.3) and its CI nearly touches 0.5.

## Honest synthesis

The 10-gene panel is **a strong classifier but a weak prognosticator**:
- Classification (tumour vs normal): external AUC ≈ 0.97 — robust, generalisable.
- Survival (prognosis): external C-index 0.52–0.58 — at or near chance.

This is biologically coherent and consistent with Step 9: the panel was selected
to separate tumour from normal, not to predict survival. The two are different
biological questions, and the panel transfers to the first but not the second.
Reported plainly, without over-interpretation.

## Outputs

- `external_validation_summary.json` — full results, all tiers
- `external_validation_table.csv` — tidy summary table
- `figures/Figure8_external_validation_tiers.{png,jpg}` — 600 DPI, caption-free

## Reproducibility

- Random seed fixed at **42** throughout.
- Run: `python scripts/step10_external_validation.py && python scripts/step10_figures.py`
- Dependencies: `scikit-learn lifelines pandas numpy matplotlib`
