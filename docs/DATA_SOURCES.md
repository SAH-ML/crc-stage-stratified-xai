# Data Sources — Complete Provenance

Every dataset used in this study is public. This document gives the exact source,
access URL, platform, role, and the processing applied. No controlled-access or
private data are used.

---

## 1. TCGA-COAD (training / discovery cohort)

- **What:** The Cancer Genome Atlas — Colon Adenocarcinoma. RNA-seq gene-expression
  with matched clinical annotation including AJCC pathologic stage.
- **Source:** UCSC Xena platform (mirrors TCGA, no data-use agreement required).
  - Browser: https://xenabrowser.net/
  - Download hub: https://tcga.xenahubs.net/download/
- **Exact files used:**
  - Expression: `TCGA.COAD.sampleMap/HiSeqV2` — gene-level RNA-seq, values are
    `log2(norm_count + 1)`, 20,530 genes.
  - Clinical/phenotype: `TCGA.COAD.sampleMap/COAD_clinicalMatrix`.
  - Survival: `survival/COAD_survival.txt` (overall-survival status + time).
- **Cohort after filtering (confirmed in this study):**
  - 286 primary tumours with usable expression + documented AJCC stage:
    Stage I = 45, II = 110, III = 80, IV = 39 (+12 indeterminate, excluded from
    stage-stratified analyses).
  - 41 adjacent solid-tissue normal samples (tumour-vs-normal comparison group).
  - Survival: 545 patients with annotation, 123 death events.
- **Processing:** sub-stages collapsed to parent AJCC stage (IIA/B/C→II; IIIA/B/C→III);
  3 "[Discrepancy]"-annotated samples treated as unstaged. See `code/01_preprocessing`.

## 2. GSE50760 (external classification cohort)

- **What:** RNA-seq, colorectal tumour vs normal (independent of TCGA).
- **Source:** NCBI GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760
- **Role:** Tier-1 external validation of the classifier (external AUC 0.969).
- **Why RNA-seq:** kept on the same platform family as training to avoid cross-platform
  feature collapse in the *classification* test.

## 3. GSE39582 (external stage + survival cohort)

- **What:** Affymetrix microarray, large CRC cohort with full AJCC staging + overall survival.
- **Source:** NCBI GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582
- **Role:** (a) external replication of the stage-stratified SHAP importance structure
  (n=579 with staging); (b) Tier-2 survival validation (EPV 19.4 — adequately powered).

## 4. GSE17536 (external survival cohort)

- **What:** Microarray CRC cohort with overall-survival annotation.
- **Source:** NCBI GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17536
- **Role:** Tier-3 survival validation (EPV 7.3 — underpowered, reported as exploratory).

## 5. Saudi molecular data (literature comparison ONLY — not training/validation)

Used solely as a pathway-level comparison point in the Discussion. No individual-level
or expression data enter the pipeline.
- **Alsolme et al. 2023**, *Diagnostics* 13, 2993. Saudi CRC genomic panel (n=107).
  https://doi.org/10.3390/diagnostics13182993
- **Alatwi et al. 2025**, *Front. Oncol.* 15, 1679528. Saudi CRC whole-exome pathway
  frequencies (Wnt 65%, PI3K 48%, homologous recombination 61%, RTK/RAS 43%).
  https://doi.org/10.3389/fonc.2025.1679528

## 6. Web-tool resources (Step 3 — network & enrichment)

Run interactively; outputs are included under `results/03_network_enrichment/`.
- **STRING** (PPI network): https://string-db.org/  (confidence > 0.9 primary; 0.7 sensitivity)
- **Cytoscape + MCODE** (hub modules): https://cytoscape.org/
- **DAVID** (enrichment): https://david.ncifcrf.gov/
- **KEGG** (via Enrichr): https://maayanlab.cloud/Enrichr/
- **Reactome**: https://reactome.org/
- **ClueGO** (Cytoscape plugin): http://www.ici.upmc.fr/cluego/
- **GEPIA2** (expression corroboration): http://gepia2.cancer-pku.cn/

## 7. ICGC note (transparency)

The original protocol named ICGC COAD-READ for survival validation. The ICGC data
portal was retired in June 2024 and now requires DACO controlled-access approval, so it
was replaced with the fully open-access GEO cohorts (GSE39582, GSE17536) above. This
substitution is documented in the paper (§2.9, Issue 6).
