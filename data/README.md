# Data Folder

This folder holds small processed inputs and pointers to the public raw data.

## Included here
- `TCGA_labels.csv` — sample-level tumour/normal + AJCC stage labels for TCGA-COAD,
  aligned to the expression matrix (small, committed directly).

## Not committed here (download from source — see ../docs/DATA_SOURCES.md)
- `TCGA_expression_qc.csv.gz` — QC'd expression matrix (~17 MB). Rebuilt by
  `code/01_preprocessing`, or attach/download as a GitHub **Release asset** (recommended).
- GEO series matrices for GSE50760, GSE39582, GSE17536 — download from GEO.

## Raw data origin (full detail in ../docs/DATA_SOURCES.md)
| Dataset | Source | URL |
|---------|--------|-----|
| TCGA-COAD | UCSC Xena | https://tcga.xenahubs.net/download/ |
| GSE50760 | NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760 |
| GSE39582 | NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582 |
| GSE17536 | NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17536 |

All data are public; no controlled-access or private data are used.
