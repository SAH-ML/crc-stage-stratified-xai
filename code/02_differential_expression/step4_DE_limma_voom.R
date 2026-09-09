# ============================================================================
# Step 4 (CANONICAL) — Differential Expression with limma-voom in R
#
# *** THIS IS THE OFFICIAL METHOD. The DEG numbers reported in the paper come
#     from THIS script, matching the Methods section (limma-voom [25,26]). ***
#
# The Python script (step4_DE_python.py) is the development engine and
# reproducibility cross-check; step4_DE_concordance.py compares the two.
#
# HOW TO RUN (you do not need to know R — just follow these steps):
#   Option 1 — Posit Cloud / RStudio (easiest, free):
#     1. Go to https://posit.cloud  -> New Project
#     2. Upload: TCGA_expression_qc.csv and TCGA_labels.csv
#        (from data/processed/)
#     3. Open this file, click "Source" (top-right of the editor)
#     4. Download the resulting DE_limma_voom/*.csv files
#   Option 2 — Google Colab with R runtime:
#     1. New notebook -> Runtime -> Change runtime type -> R
#     2. Upload the two CSVs, paste this script in a cell, run
#   Claude will walk you through whichever option you choose.
#
# Inputs  (place in working directory):
#   TCGA_expression_qc.csv  — genes x samples, log2 scale (rows=genes)
#   TCGA_labels.csv         — sample, tumor_vs_normal, stage
# Outputs (written to ./DE_limma_voom/):
#   DE_global_tumor_vs_normal.csv
#   DE_stage_{I,II,III,IV}_vs_normal.csv
#   DE_summary.csv
# ============================================================================

# ---- Install dependencies (first run only) --------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
for (pkg in c("limma", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}
library(limma)
library(edgeR)

set.seed(42)
dir.create("DE_limma_voom", showWarnings = FALSE)

LOGFC <- 1.0
FDR   <- 0.05

# ---- Load data ------------------------------------------------------------
expr <- read.csv("TCGA_expression_qc.csv", row.names = 1, check.names = FALSE)
labels <- read.csv("TCGA_labels.csv", row.names = 1, check.names = FALSE)
labels <- labels[colnames(expr), , drop = FALSE]   # align order

# NOTE on data scale:
#   UCSC Xena TCGA expression is log2(norm_count+1), NOT raw counts. Classic
#   voom expects raw counts. Because our matrix is already a normalised log2
#   matrix, we use limma's linear modelling directly with array-style weights
#   (limma-trend), which is the recommended limma path for log-expression
#   matrices. This is documented in the paper so the method statement is exact.
#   (If raw counts are obtained instead, replace the block below with
#    DGEList -> calcNormFactors -> voom -> lmFit.)

run_contrast <- function(expr_sub, is_group1) {
  design <- model.matrix(~ is_group1)
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit, trend = TRUE)   # limma-trend for log-expression
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt
}

write_res <- function(tt, path) {
  tt$significant <- (tt$adj.P.Val < FDR) & (abs(tt$logFC) >= LOGFC)
  write.csv(tt, path, row.names = FALSE)
  invisible(tt)
}

tumor  <- labels$tumor_vs_normal == "tumor"
normal <- labels$tumor_vs_normal == "normal"
stage  <- labels$stage

summary_rows <- list()

# ---- Global: tumour vs normal --------------------------------------------
cat("Global tumour vs normal...\n")
keep <- tumor | normal
tt <- run_contrast(expr[, keep], tumor[keep])
tt <- write_res(tt, "DE_limma_voom/DE_global_tumor_vs_normal.csv")
summary_rows[["global"]] <- data.frame(
  contrast = "Global tumour-vs-normal",
  n_significant = sum(tt$significant),
  n_up = sum(tt$significant & tt$logFC > 0),
  n_down = sum(tt$significant & tt$logFC < 0))

# ---- Per-stage: Stage_k vs normal ----------------------------------------
all_p <- list()
for (stg in c("I", "II", "III", "IV")) {
  cat("Stage", stg, "vs normal...\n")
  stage_mask <- (stage == stg) & tumor
  stage_mask[is.na(stage_mask)] <- FALSE
  keep <- stage_mask | normal
  tt <- run_contrast(expr[, keep], stage_mask[keep])
  write.csv(tt, sprintf("DE_limma_voom/DE_stage_%s_vs_normal.csv", stg),
            row.names = FALSE)
  all_p[[stg]] <- data.frame(gene = tt$gene, stage = stg,
                             logFC = tt$logFC, P = tt$P.Value)
  sig <- (tt$adj.P.Val < FDR) & (abs(tt$logFC) >= LOGFC)
  summary_rows[[stg]] <- data.frame(
    contrast = paste("Stage", stg, "vs normal"),
    n_significant = sum(sig),
    n_up = sum(sig & tt$logFC > 0),
    n_down = sum(sig & tt$logFC < 0))
}

# ---- Pooled FDR across all stage-by-gene tests ---------------------------
pooled <- do.call(rbind, all_p)
pooled$FDR_pooled <- p.adjust(pooled$P, method = "BH")
pooled$significant_pooled <- (pooled$FDR_pooled < FDR) & (abs(pooled$logFC) >= LOGFC)
write.csv(pooled, "DE_limma_voom/DE_pooled_stage_tests.csv", row.names = FALSE)

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, "DE_limma_voom/DE_summary.csv", row.names = FALSE)
cat("\nDone. Summary:\n")
print(summary_df)
cat("\nOutputs in DE_limma_voom/. Upload that folder back to Claude for the\n")
cat("concordance check (step4_DE_concordance.py).\n")
