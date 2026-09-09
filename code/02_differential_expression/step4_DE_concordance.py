#!/usr/bin/env python3
"""
Step 4 (concordance) — Compare Python and canonical R limma-voom DE results.

Run this in PHASE 2, after the R limma-voom script has produced its outputs
(results/DE_limma_voom/*.csv) alongside the Python outputs
(results/DE_python/*.csv).

It quantifies and reports, for the Supplementary Material:
  * Overlap of significant-DEG sets (Jaccard, % shared of top-N).
  * Spearman correlation of gene rankings (by p-value / by log2FC).
  * A concordance summary table.

This is the trust anchor demonstrating that the two independent DE
implementations agree, justifying use of the Python pipeline as a faithful
reproducibility engine while reporting canonical limma-voom numbers.

License: MIT
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PY = Path("results/DE_python")
R = Path("results/DE_limma_voom")
OUT = Path("results")


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else np.nan


def compare_contrast(py_file, r_file, py_sig_col, name):
    py = pd.read_csv(PY / py_file)
    r = pd.read_csv(R / r_file)

    # Standardise column names between implementations
    # Python: gene, log2FC, p_value, FDR_per_contrast, <sig col>
    # R (limma topTable): gene, logFC, P.Value, adj.P.Val, significant
    r = r.rename(columns={"logFC": "log2FC", "P.Value": "p_value",
                          "adj.P.Val": "FDR_per_contrast"})

    merged = py.merge(r, on="gene", suffixes=("_py", "_r"))

    # ranking correlation
    rho_p, _ = spearmanr(merged["p_value_py"], merged["p_value_r"])
    rho_fc, _ = spearmanr(merged["log2FC_py"], merged["log2FC_r"])

    # significant set overlap
    py_sig = set(py.loc[py[py_sig_col], "gene"])
    r_sig = set(r.loc[r["significant"], "gene"])
    jac = jaccard(py_sig, r_sig)
    shared = len(py_sig & r_sig)
    pct_py = shared / len(py_sig) * 100 if py_sig else np.nan
    pct_r = shared / len(r_sig) * 100 if r_sig else np.nan

    return {
        "contrast": name,
        "n_sig_python": len(py_sig),
        "n_sig_R_limma": len(r_sig),
        "n_shared": shared,
        "jaccard": round(jac, 3),
        "pct_of_python_shared": round(pct_py, 1),
        "pct_of_R_shared": round(pct_r, 1),
        "spearman_pvalue_rank": round(rho_p, 3),
        "spearman_log2FC": round(rho_fc, 3),
    }


def main():
    rows = []
    rows.append(compare_contrast(
        "DE_global_tumor_vs_normal.csv", "DE_global_tumor_vs_normal.csv",
        "significant", "Global tumour-vs-normal"))
    for s in ["I", "II", "III", "IV"]:
        rows.append(compare_contrast(
            f"DE_stage_{s}_vs_normal.csv", f"DE_stage_{s}_vs_normal.csv",
            "significant_per_contrast", f"Stage {s} vs normal"))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "DE_concordance_report.csv", index=False)
    print("=== DE concordance: Python vs canonical R limma-voom ===")
    print(df.to_string(index=False))
    print(f"\nWritten to {OUT/'DE_concordance_report.csv'}")
    print("\nInterpretation: high Jaccard / high % shared / Spearman near 1.0")
    print("indicate the two implementations agree; the paper reports the R")
    print("limma-voom numbers, with this table as Supplementary evidence of")
    print("reproducibility across implementations.")


if __name__ == "__main__":
    main()
