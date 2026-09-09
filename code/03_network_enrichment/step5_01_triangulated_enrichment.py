#!/usr/bin/env python3
"""
Step 5 (script 1/3) - Triangulated enrichment across four platforms (Methods 2.4.3).

Reads the six enrichment outputs produced by the web tools (DAVID across its
KEGG / GO-BP / Reactome ontologies; Enrichr-KEGG; standalone Reactome; ClueGO),
applies an FDR<0.05 significance filter within each, maps them onto the FOUR
methodology platforms (DAVID, KEGG, Reactome, ClueGO), and flags a pathway as
"robust" if it is significant in >=2 of the 4 platforms.

Platform-mapping rationale (documented, faithful, no double-counting):
  * DAVID    = union of DAVID's own KEGG + GO_BP + Reactome ontology results
               (DAVID as a single tool).
  * KEGG     = Enrichr "KEGG 2021 Human" (KEGG database via an enrichment tool).
  * Reactome = DAVID-Reactome UNION standalone reactome.org. The standalone
               reactome.org applies an ultra-conservative whole-hierarchy FDR
               (here min FDR ~0.68 -> 0 terms at FDR<0.05), whereas DAVID's
               Reactome ontology (same database, standard Benjamini-Hochberg)
               returns 71. Reactome is credited where significant in EITHER
               rendering. This is why multi-platform triangulation is used: no
               single tool/correction is treated as ground truth.
  * ClueGO   = ClueGO (GO-BP + KEGG + Reactome, Bonferroni step-down).
  The same underlying database is never counted as two independent platforms,
  so the >=2-of-4 robustness rule is not inflated.

Inputs (web-tool outputs, in --uploads dir):
  DAVID_KEGG.csv, DAVID_GO_BP.csv, DAVID_Reactome.csv  (cols: Term, FDR, User Ids)
  Enrichr_for_KEGG_2026_table.txt  (cols: Term, Adjusted P-value, Genes)
  Reactome_results.csv             (cols: Pathway name, Entities FDR, Submitted entities found)
  ClueGOResultTable-1.xls          (cols: Term, 'Term PValue Corrected with Bonferroni step down', Associated Genes Found)

Outputs (results/step5/):
  triangulated_enrichment_all.csv     - every pathway, # platforms, which platforms
  triangulated_enrichment_robust.csv  - pathways significant in >=2 of 4

License: MIT   |   Random seed: 42 (no stochastic step here)
"""
import argparse, re
from pathlib import Path
import pandas as pd

FDR_THRESHOLD = 0.05


def norm(name: str) -> str:
    """Normalise a pathway/term name for cross-platform matching."""
    n = str(name).lower().strip()
    n = re.sub(r'\br-hsa-\d+\b', '', n)        # strip Reactome IDs
    n = re.sub(r'[^a-z0-9 ]', ' ', n)
    n = re.sub(r'\bhomo sapiens\b', '', n)
    n = re.sub(r'\bpathway\b', '', n)
    n = re.sub(r'\s+', ' ', n).strip()
    return n


def sig_terms(path, term_col, fdr_col, sep=','):
    """Return {normalised_term: original_term} for rows with FDR < threshold."""
    if str(path).endswith(('.csv', '.txt')):
        d = pd.read_csv(path, sep=sep)
    else:
        d = pd.read_excel(path)
    d = d[pd.to_numeric(d[fdr_col], errors='coerce') < FDR_THRESHOLD]
    return {norm(t): t for t in d[term_col]}


def main(uploads, outdir):
    uploads, outdir = Path(uploads), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    david_kegg  = sig_terms(uploads / 'DAVID_KEGG.csv',     'Term', 'FDR')
    david_gobp  = sig_terms(uploads / 'DAVID_GO_BP.csv',    'Term', 'FDR')
    david_react = sig_terms(uploads / 'DAVID_Reactome.csv', 'Term', 'FDR')
    enrichr     = sig_terms(uploads / 'Enrichr_for_KEGG_2026_table.txt',
                            'Term', 'Adjusted P-value', sep='\t')
    react_std   = sig_terms(uploads / 'Reactome_results.csv',
                            'Pathway name', 'Entities FDR')
    clue        = sig_terms(uploads / 'ClueGOResultTable-1.xls',
                            'Term', 'Term PValue Corrected with Bonferroni step down')

    four = {
        'DAVID':    set(david_kegg) | set(david_gobp) | set(david_react),
        'KEGG':     set(enrichr),
        'Reactome': set(david_react) | set(react_std),
        'ClueGO':   set(clue),
    }
    print("Significant terms per methodology platform (FDR<0.05):")
    for p, s in four.items():
        print(f"  {p}: {len(s)}")

    # human-readable label per normalised term
    label = {}
    for src in [david_kegg, david_gobp, david_react, enrichr, react_std, clue]:
        for nt, orig in src.items():
            if nt and nt not in label:
                label[nt] = orig

    rows = []
    for nt, lab in label.items():
        hits = [p for p, s in four.items() if nt in s]
        rows.append({'pathway': lab, 'n_platforms': len(hits),
                     'platforms': ';'.join(hits), 'robust_ge2': len(hits) >= 2})
    res = pd.DataFrame(rows).sort_values(['n_platforms', 'pathway'],
                                         ascending=[False, True])
    robust = res[res['robust_ge2']]
    res.to_csv(outdir / 'triangulated_enrichment_all.csv', index=False)
    robust.to_csv(outdir / 'triangulated_enrichment_robust.csv', index=False)

    print(f"\nTotal significant (>=1 platform): {len(res)}")
    print(f"ROBUST (>=2 of 4 platforms):      {len(robust)}")
    print(f"Significant in >=3 platforms:     {(res['n_platforms'] >= 3).sum()}")
    print(f"Saved to {outdir}/")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--uploads', default='/mnt/user-data/uploads',
                    help='directory containing the six web-tool output files')
    ap.add_argument('--outdir', default='results/step5')
    a = ap.parse_args()
    main(a.uploads, a.outdir)
