#!/usr/bin/env python3
"""
Step 5 / Step 8c (script 3/3) - Saudi CRC pathway cross-reference (Methods 2.6).

Cross-references the robust enriched pathways (from script 1) against published
Saudi-specific CRC molecular alteration frequencies. This contextualises the
findings for the Saudi population WITHOUT using Saudi data for training or
validation (it is an independent comparison only). The outcome is reported
regardless of direction; non-convergences are stated honestly.

Saudi CRC molecular features are hard-coded below with their literature sources
(verified via the published studies cited). Each is matched against our robust
pathways by keyword.

Inputs:
  results/step5/triangulated_enrichment_robust.csv (from script 1)
Outputs:
  results/step5/step8c_saudi_pathway_crossref.csv

License: MIT
"""
import argparse
from pathlib import Path
import pandas as pd

# Saudi CRC altered pathways: (label, frequency/detail, source, keyword matchers)
SAUDI = [
    ("WNT/beta-catenin signaling",
     "65% of Saudi CRC (WES 2025); APC mutation 47.8%; CTNNB1 variances",
     "PMC12540182; Frontiers fonc.2026.1792016", ["wnt"]),
    ("PI3K signaling",
     "48% of Saudi CRC (WES 2025)",
     "PMC12540182", ["pi3k", "pathways in cancer"]),
    ("Homologous recombination / DNA repair",
     "61% of Saudi CRC; BRCA2 79%, ATM 76%, CHEK1 78%",
     "PMC12540182; ResearchGate 311550881", ["dna repair", "homologous"]),
    ("TP53 pathway",
     "TP53 mutated 33.7-50% in Saudi CRC",
     "PMC2981835; Frontiers fonc.2026.1792016", ["p53", "apoptosis", "pathways in cancer"]),
    ("KRAS/RAS-MAPK signaling",
     "KRAS 6.5-45.2% (variable across Saudi studies)",
     "PMC9124608; Frontiers fonc.2026.1792016",
     ["mapk", "ras", "pathways in cancer", "signal transduction"]),
]


def main(outdir):
    outdir = Path(outdir)
    robust = pd.read_csv(outdir / 'triangulated_enrichment_robust.csv')
    rows = []
    for label, freq, src, keys in SAUDI:
        found = [p for p in robust['pathway']
                 if any(k in p.lower() for k in keys)]
        rows.append({"saudi_pathway": label, "saudi_frequency": freq, "source": src,
                     "convergent_with_our_enrichment": "yes" if found else "no",
                     "our_matching_pathways": "; ".join(found[:5])})
        print(f"{label}: {'CONVERGENT' if found else 'not convergent'}"
              + (f" -> {found[:3]}" if found else ""))
    pd.DataFrame(rows).to_csv(outdir / "step8c_saudi_pathway_crossref.csv", index=False)
    n_conv = sum(1 for r in rows if r['convergent_with_our_enrichment'] == 'yes')
    print(f"\nConvergent: {n_conv}/{len(SAUDI)} Saudi-altered pathways.")
    print(f"Saved step8c_saudi_pathway_crossref.csv to {outdir}/")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='results/step5')
    a = ap.parse_args()
    main(a.outdir)
