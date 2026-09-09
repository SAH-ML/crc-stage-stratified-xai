# Step 8b — GEPIA2 Expression-Direction Confirmation (USER web step)

GEPIA2 (http://gepia2.cancer-pku.cn/) is a browser tool the analysis sandbox
cannot reach, so this step must be run by you in a browser. It independently
confirms the tumour-vs-normal expression direction of each panel gene in an
independent TCGA+GTEx-based reference, complementing our own DE analysis.

## What to do (≈10 minutes)

For EACH of the 10 panel genes below, do the following on GEPIA2:

1. Go to http://gepia2.cancer-pku.cn/#analysis
2. Enter the gene symbol in the "Gene" box.
3. Under "Expression DIY" → "Box Plot":
   - Dataset: select **COAD** (and optionally READ)
   - Log Scale: tick (log2 TPM)
   - Match TCGA normal AND GTEx normal: leave default (adds normal reference)
   - Jitter Size: default
4. Click "Plot". Note whether the **tumour (red) median is higher or lower than
   normal (grey)**, and whether the difference is marked significant (*).
5. Record the direction in the table below.

## Genes and our expected direction (from our DE analysis)

| Gene | Our direction (tumour vs normal) | GEPIA2 direction (you fill) | Matches? |
|------|----------------------------------|------------------------------|----------|
| CDH3 | UP | | |
| KRT80 | UP | | |
| ETV4 | UP | | |
| ESM1 | UP | | |
| OTOP2 | DOWN | | |
| DHRS7C | DOWN | | |
| AADACL2 | DOWN | | |
| LGI1 | DOWN | | |
| TMEFF2 | DOWN | | |
| KHDRBS2 | DOWN | | |

## Optional but useful (stage correlation for ESM1)

For **ESM1** specifically (our key replicated late-stage gene), also run:
- "Expression DIY" → "Stage Plot" → Dataset COAD → Plot.
- This shows ESM1 expression across pathologic stages I–IV. Note whether it
  trends upward with stage (our finding predicts higher in late stage).

## What happens next

Save the GEPIA2 box plots (right-click → save image) and paste your filled
table back to me. I will:
- tabulate concordance between GEPIA2 and our DE directions,
- add it to the Supplementary material and reference it in the Discussion
  (NOT the abstract), reporting the outcome regardless of direction.

This is a confirmatory cross-check; even partial concordance is reported honestly.
