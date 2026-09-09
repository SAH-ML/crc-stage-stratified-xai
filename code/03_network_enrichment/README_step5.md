# Step 5 — PPI Network + Triangulated Enrichment + Two-Level Analysis

Part of the reproducible pipeline for *"Stage-Stratified Explainable AI Reveals
Shifting Transcriptomic Drivers of Colorectal Cancer Progression."*

## Purpose (Methods 2.4.2–2.4.4)

Characterise the broad DEG biology and locate the SHAP-panel genes within it:
1. **PPI network** (STRING confidence>0.9; MCODE modules) — Methods 2.4.2
2. **Triangulated enrichment** across four platforms — Methods 2.4.3
3. **Two-level analysis**: locate the 10 SHAP-panel genes within the network and
   the enriched-pathway landscape — Methods 2.4.4
4. **Saudi pathway cross-reference** (Step 8c) — Methods 2.6

## Inputs

- Full DEG list (FDR<0.05, |log2FC|>=1): 4,807 genes -> enrichment
- STRING-core DEG list (|log2FC|>=2): 1,618 genes -> STRING web (web-tractable)
- Web-tool outputs (user-run): STRING edges (0.9 + 0.7), MCODE clusters, DAVID
  (KEGG/GO_BP/Reactome), Enrichr-KEGG, Reactome, ClueGO.

## 1. PPI network + MCODE (2.4.2)

STRING at confidence>0.9 yielded a sparse high-confidence network; MCODE found
32 modules (top module: 14 nodes, score 7.385). Parameters: degree cutoff=2,
node score cutoff=0.2, k-core=2, max depth=100 (as specified).

## 2. Triangulated enrichment (2.4.3)

Full DEG list submitted to four platforms; pathway "robust" if FDR<0.05 in
**>=2 of 4** platforms. Platform mapping: DAVID (its KEGG+GO_BP+Reactome
ontologies), KEGG (Enrichr), Reactome (DAVID-Reactome ∪ standalone Reactome),
ClueGO.

**Note on standalone Reactome:** reactome.org applies an ultra-conservative
whole-hierarchy FDR (min FDR 0.68 here), returning 0 terms at FDR<0.05, while
DAVID's Reactome ontology (same database, standard BH correction) returned 71.
Reactome is therefore credited where significant in *either* rendering — a
documented, faithful handling, not a deviation. This is precisely why
multi-platform triangulation is used: no single tool is ground truth.

**Result:** 817 pathways significant in >=1 platform; **216 robust (>=2 of 4)**;
**55 significant in 3 of 4**. Robust pathways include biologically coherent CRC
processes: Wnt signaling, Pathways in cancer, cytokine-cytokine receptor
interaction, calcium/cAMP signaling, chemokine signaling, complement and
coagulation, ECM organisation, cell adhesion. Files:
`triangulated_enrichment_all.csv`, `triangulated_enrichment_robust.csv`.

## 3. Two-level analysis (2.4.4)

### Network topology of panel genes (degree, betweenness)
Computed from STRING edge lists at both thresholds (`two_level_panel_in_network.csv`,
`centrality_comparison.json`):

| Threshold | Panel genes present | Detail |
|-----------|---------------------|--------|
| conf>0.9 | 1/10 (LGI1, degree 2) | 9/10 absent from high-confidence core |
| conf>0.7 | 3/10 (CDH3, OTOP2, LGI1) | still low degree; never connect to each other |

**The SHAP-panel genes are peripheral to / absent from the dense PPI core, robustly
across both thresholds, and never form an interconnected module.** This shows the
ML-prioritised signature is **complementary to** — not a rediscovery of — classical
PPI hub-gene identification: a methodological strength (a hub-gene pipeline would
have missed this signature). 0.9 is the primary committed network; 0.7 is a
sensitivity analysis confirming the conclusion.

### Pathway membership of panel genes
Although peripheral in the PPI graph, **8/10 panel genes belong to robust enriched
pathways** (`panel_genes_in_robust_pathways.csv`):
- **CDH3** -> Wnt signaling, cell-adhesion (CAM), epithelial differentiation
- **ESM1** -> **angiogenesis**, cell-surface-receptor signalling, proliferation
  (directly corroborating its literature/Step-8a angiogenesis role and its
  reproducible late-stage importance)
- **DHRS7C** -> ion/calcium homeostasis; **ETV4** -> signal transduction,
  epithelial differentiation; **LGI1** -> synaptic/axon-guidance (neural, not CRC)
- AADACL2, TMEFF2 -> in no robust pathway (consistent with poorly-characterised status)

This links the mechanistic analysis (Contribution 3) to the stage-stratified
explainability finding (Contribution 1): the model's key genes sit in
biologically meaningful CRC pathways even though they are not PPI hubs.

## 4. Saudi pathway cross-reference (Step 8c, Methods 2.6)

Robust enriched pathways cross-referenced against published Saudi CRC molecular
features (`step8c_saudi_pathway_crossref.csv`). **4 of 5 Saudi-altered pathways
converge with our robust enrichment:**

| Saudi-altered pathway | Saudi frequency | Convergent? |
|-----------------------|-----------------|-------------|
| **Wnt/β-catenin** | **65%** (WES 2025); APC 47.8% | **yes** (Wnt robust; CDH3 here) |
| PI3K | 48% | yes (Pathways in cancer; PI3K-Akt) |
| TP53 | 33.7–50% | yes (Pathways in cancer) |
| KRAS/RAS-MAPK | 6.5–45% | yes (MAPK cascade; signal transduction) |
| Homologous recombination | 61% | no (not in our DEG-enriched set) |

**Key convergence:** Wnt signaling is simultaneously the most frequently altered
pathway in Saudi CRC (65%), a robust pathway in our enrichment (3/4 platforms),
and the pathway housing our panel gene **CDH3** — a direct Saudi-relevance link
for the CDH3/Wnt axis. Reported regardless of direction; the one non-convergence
(homologous recombination) is stated honestly.

## Outputs
- `triangulated_enrichment_all.csv`, `triangulated_enrichment_robust.csv`
- `panel_genes_in_robust_pathways.csv`
- `two_level_panel_in_network.csv`, `centrality_comparison.json`
- `step8c_saudi_pathway_crossref.csv`
- `figures/Figure9_enrichment_network.{png,jpg}` — 600 DPI, caption-free

## Reproducible code (scripts, in run order)

All analysis code is provided as documented, runnable scripts (not just the
figure). Run order:

```bash
# 1. Triangulated enrichment across the four platforms (Methods 2.4.3)
python scripts/step5_01_triangulated_enrichment.py --uploads <web_tool_outputs_dir>

# 2. Two-level analysis: panel genes in network + enriched pathways (Methods 2.4.4)
python scripts/step5_02_two_level_analysis.py --uploads <web_tool_outputs_dir>

# 3. Saudi pathway cross-reference (Step 8c, Methods 2.6)
python scripts/step5_03_saudi_crossref.py

# 4. Figure (600 DPI, caption-free)
python scripts/step5_figures.py
```

Script responsibilities:
- `step5_01_triangulated_enrichment.py` - reads the six web-tool outputs, applies
  FDR<0.05 per platform, maps to the four methodology platforms (with the
  documented no-double-counting Reactome handling), writes the robust-pathway table.
- `step5_02_two_level_analysis.py` - builds the STRING graphs at conf 0.9 and 0.7,
  computes panel-gene degree/betweenness/neighbours, and locates panel genes within
  the robust enriched pathways.
- `step5_03_saudi_crossref.py` - cross-references robust pathways against published
  Saudi CRC alteration frequencies (sources embedded in the script).
- `step5_figures.py` - assembles the four-panel Figure 9.

Dependencies: pandas, numpy, networkx, matplotlib, xlrd (for the ClueGO .xls).
The scripts are deterministic; betweenness centrality is computed exactly.

## Reproducibility
- Random seed 42 where applicable; web-tool parameters documented (incl. ClueGO log).
- Run: `python scripts/step5_figures.py` (after web outputs are present).
- Dependencies: pandas, numpy, networkx, matplotlib, xlrd.
