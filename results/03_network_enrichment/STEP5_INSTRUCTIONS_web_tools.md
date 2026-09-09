# STEP 5 — PPI Network + Triangulated Enrichment: COMPLETE WEB-TOOL INSTRUCTIONS

This guide walks you through every web step for Step 5, matching the paper's
Methods (Sections 2.4.2, 2.4.3, 2.4.4) EXACTLY. Follow it in order. Each step
says precisely what to click, what settings to use, and what file to save.

**You will use these input files (in this folder):**
- `step5_DEG_STRINGcore_symbols.txt`  → 1,618 gene symbols, ONE PER LINE → for STRING
- `step5_DEG_full_symbols.txt`        → 4,807 gene symbols, ONE PER LINE → for DAVID/KEGG/Reactome
- `step5_DEG_full_for_enrichment.csv` → same 4,807 genes WITH log2FC + FDR (backup / reference)

**Why two lists (important — read once):**
The Methods send the full DEG list to enrichment and the DEGs to STRING at
confidence > 0.9. The full list (4,807) is correct for enrichment tools, which
handle large lists. STRING's web interface, however, becomes unstable above
~2,000 genes, so for the STRING network we use a stricter, web-tractable DEG core
(|log2FC| ≥ 2, still FDR < 0.05 → 1,618 genes), and the confidence > 0.9 filter is
then applied *inside* STRING exactly as the Methods state. All 10 SHAP-panel genes
are present in BOTH lists. This split is documented and will be reported.

Total time: ~45–60 minutes. Take your time; save every file with the EXACT name
given so I can match them when you send them back.

===============================================================================
## PART A — STRING PPI NETWORK + MCODE  (Methods 2.4.2)
===============================================================================

STRING confidence > 0.9, then MCODE hub detection
(degree cutoff=2, node score cutoff=0.2, k-core=2, max depth=100).

You can do MCODE in TWO ways. **Option A1 (Cytoscape)** matches the Methods most
exactly (MCODE via Cytoscape). **Option A2 (STRING website only)** is the fallback
if you cannot install Cytoscape. Do A1 if you can; otherwise A2.

-------------------------------------------------------------------------------
### OPTION A1 — STRING + Cytoscape + MCODE (PREFERRED, matches Methods exactly)
-------------------------------------------------------------------------------

**A1.1 — Get the high-confidence network from STRING**
1. Open https://string-db.org/
2. Top menu → click **"Multiple proteins"** (NOT "Single protein").
3. Open `step5_DEG_STRINGcore_symbols.txt`, select ALL (Ctrl+A), copy (Ctrl+C).
4. Paste into the big "List Of Names" box.
5. **Organism: type and select "Homo sapiens".**
6. Click **"SEARCH"**.
7. If a "resolve identifiers" page appears (some symbols map to multiple),
   click **"Continue"** to accept the best matches. (A few unmapped genes is normal.)
8. The network loads. Now click the **"Settings"** tab (below the network).
   - Set **"minimum required interaction score"** to **"highest confidence (0.900)"**.
     ← THIS IS THE confidence > 0.9 STEP FROM THE METHODS. Do not skip.
   - Leave other settings default.
   - Click **"UPDATE"** at the bottom of Settings.
9. Click the **"Exports"** tab.
   - Download **"as short tabular text output"** → save as
     **`STRING_network_core_conf0.9.tsv`**
   - Also click **"as image"** (high-res PNG) → save as
     **`STRING_network_image.png`** (for the supplementary figure).

**A1.2 — Open the network in Cytoscape and run MCODE**
10. Install Cytoscape (free): https://cytoscape.org/download.html (if not already).
11. EASIEST path to get the network into Cytoscape:
    - Back on the STRING result page, click the **"Exports"** tab →
      look for **"Send network to Cytoscape"** (requires Cytoscape OPEN with the
      stringApp installed). If you don't have stringApp: in Cytoscape go to
      **Apps → App Manager → search "stringApp" → Install**, then retry.
    - Alternatively: in Cytoscape, **File → Import → Network from File** and
      select the `STRING_network_core_conf0.9.tsv` you saved. Map column
      "node1" = Source, "node2" = Target.
12. Install MCODE in Cytoscape: **Apps → App Manager → search "MCODE" → Install**.
13. Run MCODE: **Apps → MCODE → Open MCODE**. In the MCODE panel set EXACTLY:
    - Network scoring: **Degree Cutoff = 2**
    - Cluster finding: **Node Score Cutoff = 0.2**
    - **K-Core = 2**
    - **Max. Depth = 100**
    - (Haircut: leave checked/default; Fluff: leave unchecked/default)
    - Click **"Analyze"** (select the whole network when prompted).
14. MCODE lists clusters ranked by score. For EACH cluster (at least the top 3),
    click it to select its nodes, then:
    - **File → Export → Table to File** (or right-click the selected nodes →
      export), OR simply note the member genes.
15. Save the MCODE results:
    - Export the node table (with the "MCODE_Cluster" / "MCODE_Score" columns) →
      **`MCODE_clusters.csv`**
    - Save a screenshot/image of the top cluster(s) →
      **`MCODE_clusters_image.png`**

-------------------------------------------------------------------------------
### OPTION A2 — STRING website only (FALLBACK if no Cytoscape)
-------------------------------------------------------------------------------
1–9. Do steps A1.1 (1–9) above exactly (gives the conf-0.9 network + tsv + image).
10. On the STRING result page, click the **"Clusters"** tab.
11. Choose **"MCL clustering"** (STRING's built-in option) — note: this is MCL,
    not MCODE. Set inflation parameter to default (e.g., 3). Click **UPDATE**.
12. Export the clustered network: **"Exports"** tab →
    "as short tabular text output" → save as **`STRING_MCL_clusters.tsv`**.
13. ⚠ Tell me you used Option A2 (MCL, website) so I report the clustering method
    accurately. (We will note this deviation from MCODE transparently; the
    Methods can be footnoted to state MCL-on-STRING was used if Cytoscape was
    unavailable. PREFER Option A1 if at all possible.)

===============================================================================
## PART B — TRIANGULATED ENRICHMENT (4 platforms)  (Methods 2.4.3)
===============================================================================

Submit the **FULL DEG list** (`step5_DEG_full_symbols.txt`, 4,807 genes) to FOUR
platforms. A pathway counts as "robust" only if FDR < 0.05 in **≥ 2 of the 4**.
Save each platform's full results table.

-------------------------------------------------------------------------------
### B1 — DAVID
-------------------------------------------------------------------------------
1. Open https://david.ncifcrf.gov/
2. Click **"Start Analysis"** (top menu).
3. On the left "Upload" panel:
   - Step 1: paste the contents of `step5_DEG_full_symbols.txt` into the box
     (or use "Choose File" to upload the .txt).
   - Step 2: **Select Identifier = "OFFICIAL_GENE_SYMBOL"**.
   - Step 3: **List Type = "Gene List"**.
   - Step 4: click **"Submit List"**.
4. If asked to select species, choose **"Homo sapiens"**.
5. Once loaded, click **"Functional Annotation Tool"**.
6. Open these annotation categories and click **"Chart"** for each:
   - **KEGG_PATHWAY**
   - **GOTERM_BP_DIRECT** (biological process)
   - (optional: Reactome within DAVID if shown)
7. In each Chart view, click **"Download File"** (top) to get the full table.
   Save as:
   - **`DAVID_KEGG.txt`**
   - **`DAVID_GO_BP.txt`**
   (These tables include the FDR/Benjamini column — keep it.)

-------------------------------------------------------------------------------
### B2 — KEGG (via Enrichr, the standard accessible route)
-------------------------------------------------------------------------------
KEGG's own site has no direct multi-gene enrichment upload; the standard
accessible route used in papers is Enrichr (which includes KEGG libraries).
1. Open https://maayanlab.cloud/Enrichr/
2. Paste the contents of `step5_DEG_full_symbols.txt` into the input box.
3. Click **"Submit"**.
4. After it loads, find the **"Pathways"** section and click **"KEGG 2021 Human"**.
5. Click the **"Table"** view, then **"Export"** (download the table).
   Save as **`KEGG_Enrichr.txt`**.
   (Enrichr reports adjusted p-value / FDR — keep that column.)

-------------------------------------------------------------------------------
### B3 — Reactome
-------------------------------------------------------------------------------
1. Open https://reactome.org/PathwayBrowser/#TOOL=AT
   (or https://reactome.org/ → "Analyze gene list / Tools → Analyse gene list").
2. Click **"Analyse gene list"**.
3. Paste the contents of `step5_DEG_full_symbols.txt`.
4. **Tick "Project to human"**. Click **"Analyse"**.
5. When results load, click **"Download → Results (CSV)"** (or the download icon).
   Save as **`Reactome_results.csv`**.
   (Reactome reports an "Entities FDR" column — keep it.)

-------------------------------------------------------------------------------
### B4 — ClueGO (Cytoscape plugin)
-------------------------------------------------------------------------------
If you installed Cytoscape for Part A, do ClueGO too; if not, this one can be
SKIPPED and we proceed with DAVID+KEGG+Reactome (still ≥2-of-N triangulation —
tell me if you skip it so I adjust the "≥2 of 4" rule to "≥2 of 3" transparently).
1. In Cytoscape: **Apps → App Manager → search "ClueGO" → Install**
   (ClueGO needs a free license key — register at the prompt; academic email).
2. **Apps → ClueGO → Start ClueGO**.
3. Number of clusters = 1. Paste `step5_DEG_full_symbols.txt` genes into Cluster 1.
4. Organism: **Homo sapiens**. Identifier type: **SymbolID**.
5. Ontologies: tick **GO Biological Process** and **KEGG** and **Reactome**.
6. Set **"Use GO Term Fusion"** on; significance **pV ≤ 0.05** with
   **Bonferroni step-down** correction.
7. Click the **"Start"** (play) button.
8. Export the results table: ClueGO results panel → **export table** →
   save as **`ClueGO_results.csv`**, and save the network image as
   **`ClueGO_image.png`**.

===============================================================================
## PART C — GEPIA2 (Step 8b, batched here)  — see separate file
===============================================================================
Also complete `../step8/GEPIA2_instructions_for_user.md` while you are in
browser-tool mode (≈10 min). Send those results back too.

===============================================================================
## WHAT TO SEND BACK TO ME
===============================================================================
Please send back whichever of these you produced (exact names help me a lot):
- STRING:   `STRING_network_core_conf0.9.tsv`, `STRING_network_image.png`
- MCODE:    `MCODE_clusters.csv`, `MCODE_clusters_image.png`
            (or `STRING_MCL_clusters.tsv` if you used Option A2)
- DAVID:    `DAVID_KEGG.txt`, `DAVID_GO_BP.txt`
- KEGG:     `KEGG_Enrichr.txt`
- Reactome: `Reactome_results.csv`
- ClueGO:   `ClueGO_results.csv`, `ClueGO_image.png`  (if done)
- GEPIA2:   your filled table + saved box plots

Then I will:
1. Build the triangulated enrichment table (pathway robust if FDR<0.05 in ≥2 tools).
2. Locate the 10 SHAP-panel genes within the STRING network (degree, betweenness,
   MCODE-module membership) — the two-level analysis (Methods 2.4.4).
3. Complete Step 8c (Saudi pathway cross-reference) using the enriched pathways.
4. Produce the Step 5 figure (600 DPI, caption-free) and README.

===============================================================================
## COMMON PITFALLS (read to avoid mistakes)
===============================================================================
- Always set organism = **Homo sapiens** everywhere.
- Always use **official gene symbols** (our files already are symbols).
- In STRING, DON'T forget to set confidence to **0.900 (highest)** and click UPDATE.
- For MCODE use the EXACT parameters: degree=2, node score=0.2, k-core=2, depth=100.
- Keep the **FDR / adjusted-p / Benjamini** column in every downloaded table —
  that is what the ≥2-of-4 rule is based on.
- If any tool refuses 4,807 genes, tell me; we can submit the STRING-core 1,618
  list to that tool instead and I will note it.
- Save files with the EXACT names above; if you must rename, tell me which is which.
