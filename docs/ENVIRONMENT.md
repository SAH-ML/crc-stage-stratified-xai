# Environment & Reproducibility Notes

- **Python:** 3.10-3.12. Create an isolated venv and `pip install -r requirements.txt`.
- **R:** 4.x with Bioconductor `limma` + `edgeR` (canonical DE only, Step 2). Posit Cloud
  (https://posit.cloud) is the easiest way to run the R step if R is unfamiliar.
- **Random seed:** 42, fixed in every stochastic step (model training, bootstrap, SHAP
  sampling) for bit-reproducible results.
- **Determinism:** bootstrap CIs use 1,000 iterations; betweenness centrality is exact.
- **Hardware:** runs on a standard laptop; no GPU required (DeepSurv trains on CPU in minutes).
- **Web tools** (Steps 3 & 6): STRING, Cytoscape/MCODE, DAVID, KEGG/Enrichr, Reactome,
  ClueGO, GEPIA2 — run interactively; their exported outputs are committed under
  `results/` so the downstream Python is fully reproducible without re-running them.

## Expected run time (full pipeline, CPU laptop)
| Step | Approx time |
|------|-------------|
| 1 Preprocessing | 2-3 min |
| 2 DE (Python) | 1 min |
| 2 DE (R limma) | ~12 min (first-time Bioconductor install) |
| 3 Enrichment (Python part) | 1 min |
| 4 Model benchmark + bootstrap | 5-8 min |
| 5 Stage SHAP + external | 5-8 min |
| 7 Survival benchmark | 3-5 min |
| 8 External validation | 2 min |
