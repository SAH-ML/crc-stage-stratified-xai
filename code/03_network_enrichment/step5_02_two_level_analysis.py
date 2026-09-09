#!/usr/bin/env python3
"""
Step 5 (script 2/3) - Two-level network analysis (Methods 2.4.4).

Locates the 10 SHAP-prioritised panel genes within (a) the STRING PPI network and
(b) the robust enriched-pathway landscape.

Part A - network topology:
  For each STRING confidence threshold (0.9 primary; 0.7 sensitivity), build the
  PPI graph from the exported edge list and compute, for each panel gene:
  presence, degree, betweenness centrality, and neighbours. This quantifies how
  peripheral/central the ML-prioritised genes are in the PPI core. (Primary =
  0.9 as committed in Methods 2.4.2; 0.7 reported as a sensitivity analysis.)

Part B - pathway membership:
  Using the per-pathway gene memberships parsed from the enrichment outputs,
  identify which ROBUST enriched pathways (from script 1) contain each panel gene.

Inputs:
  STRING_network_core_conf0_9.tsv, STRING_network_core_conf0_7.tsv
      (STRING "short tabular text output": cols #node1, node2, ...)
  the six enrichment files (for pathway->gene membership)
  results/step5/triangulated_enrichment_robust.csv (from script 1)

Outputs (results/step5/):
  two_level_panel_in_network.csv        - per gene x threshold: degree, betweenness
  centrality_comparison.json            - summary of 0.9 vs 0.7
  panel_genes_in_robust_pathways.csv    - per gene: # and names of robust pathways

License: MIT   |   Random seed: 42 (betweenness is deterministic here)
"""
import argparse, json, re
from pathlib import Path
import pandas as pd
import networkx as nx

PANEL = ['CDH3', 'OTOP2', 'DHRS7C', 'AADACL2', 'KRT80',
         'ETV4', 'ESM1', 'TMEFF2', 'LGI1', 'KHDRBS2']


def norm(name: str) -> str:
    n = str(name).lower().strip()
    n = re.sub(r'\br-hsa-\d+\b', '', n)
    n = re.sub(r'[^a-z0-9 ]', ' ', n)
    n = re.sub(r'\bhomo sapiens\b', '', n)
    n = re.sub(r'\bpathway\b', '', n)
    return re.sub(r'\s+', ' ', n).strip()


def network_topology(uploads, outdir):
    rows = []
    summary = {}
    for conf, fn in [("0.9", "STRING_network_core_conf0_9.tsv"),
                     ("0.7", "STRING_network_core_conf0_7.tsv")]:
        path = uploads / fn
        if not path.exists():
            print(f"  [skip] {fn} not found"); continue
        edges = pd.read_csv(path, sep='\t')
        G = nx.from_pandas_edgelist(edges, '#node1', 'node2')
        degs = dict(G.degree())
        betw = nx.betweenness_centrality(G)          # deterministic (exact)
        for g in PANEL:
            if g in G:
                rows.append({"confidence": conf, "gene": g, "in_network": "yes",
                             "degree": degs[g], "betweenness": round(betw[g], 5),
                             "n_neighbors": G.degree(g),
                             "neighbors": ";".join(list(G.neighbors(g)))})
            else:
                rows.append({"confidence": conf, "gene": g, "in_network": "no",
                             "degree": 0, "betweenness": 0,
                             "n_neighbors": 0, "neighbors": ""})
        summary[f"conf_{conf}"] = {
            "nodes": G.number_of_nodes(), "edges": G.number_of_edges(),
            "panel_present": sum(1 for g in PANEL if g in G)}
    pd.DataFrame(rows).to_csv(outdir / "two_level_panel_in_network.csv", index=False)
    json.dump(summary, open(outdir / "centrality_comparison.json", "w"), indent=2)
    return summary


def pathway_membership(uploads, outdir):
    # build pathway(normalised) -> set(genes) from all enrichment files
    pathway_genes = {}

    def add(path, term_col, gene_col, sep=','):
        if str(path).endswith(('.csv', '.txt')):
            d = pd.read_csv(path, sep=sep)
        else:
            d = pd.read_excel(path)
        for _, r in d.iterrows():
            nt = norm(r[term_col])
            genes = str(r[gene_col]).replace(';', ',').replace('|', ',').split(',')
            genes = [x.strip() for x in genes if x.strip()]
            pathway_genes.setdefault(nt, set()).update(genes)

    add(uploads / 'DAVID_KEGG.csv', 'Term', 'User Ids')
    add(uploads / 'DAVID_GO_BP.csv', 'Term', 'User Ids')
    add(uploads / 'DAVID_Reactome.csv', 'Term', 'User Ids')
    add(uploads / 'Enrichr_for_KEGG_2026_table.txt', 'Term', 'Genes', sep='\t')
    add(uploads / 'Reactome_results.csv', 'Pathway name', 'Submitted entities found')
    add(uploads / 'ClueGOResultTable-1.xls', 'Term', 'Associated Genes Found')

    robust = pd.read_csv(outdir / 'triangulated_enrichment_robust.csv')
    rows = []
    for g in PANEL:
        pws = []
        for _, r in robust.iterrows():
            if g in pathway_genes.get(norm(r['pathway']), set()):
                pws.append((r['pathway'], r['n_platforms']))
        top = sorted(pws, key=lambda x: -x[1])[:5]
        rows.append({"gene": g, "n_robust_pathways": len(pws),
                     "top_pathways": "; ".join(f"{p} [{k}/4]" for p, k in top) or "(none)"})
    pd.DataFrame(rows).to_csv(outdir / "panel_genes_in_robust_pathways.csv", index=False)


def main(uploads, outdir):
    uploads, outdir = Path(uploads), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print("Part A - network topology (0.9 primary, 0.7 sensitivity):")
    summ = network_topology(uploads, outdir)
    for k, v in summ.items():
        print(f"  {k}: {v}")
    print("Part B - pathway membership of panel genes:")
    pathway_membership(uploads, outdir)
    print(f"Saved two-level outputs to {outdir}/")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--uploads', default='/mnt/user-data/uploads')
    ap.add_argument('--outdir', default='results/step5')
    a = ap.parse_args()
    main(a.uploads, a.outdir)
