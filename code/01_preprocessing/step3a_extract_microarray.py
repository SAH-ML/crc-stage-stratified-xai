#!/usr/bin/env python3
"""
Step 3a — Extract expression matrices from GEO SOFT files (microarray cohorts).

GSE39582 and GSE17536 are Affymetrix GPL570 microarray series whose expression
tables are embedded inside the family.soft.gz file (one table per sample). This
script parses them into gene-level expression matrices and maps Affymetrix probe
IDs to gene symbols, so all cohorts share a common gene-symbol index.

Part of the reproducible pipeline for:
  "Stage-Stratified Explainable AI Reveals Shifting Transcriptomic Drivers
   of Colorectal Cancer Progression"

Author: [Author list]
License: MIT
"""

import gzip
import re
from pathlib import Path
import pandas as pd
import numpy as np

RAW = Path("data/raw")
PROC = Path("data/processed")
PROC.mkdir(parents=True, exist_ok=True)


def parse_soft_expression(soft_path: Path) -> pd.DataFrame:
    """
    Parse a GEO family.soft.gz file into a probe x sample expression matrix.

    Memory-efficient: collects each sample as a list of (probe, value) pairs,
    converts to a float32 Series keyed by probe, and concatenates once at the
    end. This avoids holding a giant dict-of-dicts in memory for large series
    (e.g., GSE39582 with 585 samples).

    Each sample block is delimited by:
        ^SAMPLE = GSMxxxxx
        !sample_table_begin
        ID_REF<TAB>VALUE
        <probe><TAB><value>
        !sample_table_end
    """
    sample_series = {}
    current_gsm = None
    in_table = False
    header_seen = False
    probes = []
    values = []

    def _flush():
        """Store the current sample's collected probe/value pairs as a Series."""
        if current_gsm is not None and probes:
            s = pd.Series(values, index=probes, dtype="float32")
            # collapse duplicate probe ids within a sample (rare) by mean
            s = s.groupby(level=0).mean()
            sample_series[current_gsm] = s

    with gzip.open(soft_path, "rt") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                _flush()
                current_gsm = line.split("=")[1].strip()
                probes, values = [], []
                in_table = False
                header_seen = False
            elif line.startswith("!sample_table_begin"):
                in_table = True
                header_seen = False
            elif line.startswith("!sample_table_end"):
                in_table = False
            elif in_table:
                if not header_seen:
                    header_seen = True  # skip header row (ID_REF  VALUE)
                    continue
                tab = line.find("\t")
                if tab == -1:
                    continue
                probe = line[:tab]
                rest = line[tab + 1:]
                # value is up to the next tab or newline
                tab2 = rest.find("\t")
                vstr = rest[:tab2] if tab2 != -1 else rest.rstrip("\n")
                try:
                    values.append(float(vstr))
                    probes.append(probe)
                except ValueError:
                    pass
    _flush()  # last sample

    expr = pd.DataFrame(sample_series)
    return expr


def parse_platform_annotation(soft_path: Path) -> dict:
    """
    Extract probe -> gene symbol mapping from the platform (^PLATFORM) table
    embedded in the SOFT file. For GPL570, the gene symbol column is
    'Gene Symbol'.
    """
    probe_to_gene = {}
    in_platform_table = False
    header = None
    symbol_idx = None
    id_idx = None

    with gzip.open(soft_path, "rt") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_platform_table = True
                header = None
                continue
            if line.startswith("!platform_table_end"):
                in_platform_table = False
                continue
            if in_platform_table:
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = parts
                    # find the ID and Gene Symbol columns
                    for i, col in enumerate(header):
                        if col.strip().upper() == "ID":
                            id_idx = i
                        if "GENE SYMBOL" in col.strip().upper():
                            symbol_idx = i
                    continue
                if id_idx is not None and symbol_idx is not None and len(parts) > symbol_idx:
                    probe = parts[id_idx]
                    symbol = parts[symbol_idx].strip()
                    # GPL570 sometimes lists multiple symbols "A /// B"; take first
                    if symbol and symbol != "":
                        symbol = symbol.split("///")[0].strip()
                        probe_to_gene[probe] = symbol
    return probe_to_gene


def collapse_probes_to_genes(expr: pd.DataFrame, probe_to_gene: dict) -> pd.DataFrame:
    """
    Map probes to gene symbols and collapse multiple probes per gene by taking
    the probe with the highest mean expression (standard 'maxMean' approach).
    """
    expr = expr.copy()
    expr["gene"] = expr.index.map(probe_to_gene)
    expr = expr.dropna(subset=["gene"])
    expr = expr[expr["gene"] != ""]

    # maxMean collapse: for each gene keep the probe row with the highest mean
    expr["_mean"] = expr.drop(columns=["gene"]).mean(axis=1, numeric_only=True)
    expr = expr.sort_values("_mean", ascending=False)
    expr = expr.drop_duplicates(subset="gene", keep="first")
    expr = expr.set_index("gene").drop(columns=["_mean"])
    return expr


def main():
    for gse in ["GSE39582", "GSE17536"]:
        soft = RAW / f"{gse}_family.soft.gz"
        print(f"\n=== {gse} ===")
        print("  Parsing expression tables ...")
        expr = parse_soft_expression(soft)
        print(f"  Raw probe matrix: {expr.shape[0]} probes x {expr.shape[1]} samples")

        print("  Parsing platform annotation (probe -> gene) ...")
        mapping = parse_platform_annotation(soft)
        print(f"  Probe-to-gene mappings: {len(mapping)}")

        print("  Collapsing probes to genes (maxMean) ...")
        gene_expr = collapse_probes_to_genes(expr, mapping)
        print(f"  Gene-level matrix: {gene_expr.shape[0]} genes x {gene_expr.shape[1]} samples")

        out = PROC / f"{gse}_expression_gene.csv"
        gene_expr.to_csv(out)
        print(f"  Saved -> {out}")


if __name__ == "__main__":
    main()
