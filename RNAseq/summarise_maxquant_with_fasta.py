#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_maxquant_with_fasta.py

Convert MaxQuant proteinGroups.txt into a gene-level proteomics
evidence table using the FASTA referenced in parameters.txt.

Outputs
-------
TSV gene-level table with:
- prot_present_any
- prot_present_fraction
- prot_n_samples_present
- prot_n_samples_total
- prot_peptide_counts_unique_max
- prot_peptide_counts_razor_unique_max
- prot_sequence_coverage_pct_max

All arguments are named. Output is TSV.
"""

from __future__ import annotations
import argparse
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import numpy as np
import csv
import sys

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--protein_groups_tsv",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--parameters_txt",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
    )
    return parser.parse_args()


def _normalise_gene(value: str) -> str:
    """
    Normalise gene symbol.
    """
    return re.sub(r"\s+", "", str(value).strip()).upper()


def extract_fasta_path(parameters_txt: Path) -> Optional[str]:
    """
    Extract FASTA path from parameters.txt.
    """
    with parameters_txt.open() as f:
        for line in f:
            if "fasta file" in line.lower():
                return line.strip().split("\t")[-1]
    return None


def build_accession_to_gene_map(fasta_path: Path) -> Dict[str, str]:
    """
    Build accession → gene symbol mapping from FASTA.
    """
    mapping = {}

    with fasta_path.open(errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue

            acc_match = re.search(r"\|([A-Z0-9]+)\|", line)
            gn_match = re.search(r"GN=([A-Za-z0-9\-]+)", line)
            gene_symbol_match = re.search(r"Gene_Symbol=([A-Za-z0-9\-]+)", line)

            acc = acc_match.group(1) if acc_match else None

            gene = None
            if gn_match:
                gene = gn_match.group(1)
            elif gene_symbol_match:
                gene = gene_symbol_match.group(1)

            if acc and gene:
                mapping[acc.split("-")[0]] = _normalise_gene(gene)

    return mapping


def extract_gene_from_protein_ids(protein_ids: str, mapping: Dict[str, str]) -> Optional[str]:
    """
    Map Majority protein IDs → gene.
    """
    if protein_ids is None:
        return None

    for token in str(protein_ids).split(";"):
        acc = token.strip().split("-")[0]
        if acc in mapping:
            return mapping[acc]

    return None


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    fasta_path = extract_fasta_path(parameters_txt=args.parameters_txt)

    if fasta_path is None:
        raise RuntimeError("Could not find FASTA path in parameters.txt")

    fasta_path = Path(fasta_path)
    print("Using FASTA:", fasta_path)

    mapping = build_accession_to_gene_map(fasta_path=fasta_path)
    print("Mapped accessions:", len(mapping))

    df = pd.read_csv(args.protein_groups_tsv, sep="\t", dtype=str)

    # Remove contaminants / reverse
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].fillna("") != "+"]

    df["gene_key"] = df["Majority protein IDs"].map(
        lambda x: extract_gene_from_protein_ids(x, mapping)
    )

    df = df[df["gene_key"].notna()].copy()

    peptide_cols = [c for c in df.columns if c.startswith("Peptides ")]
    peptides = df[peptide_cols].apply(pd.to_numeric, errors="coerce")

    present = peptides > 0

    df["prot_n_samples_present"] = present.sum(axis=1)
    df["prot_n_samples_total"] = len(peptide_cols)
    df["prot_present_fraction"] = (
        df["prot_n_samples_present"] / df["prot_n_samples_total"]
    )
    df["prot_present_any"] = df["prot_n_samples_present"] > 0

    df["prot_peptide_counts_unique_max"] = pd.to_numeric(
        df.get("Peptide counts (unique)", np.nan), errors="coerce"
    )

    df["prot_peptide_counts_razor_unique_max"] = pd.to_numeric(
        df.get("Peptide counts (razor+unique)", np.nan), errors="coerce"
    )

    df["prot_sequence_coverage_pct_max"] = pd.to_numeric(
        df.get("Sequence coverage [%]", np.nan), errors="coerce"
    )

    out = df.groupby("gene_key").agg(
        {
            "prot_present_any": "max",
            "prot_present_fraction": "max",
            "prot_n_samples_present": "max",
            "prot_n_samples_total": "max",
            "prot_peptide_counts_unique_max": "max",
            "prot_peptide_counts_razor_unique_max": "max",
            "prot_sequence_coverage_pct_max": "max",
        }
    ).reset_index()

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print("Genes with protein evidence:", int(out["prot_present_any"].sum()))
    print("Wrote:", args.out_tsv)


if __name__ == "__main__":
    main()
