#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_maxquant_proteingroups_gene_level.py

Convert MaxQuant proteinGroups.txt into a gene-level proteomics evidence TSV.

This script supports mapping protein identifiers to gene symbols using a local
FASTA file (recommended for reproducibility and offline use).

The FASTA referenced in MaxQuant parameters.txt often contains an absolute path
from the original analysis environment (e.g. Windows drive letter). Therefore,
this script allows passing a local FASTA path explicitly.

Mapping logic
-------------
- Extract identifiers from 'Majority protein IDs' (semicolon-separated).
- Build a mapping dictionary from FASTA headers:
  - UniProt style: >sp|P12345|... GN=GENE ...
  - Ensembl pep style: >ENSP... gene:ENSG... gene_symbol:ABC1 ...
  - Other styles: attempt Gene_Symbol=ABC1

Presence logic
--------------
Presence per sample is called using per-sample peptide count columns:
'Peptides <sample>' > 0

Filtering
---------
Removes entries flagged as:
- Reverse
- Potential contaminant
- Only identified by site

Outputs
-------
TSV with per-gene evidence suitable for downstream merges on gene_key.

All arguments are named. Output is TSV (not comma-separated).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise MaxQuant proteinGroups.txt into a gene-level TSV using a local FASTA."
    )
    parser.add_argument(
        "--protein_groups_tsv",
        required=True,
        type=Path,
        help="Path to MaxQuant proteinGroups.txt (tab-delimited).",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Output gene-level TSV path.",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Local FASTA used for identifier->gene mapping (UniProt or Ensembl pep FASTA).",
    )
    parser.add_argument(
        "--parameters_txt",
        required=False,
        type=Path,
        default=None,
        help="Optional MaxQuant parameters.txt (for provenance/logging).",
    )
    parser.add_argument(
        "--presence_threshold",
        required=False,
        default=0.0,
        type=float,
        help="Threshold for per-sample peptide counts to call presence (default: 0).",
    )
    parser.add_argument(
        "--write_per_sample_columns",
        action="store_true",
        help="If set, include per-sample presence columns in the output TSV.",
    )
    return parser.parse_args()


def _normalise_gene(value: str) -> str:
    """
    Normalise a gene symbol for joins.

    Parameters
    ----------
    value : str
        Gene symbol.

    Returns
    -------
    str
        Normalised gene symbol (uppercase, no whitespace).
    """
    return re.sub(r"\s+", "", str(value).strip()).upper()


def _strip_isoform(identifier: str) -> str:
    """
    Strip isoform suffix for identifiers like P12345-2.

    Parameters
    ----------
    identifier : str
        Identifier.

    Returns
    -------
    str
        Base identifier.
    """
    return str(identifier).split("-")[0].upper()


def _safe_read_tsv(path: Path) -> pd.DataFrame:
    """
    Read a TSV safely, allowing very large fields.

    Parameters
    ----------
    path : Path
        Input TSV path.

    Returns
    -------
    pd.DataFrame
        Loaded DataFrame.
    """
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, engine="c")
    except Exception:
        df = pd.read_csv(path, sep="\t", dtype=str, engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _parse_fasta_header(header: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse an identifier and gene symbol from a FASTA header.

    Supports:
    - UniProt: >sp|P12345|... GN=GENE ...
    - Ensembl pep: >ENSP0000... gene_symbol:ABC1 ...
    - Custom: Gene_Symbol=ABC1

    Parameters
    ----------
    header : str
        FASTA header line.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (identifier, gene_symbol) or (None, None) if not parsable.
    """
    h = header.strip()
    if not h.startswith(">"):
        return (None, None)

    ident = None
    gene = None

    # UniProt ID
    m = re.search(r"\b(?:sp|tr)\|([A-Z0-9]+(?:-[0-9]+)?)\|", h)
    if m:
        ident = m.group(1).upper()

    # Ensembl peptide ID
    if ident is None:
        m = re.match(r"^>(ENSP[0-9A-Za-z]+)", h)
        if m:
            ident = m.group(1).upper()

    # Gene symbol patterns
    m = re.search(r"\bGN=([A-Za-z0-9\-]+)\b", h)
    if m:
        gene = m.group(1)
    else:
        m = re.search(r"\bgene_symbol:([A-Za-z0-9\-]+)\b", h)
        if m:
            gene = m.group(1)
        else:
            m = re.search(r"\bGene_Symbol=([A-Za-z0-9\-]+)\b", h)
            if m:
                gene = m.group(1)

    if ident is None or gene is None:
        return (None, None)

    return (_strip_isoform(ident), _normalise_gene(gene))


def build_id_to_gene_map(fasta_path: Path) -> Dict[str, str]:
    """
    Build identifier -> gene_key mapping from a FASTA file.

    Parameters
    ----------
    fasta_path : Path
        Local FASTA path.

    Returns
    -------
    Dict[str, str]
        Mapping of identifier (base) to gene_key.
    """
    mapping: Dict[str, str] = {}
    with fasta_path.open("r", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            ident, gene = _parse_fasta_header(header=line)
            if ident is None or gene is None:
                continue
            mapping[ident] = gene
    return mapping


def _extract_ids(value: str) -> List[str]:
    """
    Extract identifiers from MaxQuant protein ID fields.

    Parameters
    ----------
    value : str
        Field containing IDs (semicolon-separated).

    Returns
    -------
    List[str]
        List of base identifiers.
    """
    if value is None:
        return []
    ids = []
    for token in str(value).split(";"):
        t = token.strip()
        if t == "":
            continue
        if t.startswith("CON__"):
            continue
        ids.append(_strip_isoform(t))
    return ids


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    if args.parameters_txt is not None and args.parameters_txt.exists():
        # purely informational (do not fail if path is foreign)
        for line in args.parameters_txt.read_text(errors="replace").splitlines():
            if "fasta" in line.lower() and "\t" in line:
                print("MaxQuant parameters FASTA entry:", line.strip())
                break

    if not args.fasta.exists():
        raise FileNotFoundError(f"FASTA not found: {args.fasta}")

    id_to_gene = build_id_to_gene_map(fasta_path=args.fasta)
    print(f"Built mapping entries: {len(id_to_gene)}")

    df = _safe_read_tsv(path=args.protein_groups_tsv)

    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df[df[col].fillna("") != "+"].copy()

    if "Majority protein IDs" not in df.columns:
        raise ValueError("Expected column 'Majority protein IDs' not found.")

    # Map each protein group to a gene_key using the first ID that maps
    gene_keys: List[Optional[str]] = []
    for v in df["Majority protein IDs"].fillna("").astype(str).tolist():
        g = None
        for ident in _extract_ids(value=v):
            if ident in id_to_gene:
                g = id_to_gene[ident]
                break
        gene_keys.append(g)

    df["gene_key"] = gene_keys
    df = df[df["gene_key"].notna()].copy()

    peptide_cols = [c for c in df.columns if str(c).startswith("Peptides ")]
    if not peptide_cols:
        raise ValueError("No per-sample peptide count columns found with prefix 'Peptides '.")

    peptides = df[peptide_cols].apply(pd.to_numeric, errors="coerce")
    present = peptides.notna() & (peptides > float(args.presence_threshold))

    df["_n_samples_present"] = present.sum(axis=1).astype(int)
    df["_n_samples_total"] = int(len(peptide_cols))
    df["_present_any"] = df["_n_samples_present"] > 0
    df["_present_fraction"] = df["_n_samples_present"] / df["_n_samples_total"]

    df["_pep_unique"] = pd.to_numeric(df.get("Peptide counts (unique)", np.nan), errors="coerce")
    df["_pep_razor_unique"] = pd.to_numeric(df.get("Peptide counts (razor+unique)", np.nan), errors="coerce")
    df["_coverage_pct"] = pd.to_numeric(df.get("Sequence coverage [%]", np.nan), errors="coerce")
    df["_q_value"] = pd.to_numeric(df.get("Q-value", np.nan), errors="coerce")

    gene = (
        df.groupby("gene_key", as_index=False)
        .agg(
            prot_present_any=("_present_any", "max"),
            prot_present_fraction=("_present_fraction", "max"),
            prot_n_samples_present=("_n_samples_present", "max"),
            prot_n_samples_total=("_n_samples_total", "max"),
            prot_peptide_counts_unique_max=("_pep_unique", "max"),
            prot_peptide_counts_razor_unique_max=("_pep_razor_unique", "max"),
            prot_sequence_coverage_pct_max=("_coverage_pct", "max"),
            prot_min_q_value=("_q_value", "min"),
        )
        .sort_values(
            by=["prot_present_any", "prot_present_fraction", "prot_peptide_counts_unique_max"],
            ascending=[False, False, False],
            na_position="last",
        )
        .reset_index(drop=True)
    )

    if args.write_per_sample_columns:
        present_bool = present.copy()
        present_bool.columns = [f"prot_present__{c.replace('Peptides ', '').strip()}" for c in present_bool.columns]
        tmp = pd.concat([df[["gene_key"]].reset_index(drop=True), present_bool.reset_index(drop=True)], axis=1)
        per_sample = tmp.groupby("gene_key", as_index=False).max()
        gene = gene.merge(per_sample, on="gene_key", how="left")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    gene.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Genes in output: {len(gene)}")
    print(f"Genes present in >=1 sample: {int(gene['prot_present_any'].sum())}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
