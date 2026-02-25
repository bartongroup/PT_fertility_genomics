#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_pxd037531_maxquant_proteingroups.py

Summarise MaxQuant proteinGroups.txt into a gene-level proteomics evidence TSV,
using a local Ensembl/UniProt FASTA to map protein identifiers to gene symbols.

This is designed for PXD037531 (Greither et al., 2023), where MaxQuant protein
IDs are Ensembl peptide identifiers with version suffixes (e.g. ENSP... .5).

Key features
------------
- Robust identifier normalisation to ensure MaxQuant IDs match FASTA IDs:
  - strips Ensembl version suffix (e.g. .5)
  - strips isoform suffix (e.g. -1)
  - strips contaminants prefix (CON__)
  - supports UniProt pipe headers (sp|P12345|...)
- Filters MaxQuant rows flagged as Reverse / Potential contaminant /
  Only identified by site.
- Calls per-sample presence using per-sample peptide count columns:
  'Peptides <sample>' > presence_threshold (default 0).
- Produces a gene-level TSV suitable for merging into downstream summaries.

Outputs
-------
TSV with columns:
- gene_key
- prot_present_any
- prot_present_fraction
- prot_n_samples_present
- prot_n_samples_total
- prot_peptide_counts_razor_unique_max
- prot_peptide_counts_unique_max
- prot_sequence_coverage_pct_max
- prot_min_q_value
- plus optional per-sample presence columns (prot_present__<sample>)

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
        description="Summarise MaxQuant proteinGroups.txt into a gene-level TSV."
    )
    parser.add_argument(
        "--protein_groups_tsv",
        required=True,
        type=Path,
        help="Path to MaxQuant proteinGroups.txt (tab-delimited).",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Local FASTA used for identifier->gene mapping (Ensembl pep or UniProt FASTA).",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        type=Path,
        help="Output gene-level TSV path.",
    )
    parser.add_argument(
        "--parameters_txt",
        required=False,
        type=Path,
        default=None,
        help="Optional MaxQuant parameters.txt for logging provenance only.",
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


def _normalise_gene_symbol(value: str) -> str:
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
    v = str(value).strip()
    v = re.sub(r"\s+", "", v)
    return v.upper()


def _normalise_identifier(raw_id: str) -> str:
    """
    Normalise MaxQuant/FASTA identifiers so they match between sources.

    Handles:
    - ENSP00000... .1 / .5 versions (strip after '.')
    - isoforms like -1 / -2 (strip after '-')
    - contaminants 'CON__' prefix (strip it)
    - UniProt-style pipes (sp|P12345|...) (keep accession)

    Parameters
    ----------
    raw_id : str
        Identifier token from MaxQuant or FASTA header.

    Returns
    -------
    str
        Normalised identifier.
    """
    s = str(raw_id).strip()

    if s.startswith("CON__"):
        s = s.replace("CON__", "", 1)

    # UniProt pipe format
    if "|" in s and (s.startswith("sp|") or s.startswith("tr|")):
        parts = s.split("|")
        if len(parts) >= 2:
            s = parts[1]

    # Generic pipe format
    if "|" in s:
        s = s.split("|")[0]

    # Strip isoform and version suffixes
    s = s.split("-")[0]
    s = s.split(".")[0]

    return s.upper()


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
    Parse identifier and gene symbol from a FASTA header.

    Supported examples
    ------------------
    - UniProt:
      >sp|P12345|ENTRY_HUMAN ... GN=GENE ...
    - Ensembl pep:
      >ENSP00000444533 pep:known chromosome:GRCh38:... gene:ENSG... transcript:ENST... gene_symbol:ABC1 ...
    - Custom:
      >A2A4G1 ... Gene_Symbol=Krt15 ...

    Parameters
    ----------
    header : str
        FASTA header line.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (normalised_identifier, normalised_gene_symbol) or (None, None).
    """
    h = header.strip()
    if not h.startswith(">"):
        return (None, None)

    ident: Optional[str] = None
    gene: Optional[str] = None

    # Prefer UniProt ID if present
    m = re.search(r"\b(?:sp|tr)\|([A-Z0-9]+(?:-[0-9]+)?)\|", h)
    if m:
        ident = m.group(1)

    # Otherwise, take the first token after '>'
    if ident is None:
        m = re.match(r"^>(\S+)", h)
        if m:
            ident = m.group(1)

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

    return (_normalise_identifier(ident), _normalise_gene_symbol(gene))


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
        Mapping of normalised identifier to gene_key.
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
    Extract normalised identifiers from MaxQuant 'Majority protein IDs'.

    Parameters
    ----------
    value : str
        Field containing semicolon-separated IDs.

    Returns
    -------
    List[str]
        List of normalised identifiers.
    """
    if value is None:
        return []
    ids: List[str] = []
    for token in str(value).split(";"):
        t = token.strip()
        if t == "":
            continue
        ids.append(_normalise_identifier(t))
    return ids


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    if args.parameters_txt is not None and args.parameters_txt.exists():
        txt = args.parameters_txt.read_text(errors="replace").splitlines()
        for line in txt:
            if "fasta" in line.lower():
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

    # Map each protein group to a gene key using the first identifier that maps
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
    df["_pep_razor_unique"] = pd.to_numeric(
        df.get("Peptide counts (razor+unique)", np.nan), errors="coerce"
    )
    df["_coverage_pct"] = pd.to_numeric(df.get("Sequence coverage [%]", np.nan), errors="coerce")
    df["_q_value"] = pd.to_numeric(df.get("Q-value", np.nan), errors="coerce")

    gene = (
        df.groupby("gene_key", as_index=False)
        .agg(
            prot_present_any=("_present_any", "max"),
            prot_present_fraction=("_present_fraction", "max"),
            prot_n_samples_present=("_n_samples_present", "max"),
            prot_n_samples_total=("_n_samples_total", "max"),
            prot_peptide_counts_razor_unique_max=("_pep_razor_unique", "max"),
            prot_peptide_counts_unique_max=("_pep_unique", "max"),
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
        present_bool.columns = [
            f"prot_present__{c.replace('Peptides ', '').strip()}" for c in present_bool.columns
        ]
        tmp = pd.concat(
            [df[["gene_key"]].reset_index(drop=True), present_bool.reset_index(drop=True)],
            axis=1,
        )
        per_sample = tmp.groupby("gene_key", as_index=False).max()
        gene = gene.merge(per_sample, on="gene_key", how="left")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    gene.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Genes in output: {len(gene)}")
    print(f"Genes present in >=1 sample: {int(gene['prot_present_any'].sum())}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()