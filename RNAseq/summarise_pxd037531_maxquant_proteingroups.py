#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_pxd037531_maxquant_proteingroups.py

Summarise MaxQuant proteinGroups.txt into a gene-level proteomics evidence TSV.

This version maps UniProt accessions to gene symbols via either:
1) A FASTA file used by MaxQuant (preferred; fully offline), or
2) A prebuilt UniProt accession → gene_key TSV mapping file.

Filtering
---------
Removes rows flagged as Reverse, Potential contaminant, Only identified by site.

Outputs
-------
Tab-separated gene-level evidence table for merging into master summary tables.
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

import sys

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise MaxQuant proteinGroups.txt into a gene-level proteomics TSV."
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
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional FASTA used by MaxQuant. If provided, build UniProt accession → gene mapping "
            "from FASTA headers."
        ),
    )
    parser.add_argument(
        "--uniprot_to_gene_tsv",
        required=False,
        type=Path,
        default=None,
        help=(
            "Optional TSV mapping file with columns: uniprot_accession, gene_key. "
            "Used if --fasta is not provided."
        ),
    )
    parser.add_argument(
        "--use_intensity",
        action="store_true",
        help="If set, use per-sample 'Intensity ' columns instead of 'LFQ intensity '.",
    )
    parser.add_argument(
        "--presence_threshold",
        required=False,
        default=0.0,
        type=float,
        help="Threshold for calling presence from per-sample quant fields (default: 0).",
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


def _read_tsv(path: Path) -> pd.DataFrame:
    """
    Read a TSV safely, allowing very large fields.

    Parameters
    ----------
    path : Path
        Input file path.

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


def _extract_uniprot_accessions(value: str) -> List[str]:
    """
    Extract UniProt accessions from a MaxQuant protein ID field.

    Parameters
    ----------
    value : str
        Protein IDs field (often ';' separated accessions).

    Returns
    -------
    List[str]
        List of accessions (uppercased), isoform suffix preserved.
    """
    if value is None:
        return []
    items = [x.strip() for x in str(value).split(";") if x.strip() != ""]
    return [x.upper() for x in items]


def _strip_isoform(accession: str) -> str:
    """
    Strip UniProt isoform suffix (e.g. P12345-2 -> P12345).

    Parameters
    ----------
    accession : str
        UniProt accession.

    Returns
    -------
    str
        Base accession.
    """
    return accession.split("-")[0].upper()


def _parse_fasta_accession_and_gene(header: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse accession and gene symbol from a FASTA header.

    Supports common UniProt-style headers:
    >sp|P12345|ENTRY_HUMAN ... GN=GENE ...

    Parameters
    ----------
    header : str
        FASTA header line (without newline).

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (accession, gene_symbol) if detected, otherwise (None, None).
    """
    h = header.strip()
    if not h.startswith(">"):
        return (None, None)

    # Accession from UniProt-style header
    m_acc = re.search(r"\b(?:sp|tr)\|([A-Z0-9]+(?:-[0-9]+)?)\|", h)
    acc = m_acc.group(1).upper() if m_acc else None

    # Gene symbol token (if present)
    m_gn = re.search(r"\bGN=([A-Za-z0-9\-]+)\b", h)
    gene = _normalise_gene_symbol(m_gn.group(1)) if m_gn else None

    return (acc, gene)


def build_uniprot_to_gene_map_from_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Build UniProt accession → gene_key mapping from a FASTA file.

    Parameters
    ----------
    fasta_path : Path
        FASTA path.

    Returns
    -------
    Dict[str, str]
        Mapping of base accession (no isoform suffix) to gene_key.
    """
    mapping: Dict[str, str] = {}
    with fasta_path.open("r", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            acc, gene = _parse_fasta_accession_and_gene(header=line)
            if acc is None or gene is None:
                continue
            mapping[_strip_isoform(acc)] = gene
    return mapping


def load_uniprot_to_gene_map_from_tsv(mapping_tsv: Path) -> Dict[str, str]:
    """
    Load UniProt accession → gene_key mapping from TSV.

    Required columns: uniprot_accession, gene_key

    Parameters
    ----------
    mapping_tsv : Path
        TSV path.

    Returns
    -------
    Dict[str, str]
        Mapping of base accession to gene_key.
    """
    df = pd.read_csv(mapping_tsv, sep="\t", dtype=str)
    cols = [c.strip() for c in df.columns]
    df.columns = cols

    if "uniprot_accession" not in df.columns or "gene_key" not in df.columns:
        raise ValueError(
            "Mapping TSV must contain columns: uniprot_accession, gene_key"
        )

    out: Dict[str, str] = {}
    for _, row in df.iterrows():
        acc = str(row["uniprot_accession"]).strip().upper()
        gene = str(row["gene_key"]).strip()
        if acc == "" or gene == "":
            continue
        out[_strip_isoform(acc)] = _normalise_gene_symbol(gene)
    return out


def _to_numeric(series: pd.Series) -> pd.Series:
    """
    Convert a Series to numeric floats safely.

    Parameters
    ----------
    series : pd.Series
        Input series.

    Returns
    -------
    pd.Series
        Float series (NaN for non-numeric).
    """
    return pd.to_numeric(series, errors="coerce")


def build_gene_level_table(
    protein_groups_path: Path,
    uniprot_to_gene: Dict[str, str],
    use_intensity: bool,
    presence_threshold: float,
    write_per_sample_columns: bool,
) -> pd.DataFrame:
    """
    Build gene-level proteomics evidence table from MaxQuant proteinGroups.txt.

    Parameters
    ----------
    protein_groups_path : Path
        Path to proteinGroups.txt.
    uniprot_to_gene : Dict[str, str]
        Mapping from UniProt base accession to gene_key.
    use_intensity : bool
        If True, use 'Intensity ' columns; else use 'LFQ intensity ' columns.
    presence_threshold : float
        Threshold above which a sample is considered present.
    write_per_sample_columns : bool
        If True, include per-sample presence columns.

    Returns
    -------
    pd.DataFrame
        Gene-level table.
    """
    df = _read_tsv(path=protein_groups_path)

    # Standard MaxQuant filter flags
    for flag_col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if flag_col in df.columns:
            df = df[df[flag_col].fillna("") != "+"].copy()

    if "Majority protein IDs" not in df.columns:
        raise ValueError("Expected column 'Majority protein IDs' not found.")

    # Map protein groups → gene_key (use first mapped accession)
    maj = df["Majority protein IDs"].fillna("").astype(str)
    acc_lists = maj.map(_extract_uniprot_accessions)

    gene_keys: List[Optional[str]] = []
    for accs in acc_lists.tolist():
        g = None
        for a in accs:
            base = _strip_isoform(a)
            if base in uniprot_to_gene:
                g = uniprot_to_gene[base]
                break
        gene_keys.append(g)

    df["gene_key"] = gene_keys
    df = df[df["gene_key"].notna()].copy()
    df["gene_key"] = df["gene_key"].astype(str)

    # Quant columns for presence calls
    quant_prefix = "Intensity " if use_intensity else "LFQ intensity "
    sample_cols = [c for c in df.columns if str(c).startswith(quant_prefix)]
    if not sample_cols:
        raise ValueError(f"No per-sample quant columns found with prefix '{quant_prefix}'.")

    quant = df[sample_cols].apply(_to_numeric)
    present = quant.notna() & (quant > float(presence_threshold))

    df["_n_samples_present"] = present.sum(axis=1).astype(int)
    df["_n_samples_total"] = int(len(sample_cols))
    df["_present_any"] = df["_n_samples_present"] >= 1
    df["_present_fraction"] = df["_n_samples_present"] / df["_n_samples_total"]

    # Core evidence fields
    pep_razor_unique_col = "Peptide counts (razor+unique)"
    pep_unique_col = "Peptide counts (unique)"
    coverage_col = "Sequence coverage [%]"
    qval_col = "Q-value"

    df["_pep_razor_unique"] = _to_numeric(df[pep_razor_unique_col]) if pep_razor_unique_col in df.columns else np.nan
    df["_pep_unique"] = _to_numeric(df[pep_unique_col]) if pep_unique_col in df.columns else np.nan
    df["_coverage_pct"] = _to_numeric(df[coverage_col]) if coverage_col in df.columns else np.nan
    df["_q_value"] = _to_numeric(df[qval_col]) if qval_col in df.columns else np.nan

    df["_max_quant"] = quant.max(axis=1)

    agg = {
        "_present_any": "max",
        "_present_fraction": "max",
        "_n_samples_present": "max",
        "_n_samples_total": "max",
        "_pep_razor_unique": "max",
        "_pep_unique": "max",
        "_coverage_pct": "max",
        "_q_value": "min",
        "_max_quant": "max",
    }

    gene = df.groupby("gene_key", as_index=False).agg(agg)

    max_quant_col = "prot_max_intensity" if use_intensity else "prot_max_lfq_intensity"

    gene = gene.rename(
        columns={
            "_present_any": "prot_present_any",
            "_present_fraction": "prot_present_fraction",
            "_n_samples_present": "prot_n_samples_present",
            "_n_samples_total": "prot_n_samples_total",
            "_pep_razor_unique": "prot_peptide_counts_razor_unique_max",
            "_pep_unique": "prot_peptide_counts_unique_max",
            "_coverage_pct": "prot_sequence_coverage_pct_max",
            "_q_value": "prot_min_q_value",
            "_max_quant": max_quant_col,
        }
    )

    if write_per_sample_columns:
        present_bool = present.copy()
        present_bool.columns = [
            f"prot_present__{c.replace(quant_prefix, '').strip()}" for c in present_bool.columns
        ]
        tmp = pd.concat([df[["gene_key"]].reset_index(drop=True), present_bool.reset_index(drop=True)], axis=1)
        per_sample = tmp.groupby("gene_key", as_index=False).max()
        gene = gene.merge(per_sample, on="gene_key", how="left")

    gene = gene.sort_values(
        by=["prot_present_any", "prot_present_fraction", "prot_peptide_counts_razor_unique_max"],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)

    return gene


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    if args.fasta is not None:
        uniprot_to_gene = build_uniprot_to_gene_map_from_fasta(fasta_path=args.fasta)
    elif args.uniprot_to_gene_tsv is not None:
        uniprot_to_gene = load_uniprot_to_gene_map_from_tsv(mapping_tsv=args.uniprot_to_gene_tsv)
    else:
        raise ValueError("Provide either --fasta or --uniprot_to_gene_tsv for accession→gene mapping.")

    out = build_gene_level_table(
        protein_groups_path=args.protein_groups_tsv,
        uniprot_to_gene=uniprot_to_gene,
        use_intensity=bool(args.use_intensity),
        presence_threshold=float(args.presence_threshold),
        write_per_sample_columns=bool(args.write_per_sample_columns),
    )

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Genes in output: {len(out)}")
    print(f"Genes present in >=1 sample: {int(out['prot_present_any'].sum())}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()