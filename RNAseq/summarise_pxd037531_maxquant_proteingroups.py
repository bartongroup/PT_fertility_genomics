#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_pxd037531_maxquant_proteingroups.py

Convert MaxQuant 'proteinGroups.txt' (PXD037531; Greither et al., 2023) into a
gene-level proteomics evidence table suitable for merging into a master sperm /
testis candidate gene summary table.

This script is designed for datasets where MaxQuant reports protein group-level
identifications across many samples (columns such as 'LFQ intensity <sample>' or
'Intensity <sample>' and peptide count columns).

Filtering
---------
By default removes entries flagged as:
- Reverse
- Potential contaminant
- Only identified by site

Gene symbol extraction
----------------------
Primary: parse 'GN=<symbol>' from 'Fasta headers' (MaxQuant/UniProt style).
Fallback: attempt to parse gene symbol from 'Majority protein IDs' where a
'GN=' token exists (rare).
If no gene symbol can be extracted, the row is discarded (can be relaxed with
--keep_missing_gene).

Presence definition
-------------------
Presence per sample is called if the chosen per-sample quantitative field is
present and > 0 for that sample.
- Default quant field: LFQ intensity
- Alternative: Intensity (use --use_intensity)

Outputs
-------
A tab-separated file with one row per gene_key, including:
- gene_key
- prot_present_any
- prot_present_fraction
- prot_n_samples_present
- prot_n_samples_total
- prot_peptide_counts_razor_unique_max
- prot_peptide_counts_unique_max
- prot_sequence_coverage_pct_max
- prot_min_q_value
- prot_max_lfq_intensity (or prot_max_intensity)
- optional per-sample presence columns (prot_present__<sample>)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import csv
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
        description=(
            "Summarise MaxQuant proteinGroups.txt into a gene-level proteomics "
            "evidence TSV for downstream annotation."
        )
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
        "--keep_missing_gene",
        action="store_true",
        help=(
            "If set, retain rows where a gene symbol cannot be extracted, "
            "using the first Majority protein ID as gene_key fallback."
        ),
    )
    parser.add_argument(
        "--use_intensity",
        action="store_true",
        help=(
            "If set, use per-sample 'Intensity ' columns for presence calls "
            "instead of 'LFQ intensity '."
        ),
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


def _extract_gn_from_fasta_headers(headers: str) -> Optional[str]:
    """
    Extract GN= gene symbol from a MaxQuant 'Fasta headers' field.

    Parameters
    ----------
    headers : str
        'Fasta headers' field (often semicolon-separated).

    Returns
    -------
    Optional[str]
        Normalised gene symbol if found, else None.
    """
    if headers is None or str(headers).strip() == "":
        return None
    m = re.search(r"\bGN=([A-Za-z0-9\-]+)\b", str(headers))
    if not m:
        return None
    return _normalise_gene_symbol(m.group(1))


def _fallback_gene_key_from_majority_ids(majority_ids: str) -> Optional[str]:
    """
    Fallback gene key if GN= is missing: use first Majority protein ID.

    Parameters
    ----------
    majority_ids : str
        'Majority protein IDs' field (often ';' separated UniProt accessions).

    Returns
    -------
    Optional[str]
        A fallback key (uppercase) or None.
    """
    if majority_ids is None or str(majority_ids).strip() == "":
        return None
    first = str(majority_ids).split(";")[0].strip()
    if first == "":
        return None
    return first.upper()


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


def _find_sample_columns(df: pd.DataFrame, prefix: str) -> List[str]:
    """
    Identify per-sample columns by a prefix (e.g. 'LFQ intensity ').

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    prefix : str
        Column prefix.

    Returns
    -------
    List[str]
        Matching columns.
    """
    return [c for c in df.columns if str(c).startswith(prefix)]


def build_gene_level_table(
    protein_groups_path: Path,
    use_intensity: bool,
    presence_threshold: float,
    keep_missing_gene: bool,
    write_per_sample_columns: bool,
) -> pd.DataFrame:
    """
    Build gene-level proteomics evidence table from MaxQuant proteinGroups.txt.

    Parameters
    ----------
    protein_groups_path : Path
        Path to proteinGroups.txt (tab-delimited).
    use_intensity : bool
        If True, use 'Intensity ' columns; else use 'LFQ intensity ' columns.
    presence_threshold : float
        Threshold above which a sample is considered present.
    keep_missing_gene : bool
        If True, keep rows lacking GN= by using a fallback gene_key.
    write_per_sample_columns : bool
        If True, include per-sample presence columns.

    Returns
    -------
    pd.DataFrame
        Gene-level table.
    """
    df = pd.read_csv(protein_groups_path, sep="\t", dtype=str, engine="python")
    df.columns = [str(c).strip() for c in df.columns]

    # Standard MaxQuant filter flags
    for flag_col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if flag_col in df.columns:
            df = df[df[flag_col].fillna("") != "+"].copy()

    # Gene symbol extraction
    if "Fasta headers" not in df.columns:
        raise ValueError("Expected column 'Fasta headers' not found in proteinGroups.txt")

    df["gene_key"] = df["Fasta headers"].map(_extract_gn_from_fasta_headers)

    if keep_missing_gene:
        if "Majority protein IDs" in df.columns:
            mask = df["gene_key"].isna()
            df.loc[mask, "gene_key"] = df.loc[mask, "Majority protein IDs"].map(
                _fallback_gene_key_from_majority_ids
            )

    df = df[df["gene_key"].notna()].copy()
    df["gene_key"] = df["gene_key"].astype(str)

    # Quant columns for presence calls
    quant_prefix = "Intensity " if use_intensity else "LFQ intensity "
    sample_cols = _find_sample_columns(df=df, prefix=quant_prefix)

    if not sample_cols:
        raise ValueError(
            f"No per-sample quant columns found with prefix '{quant_prefix}'. "
            f"Available columns include: {df.columns[:30].tolist()} ..."
        )

    # Parse per-sample quantities and call presence
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

    # Summarise to gene-level (max evidence across protein groups for that gene)
    agg = {
        "_present_any": "max",
        "_present_fraction": "max",
        "_n_samples_present": "max",
        "_n_samples_total": "max",
        "_pep_razor_unique": "max",
        "_pep_unique": "max",
        "_coverage_pct": "max",
        "_q_value": "min",
    }

    gene = df.groupby("gene_key", as_index=False).agg(agg)

    # Rename outputs
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
        }
    )

    # Add max overall intensity metric at gene level
    df["_max_quant"] = quant.max(axis=1)
    max_quant_gene = df.groupby("gene_key", as_index=False)["_max_quant"].max()
    max_quant_col = "prot_max_intensity" if use_intensity else "prot_max_lfq_intensity"
    max_quant_gene = max_quant_gene.rename(columns={"_max_quant": max_quant_col})
    gene = gene.merge(max_quant_gene, on="gene_key", how="left")

    if write_per_sample_columns:
        # Per-sample presence columns collapsed to gene-level (any protein group present)
        present_bool = present.copy()
        present_bool.columns = [f"prot_present__{c.replace(quant_prefix, '').strip()}" for c in present_bool.columns]
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

    out = build_gene_level_table(
        protein_groups_path=args.protein_groups_tsv,
        use_intensity=bool(args.use_intensity),
        presence_threshold=float(args.presence_threshold),
        keep_missing_gene=bool(args.keep_missing_gene),
        write_per_sample_columns=bool(args.write_per_sample_columns),
    )

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Proteins table: {args.protein_groups_tsv}")
    print(f"Genes in output: {len(out)}")
    print(f"Genes present in >=1 sample: {int(out['prot_present_any'].sum())}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
