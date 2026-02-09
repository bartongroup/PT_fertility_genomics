#!/usr/bin/env python3
"""
Filter ClinVar variant_summary.txt.gz to a user-provided gene symbol list.

This script is designed for cluster use:
- Streams the gzip input (no need to decompress to disk)
- Writes tab-separated outputs (TSV)
- Keeps only rows where GeneSymbol matches your gene list

Inputs:
- variant_summary.txt.gz from ClinVar FTP (tab_delimited)
- a TSV file containing a column named 'gene_symbol' (one symbol per row)

Outputs:
- a filtered TSV containing the original ClinVar columns (or selected columns if requested)
- a small summary TSV with counts per ClinicalSignificance and ReviewStatus
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Filter ClinVar variant_summary.txt.gz by a gene symbol list."
    )
    parser.add_argument(
        "--variant_summary_gz",
        required=True,
        help="Path to ClinVar variant_summary.txt.gz",
    )
    parser.add_argument(
        "--gene_symbols_tsv",
        required=True,
        help="TSV with a 'gene_symbol' column (e.g. male_infertility_gene_symbols.tsv)",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV path for filtered ClinVar rows",
    )
    parser.add_argument(
        "--out_summary_tsv",
        required=True,
        help="Output TSV path for summary counts",
    )
    parser.add_argument(
        "--columns",
        default="",
        help=(
            "Optional comma-separated list of columns to write (must match ClinVar header exactly). "
            "If omitted, all columns are written."
        ),
    )
    parser.add_argument(
        "--require_clinsig",
        default="",
        help=(
            "Optional filter for ClinicalSignificance (case-insensitive substring match). "
            "Examples: 'Pathogenic', 'Likely_pathogenic', 'Benign'. "
            "If omitted, keep all ClinicalSignificance values."
        ),
    )
    parser.add_argument(
        "--require_review_status",
        default="",
        help=(
            "Optional filter for ReviewStatus (case-insensitive substring match). "
            "Example: 'criteria provided, multiple submitters, no conflicts'. "
            "If omitted, keep all ReviewStatus values."
        ),
    )
    return parser.parse_args()


def load_gene_symbols(gene_symbols_tsv: Path) -> Set[str]:
    """
    Load a set of gene symbols from a TSV with a 'gene_symbol' column.

    Args:
        gene_symbols_tsv: Path to TSV.

    Returns:
        Set of gene symbols (uppercased).
    """
    symbols: Set[str] = set()
    with gene_symbols_tsv.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "gene_symbol" not in reader.fieldnames or reader.fieldnames is None:
            raise ValueError(
                "gene_symbols_tsv must contain a column named 'gene_symbol'."
            )
        for row in reader:
            sym = (row.get("gene_symbol") or "").strip()
            if sym:
                symbols.add(sym.upper())
    return symbols


def _normalise_filter_value(value: str) -> str:
    """
    Normalise a filter value for case-insensitive substring matching.

    Args:
        value: Raw string.

    Returns:
        Normalised string.
    """
    return value.strip().lower()


def _passes_optional_filter(value: str, required_substring: str) -> bool:
    """
    Decide whether a value passes an optional substring filter.

    Args:
        value: Field value from input.
        required_substring: Substring filter (lowercased). If empty, pass.

    Returns:
        True if passes filter.
    """
    if not required_substring:
        return True
    return required_substring in (value or "").lower()


def filter_variant_summary(
    variant_summary_gz: Path,
    gene_symbols: Set[str],
    out_tsv: Path,
    out_summary_tsv: Path,
    columns_to_write: Optional[List[str]],
    require_clinsig: str,
    require_review_status: str,
) -> None:
    """
    Stream-filter ClinVar variant_summary.txt.gz by gene symbols and optional interpretation filters.

    Args:
        variant_summary_gz: Path to gzipped ClinVar variant_summary.
        gene_symbols: Set of gene symbols (uppercased).
        out_tsv: Output TSV of filtered variants.
        out_summary_tsv: Output TSV summary counts.
        columns_to_write: Optional list of column names to write. If None, write all.
        require_clinsig: Optional substring for ClinicalSignificance filtering (lowercased).
        require_review_status: Optional substring for ReviewStatus filtering (lowercased).
    """
    counts_clinsig: Counter[str] = Counter()
    counts_review: Counter[str] = Counter()
    kept_rows = 0
    total_rows = 0

    with gzip.open(variant_summary_gz, mode="rt", encoding="utf-8", newline="") as handle_in:
        reader = csv.DictReader(handle_in, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("ClinVar variant_summary appears to have no header.")

        fieldnames = reader.fieldnames

        # Gene symbol column name in ClinVar tab-delimited is typically 'GeneSymbol'
        if "GeneSymbol" not in fieldnames:
            raise ValueError(
                "Could not find 'GeneSymbol' column in ClinVar variant_summary header."
            )

        # Optional interpretation columns (present in standard variant_summary)
        clinsig_col = "ClinicalSignificance" if "ClinicalSignificance" in fieldnames else None
        review_col = "ReviewStatus" if "ReviewStatus" in fieldnames else None

        if columns_to_write:
            missing = [c for c in columns_to_write if c not in fieldnames]
            if missing:
                raise ValueError(
                    f"Requested columns not present in variant_summary header: {missing}"
                )
            out_fields = columns_to_write
        else:
            out_fields = fieldnames

        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_summary_tsv.parent.mkdir(parents=True, exist_ok=True)

        with out_tsv.open(mode="w", encoding="utf-8", newline="") as handle_out:
            writer = csv.DictWriter(handle_out, delimiter="\t", fieldnames=out_fields, extrasaction="ignore")
            writer.writeheader()

            for row in reader:
                total_rows += 1

                gene_symbol = (row.get("GeneSymbol") or "").strip().upper()
                if gene_symbol not in gene_symbols:
                    continue

                clinsig_val = row.get(clinsig_col, "") if clinsig_col else ""
                review_val = row.get(review_col, "") if review_col else ""

                if not _passes_optional_filter(clinsig_val, require_clinsig):
                    continue
                if not _passes_optional_filter(review_val, require_review_status):
                    continue

                writer.writerow(row)
                kept_rows += 1

                if clinsig_col:
                    counts_clinsig[clinsig_val or ""] += 1
                if review_col:
                    counts_review[review_val or ""] += 1

    with out_summary_tsv.open(mode="w", encoding="utf-8", newline="") as handle_sum:
        writer = csv.writer(handle_sum, delimiter="\t")
        writer.writerow(["metric", "category", "count"])
        writer.writerow(["rows", "total", total_rows])
        writer.writerow(["rows", "kept", kept_rows])

        for k, v in counts_clinsig.most_common():
            writer.writerow(["ClinicalSignificance", k, v])

        for k, v in counts_review.most_common():
            writer.writerow(["ReviewStatus", k, v])


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    variant_summary_gz = Path(args.variant_summary_gz)
    gene_symbols_tsv = Path(args.gene_symbols_tsv)
    out_tsv = Path(args.out_tsv)
    out_summary_tsv = Path(args.out_summary_tsv)

    if not variant_summary_gz.exists():
        raise FileNotFoundError(f"Missing input: {variant_summary_gz}")
    if not gene_symbols_tsv.exists():
        raise FileNotFoundError(f"Missing input: {gene_symbols_tsv}")

    gene_symbols = load_gene_symbols(gene_symbols_tsv=gene_symbols_tsv)

    columns_to_write = None
    if args.columns.strip():
        columns_to_write = [c.strip() for c in args.columns.split(",") if c.strip()]

    require_clinsig = _normalise_filter_value(args.require_clinsig)
    require_review_status = _normalise_filter_value(args.require_review_status)

    filter_variant_summary(
        variant_summary_gz=variant_summary_gz,
        gene_symbols=gene_symbols,
        out_tsv=out_tsv,
        out_summary_tsv=out_summary_tsv,
        columns_to_write=columns_to_write,
        require_clinsig=require_clinsig,
        require_review_status=require_review_status,
    )


if __name__ == "__main__":
    main()
