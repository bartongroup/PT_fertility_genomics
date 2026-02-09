#!/usr/bin/env python3
"""
Extract a unique gene list from HPO-derived male infertility outputs.

This script reads male_infertility_genes_summary.tsv (produced by the ontology pipeline)
and writes:
- a unique list of gene symbols
- a unique list of Entrez Gene IDs
- a combined table (Entrez ID + symbol)

Outputs are TSV (tab-separated).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Extract unique gene lists from HPO pipeline outputs.")
    parser.add_argument(
        "--genes_summary_tsv",
        required=True,
        help="Path to male_infertility_genes_summary.tsv",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for gene list TSV files",
    )
    return parser.parse_args()


def _clean_value(value: str) -> str:
    """
    Clean a string value for consistent parsing.

    Args:
        value: Raw value.

    Returns:
        Cleaned value.
    """
    return value.strip()


def load_gene_pairs(genes_summary_tsv: Path) -> Set[Tuple[str, str]]:
    """
    Load (Entrez gene_id, gene_symbol) pairs from the genes summary TSV.

    Args:
        genes_summary_tsv: Path to male_infertility_genes_summary.tsv.

    Returns:
        Set of (gene_id, gene_symbol) pairs.
    """
    pairs: Set[Tuple[str, str]] = set()

    with genes_summary_tsv.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = _clean_value(row.get("gene_id", ""))
            gene_symbol = _clean_value(row.get("gene_symbol", ""))

            if not gene_id and not gene_symbol:
                continue

            pairs.add((gene_id, gene_symbol))

    return pairs


def write_gene_lists(out_dir: Path, pairs: Set[Tuple[str, str]]) -> None:
    """
    Write gene symbol list, Entrez ID list, and combined table.

    Args:
        out_dir: Output directory.
        pairs: Set of (gene_id, gene_symbol) pairs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    symbols: Set[str] = {s for _, s in pairs if s}
    entrez_ids: Set[str] = {g for g, _ in pairs if g}

    combined_path = out_dir / "male_infertility_gene_list.tsv"
    symbols_path = out_dir / "male_infertility_gene_symbols.tsv"
    entrez_path = out_dir / "male_infertility_entrez_ids.tsv"

    with combined_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol"])
        for gene_id, gene_symbol in sorted(pairs, key=lambda x: (x[1], x[0])):
            writer.writerow([gene_id, gene_symbol])

    with symbols_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_symbol"])
        for symbol in sorted(symbols):
            writer.writerow([symbol])

    with entrez_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id"])
        for gene_id in sorted(entrez_ids, key=lambda x: int(x) if x.isdigit() else 10**12):
            writer.writerow([gene_id])


def main() -> None:
    """
    Run gene list extraction.
    """
    args = parse_args()
    genes_summary_tsv = Path(args.genes_summary_tsv)
    out_dir = Path(args.out_dir)

    if not genes_summary_tsv.exists():
        raise FileNotFoundError(f"Input file not found: {genes_summary_tsv}")

    pairs = load_gene_pairs(genes_summary_tsv=genes_summary_tsv)
    write_gene_lists(out_dir=out_dir, pairs=pairs)


if __name__ == "__main__":
    main()
