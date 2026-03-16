#!/usr/bin/env python3
"""
replace_expression_ids_with_gene_names.py

Replace the first column of an expression matrix using a mapping table with:
    original_gene_id
    mapped_gene_name

Writes a new TSV with the first column replaced.

Usage
-----
python replace_expression_ids_with_gene_names.py \
    --expression_tsv gtex_v11_tissue_medians_trinity.tsv \
    --mapping_tsv gtex_v11_tissue_medians_gene_names.tsv.mapping_table.tsv \
    --output_tsv gtex_v11_tissue_medians_gene_symbols.tsv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd


def setup_logging(verbose: bool) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Replace first column of expression matrix using mapping table."
    )
    parser.add_argument(
        "--expression_tsv",
        required=True,
        help="Original expression matrix TSV.",
    )
    parser.add_argument(
        "--mapping_tsv",
        required=True,
        help="Mapping table TSV with original_gene_id and mapped_gene_name.",
    )
    parser.add_argument(
        "--output_tsv",
        required=True,
        help="Output TSV with replaced first column.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def strip_ensembl_version(gene_id: str) -> str:
    """Remove version suffix from Ensembl gene IDs."""
    return str(gene_id).split(".")[0].strip()


def main() -> None:
    """Run the replacement."""
    args = parse_args()
    setup_logging(args.verbose)

    logging.info("Loading expression matrix: %s", args.expression_tsv)
    expr_df = pd.read_csv(args.expression_tsv, sep="\t", dtype=str)
    first_col = expr_df.columns[0]

    logging.info("Loading mapping table: %s", args.mapping_tsv)
    map_df = pd.read_csv(args.mapping_tsv, sep="\t", dtype=str).fillna("")

    required_cols = {"original_gene_id", "mapped_gene_name"}
    missing = required_cols - set(map_df.columns)
    if missing:
        raise ValueError(
            f"Mapping file is missing required columns: {', '.join(sorted(missing))}"
        )

    # Build dictionary
    mapping = {
        strip_ensembl_version(row["original_gene_id"]): row["mapped_gene_name"].strip()
        for _, row in map_df.iterrows()
        if row["original_gene_id"].strip() != ""
    }

    logging.info("Applying mapping to first column: %s", first_col)

    original_ids = expr_df[first_col].astype(str).map(strip_ensembl_version)

    replaced = []
    n_changed = 0

    for gene_id in original_ids:
        new_name = mapping.get(gene_id, gene_id)
        replaced.append(new_name)
        if new_name != gene_id:
            n_changed += 1

    expr_df[first_col] = replaced

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    expr_df.to_csv(output_path, sep="\t", index=False)

    logging.info("Wrote output matrix: %s", output_path)
    logging.info("Rows changed: %d", n_changed)
    logging.info("Total rows: %d", expr_df.shape[0])


if __name__ == "__main__":
    main()
