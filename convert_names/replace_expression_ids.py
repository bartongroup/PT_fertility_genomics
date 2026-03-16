#!/usr/bin/env python3

from __future__ import annotations

import argparse
import pandas as pd


def strip_version(value: str) -> str:
    return str(value).split(".")[0].strip()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replace Ensembl IDs in the first column of a TSV using a mapping table."
    )
    parser.add_argument("--expression_tsv", required=True)
    parser.add_argument("--mapping_tsv", required=True)
    parser.add_argument("--output_tsv", required=True)
    args = parser.parse_args()

    expr_df = pd.read_csv(args.expression_tsv, sep="\t", dtype=str)
    map_df = pd.read_csv(args.mapping_tsv, sep="\t", dtype=str).fillna("")

    required = {"original_gene_id", "mapped_gene_name"}
    missing = required - set(map_df.columns)
    if missing:
        raise ValueError(
            f"Mapping file missing required columns: {', '.join(sorted(missing))}"
        )

    mapping = {
        strip_version(row["original_gene_id"]): row["mapped_gene_name"].strip()
        for _, row in map_df.iterrows()
    }

    first_col = expr_df.columns[0]
    old_ids = expr_df[first_col].astype(str).map(strip_version)

    new_names = []
    for gene_id in old_ids:
        mapped = mapping.get(gene_id, gene_id)
        new_names.append(mapped if mapped else gene_id)

    expr_df[first_col] = new_names
    expr_df.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
