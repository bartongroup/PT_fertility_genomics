#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_ensembl_ids_to_gene_symbols.py

Convert Ensembl gene IDs (ENSG...) to HGNC-style gene symbols (gene_key)
using a GENCODE GRCh38 GTF.

This script is designed to support workflows where some gene lists are provided
as Ensembl gene IDs (optionally including version suffixes like '.7') but
downstream scripts expect HGNC gene symbols.

Inputs
------
- A TSV containing Ensembl gene IDs. By default the column is 'gene_id', but it
  can be changed via --gene_id_column.
- A GENCODE GRCh38 GTF file (chr-based recommended).

Outputs
-------
- out_gene_symbols_tsv: TSV with at least column 'gene_key'
- out_mapping_report_tsv: TSV describing mapping status per input entry

Example
-------
python convert_ensembl_ids_to_gene_symbols.py \
  --ensembl_tsv ensembl_ids.tsv \
  --gene_id_column gene_id \
  --gencode_gtf gencode.vXX.annotation.chr.gtf.gz \
  --out_gene_symbols_tsv sperm_genes_as_symbols.tsv \
  --out_mapping_report_tsv ensembl_to_symbol_mapping_report.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
from typing import Dict, Optional, Sequence, Tuple

import pandas as pd


def setup_logger(*, verbose: bool) -> None:
    """Configure logging for the script."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _open_maybe_gzip(*, path: str):
    """Open a file that may be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def parse_gtf_gene_id_to_name(*, gencode_gtf: str) -> Dict[str, str]:
    """
    Parse a GENCODE GTF and return mapping of gene_id (without version) to gene_name.

    Notes
    -----
    - GENCODE gene_id values usually include a version suffix (e.g. ENSG... .7).
      We store both the full ID and the stripped ID for robust matching.
    """
    mapping: Dict[str, str] = {}
    mapping_stripped: Dict[str, str] = {}

    with _open_maybe_gzip(path=gencode_gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "gene":
                continue

            attrs = parts[8]
            attr_map = parse_gtf_attributes(attrs=attrs)
            gene_id = str(attr_map.get("gene_id", "")).strip()
            gene_name = str(attr_map.get("gene_name", "")).strip()
            if gene_id == "" or gene_name == "":
                continue

            mapping[gene_id] = gene_name
            mapping_stripped[strip_ensembl_version(ensembl_id=gene_id)] = gene_name

    # Prefer stripped mapping for general use
    mapping.update(mapping_stripped)
    if not mapping:
        raise ValueError("No gene_id -> gene_name mappings found in the GTF.")
    return mapping


def parse_gtf_attributes(*, attrs: str) -> Dict[str, str]:
    """Parse GTF attribute column into a dictionary."""
    out: Dict[str, str] = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if item == "":
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def strip_ensembl_version(*, ensembl_id: str) -> str:
    """Strip version suffix from an Ensembl ID (e.g. ENSG... .7 -> ENSG...)."""
    s = str(ensembl_id).strip()
    if s == "":
        return ""
    return s.split(".", 1)[0]


def read_ensembl_list(*, ensembl_tsv: str, gene_id_column: str) -> pd.DataFrame:
    """Read TSV containing Ensembl gene IDs."""
    df = pd.read_csv(ensembl_tsv, sep="\t", dtype=str).fillna("")
    if gene_id_column not in df.columns:
        raise ValueError(f"Input TSV must contain column '{gene_id_column}'.")
    df[gene_id_column] = df[gene_id_column].astype(str).str.strip()
    df = df[df[gene_id_column] != ""].drop_duplicates(subset=[gene_id_column]).copy()
    df["gene_id_stripped"] = df[gene_id_column].map(lambda x: strip_ensembl_version(ensembl_id=x))
    return df


def convert_ids(
    *, df: pd.DataFrame, id_to_name: Dict[str, str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convert Ensembl gene IDs to gene symbols.

    Returns
    -------
    out_df
        Contains 'gene_key' for mapped entries.
    report_df
        One row per input ID with mapping status.
    """
    report_rows = []
    out_rows = []

    for row in df.itertuples(index=False):
        raw_id = str(getattr(row, df.columns[0]))
        stripped = str(getattr(row, "gene_id_stripped"))

        gene_name = id_to_name.get(raw_id)
        if gene_name is None:
            gene_name = id_to_name.get(stripped)

        if gene_name is None or str(gene_name).strip() == "":
            report_rows.append(
                {
                    "input_gene_id": raw_id,
                    "stripped_gene_id": stripped,
                    "status": "unmapped",
                    "gene_key": "",
                }
            )
            continue

        gene_key = str(gene_name).strip()
        report_rows.append(
            {
                "input_gene_id": raw_id,
                "stripped_gene_id": stripped,
                "status": "mapped",
                "gene_key": gene_key,
            }
        )
        out_rows.append({"gene_key": gene_key})

    out_df = pd.DataFrame(out_rows).drop_duplicates(subset=["gene_key"])
    report_df = pd.DataFrame(report_rows)
    return out_df, report_df


def write_tsv(*, df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame as TSV, creating parent directories if needed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert Ensembl gene IDs to HGNC-style gene symbols using a GENCODE GTF."
    )
    parser.add_argument("--ensembl_tsv", required=True, help="Input TSV containing Ensembl gene IDs.")
    parser.add_argument(
        "--gene_id_column",
        default="gene_id",
        help="Column name containing Ensembl IDs (default: gene_id).",
    )
    parser.add_argument("--gencode_gtf", required=True, help="GENCODE GRCh38 GTF (chr recommended).")
    parser.add_argument("--out_gene_symbols_tsv", required=True, help="Output TSV with column 'gene_key'.")
    parser.add_argument("--out_mapping_report_tsv", required=True, help="Output TSV mapping report.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv=argv)
    setup_logger(verbose=bool(args.verbose))

    logging.info("Reading Ensembl IDs: %s", args.ensembl_tsv)
    df = read_ensembl_list(ensembl_tsv=args.ensembl_tsv, gene_id_column=str(args.gene_id_column))

    logging.info("Parsing GENCODE GTF for gene_id -> gene_name mapping: %s", args.gencode_gtf)
    id_to_name = parse_gtf_gene_id_to_name(gencode_gtf=args.gencode_gtf)

    logging.info("Converting %s Ensembl IDs to gene symbols", df.shape[0])
    out_df, report_df = convert_ids(df=df, id_to_name=id_to_name)

    logging.info("Mapped %s/%s IDs to gene symbols", out_df.shape[0], df.shape[0])
    write_tsv(df=out_df, path=args.out_gene_symbols_tsv)
    write_tsv(df=report_df, path=args.out_mapping_report_tsv)

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
