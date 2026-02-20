#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_gene_universe_from_gtf.py

Extract all gene symbols from a GENCODE GRCh38 GTF and create a universe
gene list suitable for matched-set enrichment analysis.

This script:

- Reads a GENCODE GTF (gzipped supported)
- Extracts gene_name from "gene" features
- Keeps only standard chromosomes (chr1-22, chrX, chrY)
- Outputs a TSV containing a single column: gene_key

This universe list should be used to generate the genome-wide feature
table before running matched_set_enrichment.py.

Output
------
TSV with column:

gene_key

Example
-------
python build_gene_universe_from_gtf.py \
  --gencode_gtf gencode.v49.primary_assembly.annotation.gtf.gz \
  --out_tsv gencode_v49_gene_keys_universe.tsv
"""

from __future__ import annotations

import argparse
import gzip
import logging
from typing import Dict, Optional, Sequence, Set

import pandas as pd


STD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def setup_logger(verbose: bool) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def open_maybe_gzip(path: str):
    """Open plain or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def normalise_chrom(chrom: str) -> str:
    """
    Convert chromosome names to UCSC-style chr format.

    Examples:
        1  -> chr1
        chr1 -> chr1
    """
    c = str(chrom).strip()
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def parse_gtf_gene_universe(gtf_path: str) -> Set[str]:
    """
    Extract gene_name values from GTF gene features.

    Parameters
    ----------
    gtf_path : str
        Path to GENCODE GTF.

    Returns
    -------
    Set[str]
        Unique gene symbols on standard chromosomes.
    """
    gene_names: Set[str] = set()
    n_lines = 0
    n_genes = 0

    with open_maybe_gzip(gtf_path) as handle:
        for line in handle:
            n_lines += 1
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _, feature, *_rest = parts[:4]
            attrs = parts[8]

            if feature != "gene":
                continue

            chrom_norm = normalise_chrom(chrom)
            if chrom_norm not in STD_CHROMS:
                continue

            gene_name = None
            for item in attrs.split(";"):
                item = item.strip()
                if item.startswith("gene_name"):
                    gene_name = item.split(" ", 1)[1].strip().strip('"')
                    break

            if gene_name:
                gene_names.add(gene_name)
                n_genes += 1

            if n_genes % 10000 == 0:
                logging.info("Parsed %s gene records...", n_genes)

    logging.info(
        "Finished parsing: %s gene records retained from %s lines",
        len(gene_names),
        n_lines,
    )
    return gene_names


def write_universe(gene_names: Set[str], out_path: str) -> None:
    """Write gene universe TSV."""
    df = pd.DataFrame({"gene_key": sorted(gene_names)})
    df.to_csv(out_path, sep="\t", index=False)
    logging.info("Wrote universe gene list: %s genes -> %s", df.shape[0], out_path)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Build gene universe from GENCODE GTF"
    )
    parser.add_argument("--gencode_gtf", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run script."""
    args = parse_args(argv)
    setup_logger(verbose=args.verbose)

    logging.info("Parsing GTF: %s", args.gencode_gtf)
    gene_names = parse_gtf_gene_universe(gtf_path=args.gencode_gtf)

    write_universe(gene_names=gene_names, out_path=args.out_tsv)

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
