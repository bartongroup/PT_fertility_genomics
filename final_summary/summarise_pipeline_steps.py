#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_pipeline_steps.py

Summarise per-step outputs for a multi-stage gene prioritisation pipeline.

This script scans a base directory for step folders (e.g. 01_hpo, 02_...),
finds TSV files within each step, and reports basic dataset diagnostics:
- number of rows
- detected gene identifier column (heuristic)
- non-empty gene identifier count
- unique gene identifier count

Optionally, it can also summarise a nominated membership/flag column
(e.g. 'proteomics_present_any_source') as True/False counts if present.

Outputs are tab-separated (TSV).

Example
-------
python summarise_pipeline_steps.py \
  --base_dir /home/pthorpe001/data/2026_sperm_Gates/results \
  --out_dir /home/pthorpe001/data/2026_sperm_Gates/results/step_summaries \
  --step_dir_regex '^(0[1-9]_.*|testis_.*)$' \
  --membership_columns proteomics_present_any_source sperm_rnaseq_present \
  --verbose
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


@dataclass(frozen=True)
class FileSummary:
    """Container for per-file summary statistics."""
    step: str
    tsv_path: str
    rows: int
    gene_column: str
    gene_non_empty: int
    gene_unique: int
    extra_counts: str


def setup_logger(verbose: bool) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose
        If True, enable DEBUG logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_bool_tokens(series: pd.Series) -> pd.Series:
    """
    Normalise common boolean tokens.

    Parameters
    ----------
    series
        Input series.

    Returns
    -------
    pd.Series
        Lowercased string series.
    """
    return series.astype(str).str.strip().str.lower()


def detect_gene_column(df: pd.DataFrame) -> str:
    """
    Heuristically detect the most likely gene identifier column.

    Preference order:
    - gene_key
    - hgnc_symbol
    - gene_symbol
    - approvedSymbol / approved_symbol
    - GeneSymbol
    - ensembl_gene_id
    - gene_id

    If none found, returns empty string.

    Parameters
    ----------
    df
        Input DataFrame.

    Returns
    -------
    str
        Detected column name or empty string.
    """
    candidates = [
        "gene_key",
        "hgnc_symbol",
        "gene_symbol",
        "approvedSymbol",
        "approved_symbol",
        "GeneSymbol",
        "ensembl_gene_id",
        "gene_id",
    ]
    cols = {str(c): c for c in df.columns}
    for want in candidates:
        if want in cols:
            return cols[want]
    return ""


def safe_read_tsv(path: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """
    Read a TSV robustly.

    Parameters
    ----------
    path
        TSV path.
    nrows
        Optional row limit.

    Returns
    -------
    pd.DataFrame
        DataFrame with dtype=str.
    """
    return pd.read_csv(path, sep="\t", dtype=str, nrows=nrows).fillna("")


def summarise_file(
    step: str,
    tsv_path: str,
    membership_columns: Sequence[str],
) -> FileSummary:
    """
    Summarise a single TSV file.

    Parameters
    ----------
    step
        Step folder name.
    tsv_path
        TSV file path.
    membership_columns
        Optional list of boolean-like columns to count.

    Returns
    -------
    FileSummary
        Summary object.
    """
    df = safe_read_tsv(path=tsv_path)
    gene_col = detect_gene_column(df=df)

    rows = int(df.shape[0])
    gene_non_empty = 0
    gene_unique = 0

    if gene_col:
        s = df[gene_col].astype(str).str.strip()
        gene_non_empty = int((s != "").sum())
        gene_unique = int(s[s != ""].nunique())
    else:
        logging.debug("No gene column detected in: %s", tsv_path)

    extras: List[str] = []
    for col in membership_columns:
        if col in df.columns:
            v = parse_bool_tokens(df[col]).value_counts(dropna=False).to_dict()
            true_n = int(v.get("true", 0) + v.get("1", 0) + v.get("t", 0) + v.get("yes", 0) + v.get("y", 0))
            false_n = int(v.get("false", 0) + v.get("0", 0) + v.get("f", 0) + v.get("no", 0) + v.get("n", 0))
            empty_n = int(v.get("", 0))
            extras.append(f"{col}:true={true_n}|false={false_n}|empty={empty_n}")

    extra_counts = ";".join(extras)
    return FileSummary(
        step=step,
        tsv_path=tsv_path,
        rows=rows,
        gene_column=str(gene_col),
        gene_non_empty=gene_non_empty,
        gene_unique=gene_unique,
        extra_counts=extra_counts,
    )


def find_step_dirs(base_dir: str, step_dir_regex: str) -> List[str]:
    """
    Find step directories under base_dir.

    Parameters
    ----------
    base_dir
        Base directory containing step folders.
    step_dir_regex
        Regex to match step folder names.

    Returns
    -------
    list of str
        Absolute paths to step directories.
    """
    patt = re.compile(step_dir_regex)
    out: List[str] = []
    for name in sorted(os.listdir(base_dir)):
        path = os.path.join(base_dir, name)
        if os.path.isdir(path) and patt.search(name):
            out.append(path)
    return out


def find_tsvs(step_dir: str) -> List[str]:
    """
    Recursively find TSV files under a step directory.

    Parameters
    ----------
    step_dir
        Step directory.

    Returns
    -------
    list of str
        TSV paths.
    """
    tsvs: List[str] = []
    for root, _, files in os.walk(step_dir):
        for fn in files:
            if fn.lower().endswith(".tsv"):
                tsvs.append(os.path.join(root, fn))
    return sorted(tsvs)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse CLI arguments.

    Returns
    -------
    argparse.Namespace
        Parsed args.
    """
    p = argparse.ArgumentParser(description="Summarise stepwise pipeline TSV outputs.")
    p.add_argument("--base_dir", required=True, help="Base directory containing step folders.")
    p.add_argument("--out_dir", required=True, help="Output directory for summary TSVs.")
    p.add_argument(
        "--step_dir_regex",
        default=r"^(0[1-9]_.*|testis_.*)$",
        help="Regex matching step folder names.",
    )
    p.add_argument(
        "--membership_columns",
        nargs="*",
        default=[],
        help="Optional boolean-like columns to summarise if present.",
    )
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the step summary workflow."""
    args = parse_args(argv=argv)
    setup_logger(verbose=bool(args.verbose))

    if not os.path.isdir(args.base_dir):
        raise SystemExit(f"base_dir not found: {args.base_dir}")

    os.makedirs(args.out_dir, exist_ok=True)

    start = datetime.now()
    logging.info("Base dir: %s", args.base_dir)
    logging.info("Output dir: %s", args.out_dir)

    step_dirs = find_step_dirs(base_dir=args.base_dir, step_dir_regex=args.step_dir_regex)
    logging.info("Detected %s step directories", len(step_dirs))

    file_summaries: List[FileSummary] = []
    step_rollups: Dict[str, Dict[str, int]] = {}

    for step_path in step_dirs:
        step_name = os.path.basename(step_path)
        tsvs = find_tsvs(step_dir=step_path)
        logging.info("Step %s: %s TSV files", step_name, len(tsvs))

        step_rows = 0
        step_gene_unique_max = 0

        for tsv_path in tsvs:
            summ = summarise_file(
                step=step_name,
                tsv_path=tsv_path,
                membership_columns=args.membership_columns,
            )
            file_summaries.append(summ)

            step_rows += summ.rows
            step_gene_unique_max = max(step_gene_unique_max, summ.gene_unique)

        step_rollups[step_name] = {
            "n_tsv_files": len(tsvs),
            "sum_rows_across_tsvs": step_rows,
            "max_unique_genes_in_any_file": step_gene_unique_max,
        }

    file_df = pd.DataFrame([s.__dict__ for s in file_summaries])
    step_df = (
        pd.DataFrame.from_dict(step_rollups, orient="index")
        .reset_index()
        .rename(columns={"index": "step"})
        .sort_values(by=["step"], ascending=[True])
    )

    file_out = os.path.join(args.out_dir, "file_summaries.tsv")
    step_out = os.path.join(args.out_dir, "step_summaries.tsv")

    file_df.to_csv(file_out, sep="\t", index=False)
    step_df.to_csv(step_out, sep="\t", index=False)

    elapsed = datetime.now() - start
    logging.info("Wrote: %s", file_out)
    logging.info("Wrote: %s", step_out)
    logging.info("Done. Elapsed time: %s", str(elapsed).split(".")[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
