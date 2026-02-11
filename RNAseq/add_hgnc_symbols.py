#!/usr/bin/env python3
"""
Add HGNC-approved gene symbols to a GTEx-ranked table using HGNC mapping.

This script reads:
  - A ranked GTEx gene table (TSV) indexed by Ensembl gene IDs (often versioned)
  - HGNC complete set TSV containing Ensembl gene IDs and approved symbols

It writes a TSV with an added 'hgnc_symbol' column suitable for ClinVar joins.

All outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from typing import Optional

import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for HGNC symbol mapping.

    Attributes
    ----------
    ranked_genes_tsv
        Input ranked GTEx TSV (index must contain Ensembl gene IDs).
    hgnc_tsv
        HGNC complete set TSV (hgnc_complete_set.txt).
    output_tsv
        Output TSV path with added 'hgnc_symbol'.
    description_col
        Column in ranked table that may already contain symbols (default: Description).
    log_level
        Logging verbosity.
    log_path
        Optional log file path.
    """

    ranked_genes_tsv: str
    hgnc_tsv: str
    output_tsv: str
    description_col: str = "Description"
    log_level: str = "INFO"
    log_path: Optional[str] = None


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging.

    Parameters
    ----------
    log_level
        Logging level string.
    log_path
        Optional path to log file.

    Returns
    -------
    logging.Logger
        Logger.
    """
    logger = logging.getLogger("hgnc_mapper")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(fmt)
    sh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.addHandler(sh)

    if log_path:
        fh = logging.FileHandler(log_path)
        fh.setFormatter(fmt)
        fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        logger.addHandler(fh)

    return logger


def parse_args() -> Config:
    """
    Parse CLI arguments.

    Returns
    -------
    Config
        Config instance.
    """
    parser = argparse.ArgumentParser(description="Add HGNC symbols to a GTEx-ranked TSV.")
    parser.add_argument("--ranked_genes_tsv", 
                        required=True, help="Input ranked GTEx TSV.")
    parser.add_argument("--hgnc_tsv", 
                        required=True, help="HGNC complete set TSV (hgnc_complete_set.txt).")
    parser.add_argument("--output_tsv", 
                        required=True, help="Output TSV with added 'hgnc_symbol'.")
    parser.add_argument("--description_col",
                        default="Description",
                        help="Column that may already contain a symbol (default: Description).",
                        )
    parser.add_argument("--log_level", 
                        default="INFO", help="Log level (INFO, DEBUG).")
    parser.add_argument("--log_path", 
                        default=None, help="Optional log file path.")
    args = parser.parse_args()

    return Config(
        ranked_genes_tsv=args.ranked_genes_tsv,
        hgnc_tsv=args.hgnc_tsv,
        output_tsv=args.output_tsv,
        description_col=args.description_col,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def strip_ensembl_version(*, ensembl_id: str) -> str:
    """
    Strip version suffix from Ensembl gene IDs.

    Parameters
    ----------
    ensembl_id
        Ensembl gene ID possibly containing a version (e.g. ENSG... .1).

    Returns
    -------
    str
        Unversioned Ensembl gene ID.
    """
    if ensembl_id is None:
        return ""
    return str(ensembl_id).split(".")[0].strip()


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run HGNC mapping.

    Parameters
    ----------
    cfg
        Config.
    logger
        Logger.
    """
    df = pd.read_csv(cfg.ranked_genes_tsv, sep="\t", index_col=0, low_memory=False)
    logger.info("Loaded ranked table: %d rows x %d cols", df.shape[0], df.shape[1])

    if cfg.description_col not in df.columns:
        raise ValueError(
            f"Column '{cfg.description_col}' not found. Available columns: {list(df.columns)[:30]}"
        )

    df = df.copy()
    df["ensembl_gene_id"] = df.index.map(lambda x: strip_ensembl_version(ensembl_id=x))

    hgnc = pd.read_csv(cfg.hgnc_tsv, sep="\t", low_memory=False)
    logger.info("Loaded HGNC: %d rows x %d cols", hgnc.shape[0], hgnc.shape[1])

    # HGNC columns of interest are typically:
    #  - 'ensembl_gene_id'
    #  - 'symbol' (approved symbol)
    # Some rows may have empty ensembl IDs.
    required = {"ensembl_gene_id", "symbol"}
    if not required.issubset(set(hgnc.columns)):
        raise ValueError(
            "HGNC file missing required columns. "
            "Expected 'ensembl_gene_id' and 'symbol'. "
            f"Available columns: {list(hgnc.columns)[:30]}"
        )



    hgnc_sub = hgnc.loc[
        hgnc["ensembl_gene_id"].astype(str).str.startswith("ENSG", na=False),
        ["ensembl_gene_id", "symbol"],
    ].dropna()

    # Ensure one-to-one mapping: one symbol per Ensembl ID
    hgnc_sub["ensembl_gene_id"] = hgnc_sub["ensembl_gene_id"].astype(str)
    hgnc_sub["symbol"] = hgnc_sub["symbol"].astype(str)

    dup_n = int(hgnc_sub["ensembl_gene_id"].duplicated().sum())
    logger.info("HGNC duplicate Ensembl IDs (rows beyond first): %d", dup_n)


    # Group to guarantee unique Ensembl IDs (take first approved symbol encountered)
    hgnc_map = hgnc_sub.groupby("ensembl_gene_id", sort=False)["symbol"].first()

    logger.info("HGNC mapping entries (unique Ensembl IDs): %d", hgnc_map.shape[0])



    df["hgnc_symbol_from_ensembl"] = df["ensembl_gene_id"].map(hgnc_map).fillna("")

    # Prefer Description if it is a plausible symbol, otherwise use HGNC mapping.
    desc = df[cfg.description_col].astype(str).fillna("")
    desc_is_symbol_like = (~desc.str.startswith("ENSG")) & (desc.str.len() > 0)

    df["hgnc_symbol"] = ""
    df.loc[desc_is_symbol_like, "hgnc_symbol"] = desc.loc[desc_is_symbol_like]
    df.loc[~desc_is_symbol_like, "hgnc_symbol"] = df.loc[~desc_is_symbol_like, "hgnc_symbol_from_ensembl"]

    n_final = int((df["hgnc_symbol"].astype(str).str.len() > 0).sum())
    n_missing = df.shape[0] - n_final
    logger.info("Final symbols assigned: %d / %d (missing: %d)", n_final, df.shape[0], n_missing)

    df.to_csv(cfg.output_tsv, sep="\t", index=True)
    logger.info("Wrote: %s", cfg.output_tsv)


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()

