#!/usr/bin/env python3
"""
Extract high-confidence testis-specific and sperm-retained genes.

This script filters a GTEx testis-specific ranking table annotated with
sperm RNA-seq data to produce a high-confidence gene set.

Filtering criteria are configurable via command-line arguments.

All outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for high-confidence gene extraction.
    """

    input_tsv: str
    output_tsv: str
    min_tau: float
    min_testis_tpm: float
    min_sperm_tpm: float
    log_level: str
    log_path: str | None


def setup_logging(*, log_level: str, log_path: str | None) -> logging.Logger:
    """
    Set up logging.

    Parameters
    ----------
    log_level
        Logging level (INFO, DEBUG).
    log_path
        Optional log file path.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger("extract_testis_sperm")
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
    Parse command-line arguments.

    Returns
    -------
    Config
        Parsed configuration.
    """
    parser = argparse.ArgumentParser(
        description="Extract high-confidence testis-specific and sperm-retained genes."
    )
    parser.add_argument(
        "--input_tsv",
        required=True,
        help="Input TSV annotated with sperm presence (GTEx + sperm).",
    )
    parser.add_argument(
        "--output_tsv",
        required=True,
        help="Output TSV for high-confidence gene set.",
    )
    parser.add_argument(
        "--min_tau",
        type=float,
        default=0.95,
        help="Minimum tau for tissue specificity (default: 0.95).",
    )
    parser.add_argument(
        "--min_testis_tpm",
        type=float,
        default=5.0,
        help="Minimum median TPM in testis (default: 5).",
    )
    parser.add_argument(
        "--min_sperm_tpm",
        type=float,
        default=0.1,
        help="Minimum mean sperm TPM (default: 0.1).",
    )
    parser.add_argument("--log_level", default="INFO", help="Log level.")
    parser.add_argument("--log_path", default=None, help="Optional log file.")

    args = parser.parse_args()

    return Config(
        input_tsv=args.input_tsv,
        output_tsv=args.output_tsv,
        min_tau=args.min_tau,
        min_testis_tpm=args.min_testis_tpm,
        min_sperm_tpm=args.min_sperm_tpm,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run high-confidence gene extraction.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    start = time.time()
    logger.info("Loading input TSV: %s", cfg.input_tsv)
    df = pd.read_csv(cfg.input_tsv, sep="\t", index_col=0, low_memory=False)
    logger.info("Loaded table: %d rows x %d columns", df.shape[0], df.shape[1])

    required_cols = {"tau", "target_median_tpm", "sperm_tpm_mean"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    logger.info(
        "Applying filters: tau >= %.3f, testis TPM >= %.3f, sperm TPM > %.3f",
        cfg.min_tau,
        cfg.min_testis_tpm,
        cfg.min_sperm_tpm,
    )

    filtered = df.loc[
        (df["tau"] >= cfg.min_tau)
        & (df["target_median_tpm"] >= cfg.min_testis_tpm)
        & (df["sperm_tpm_mean"] > cfg.min_sperm_tpm),
        :
    ].copy()

    logger.info("Filtered genes retained: %d", filtered.shape[0])

    filtered = filtered.sort_values(
        by=["tau", "target_median_tpm", "sperm_tpm_mean"],
        ascending=[False, False, False],
    )

    logger.info("Writing output TSV: %s", cfg.output_tsv)
    filtered.to_csv(cfg.output_tsv, sep="\t", index=True)

    elapsed = time.time() - start
    logger.info("Completed in %.1f seconds", elapsed)


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()
