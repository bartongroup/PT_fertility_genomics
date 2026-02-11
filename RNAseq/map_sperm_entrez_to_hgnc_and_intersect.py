#!/usr/bin/env python3
"""
Map Entrez Gene IDs in sperm TPM tables to HGNC symbols and intersect with GTEx ranking.

Inputs
------
1) Sperm TPM table from GEO (Entrez 'GeneID' as first column).
   This may be a single-sample file (one GSM column) or multi-sample matrix.

2) HGNC complete set TSV (hgnc_complete_set.txt), containing:
   - entrez_id
   - symbol

3) GTEx ranked table (with hgnc_symbol column), optional but recommended.

Outputs
-------
A) Mapped sperm table with:
   - entrez_id
   - hgnc_symbol
   - original TPM columns

B) (Optional) GTEx ranked table annotated with sperm presence:
   - sperm_present_any: bool (TPM > threshold in at least one sample)
   - sperm_present_frac: fraction of sperm samples with TPM > threshold
   - sperm_tpm_mean: mean TPM across sperm samples
   - sperm_tpm_median: median TPM across sperm samples

All outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for mapping and intersection.

    Attributes
    ----------
    sperm_tpm_path
        Path to sperm TPM table (TSV or TSV.GZ).
    hgnc_tsv
        Path to hgnc_complete_set.txt (TSV).
    gtex_ranked_with_hgnc_tsv
        Optional path to GTEx ranked table containing 'hgnc_symbol'.
    output_sperm_mapped_tsv
        Output TSV for mapped sperm table.
    output_gtex_with_sperm_tsv
        Optional output TSV for GTEx table annotated with sperm presence.
    sperm_tpm_threshold
        TPM threshold to call a gene present in sperm.
    log_level
        Logging verbosity.
    log_path
        Optional log file path.
    """

    sperm_tpm_path: str
    hgnc_tsv: str
    gtex_ranked_with_hgnc_tsv: Optional[str]
    output_sperm_mapped_tsv: str
    output_gtex_with_sperm_tsv: Optional[str]
    sperm_tpm_threshold: float
    log_level: str
    log_path: Optional[str]


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging to stdout and optionally a file.

    Parameters
    ----------
    log_level
        Logging level string (INFO, DEBUG).
    log_path
        Optional path to log file.

    Returns
    -------
    logging.Logger
        Logger instance.
    """
    logger = logging.getLogger("sperm_entrez_to_hgnc")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.handlers = []
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(stream=sys.stdout)
    sh.setFormatter(formatter)
    sh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    logger.addHandler(sh)

    if log_path:
        fh = logging.FileHandler(log_path)
        fh.setFormatter(formatter)
        fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        logger.addHandler(fh)

    return logger


def timed(*, logger: logging.Logger, label: str):
    """
    Context manager for timing blocks.

    Parameters
    ----------
    logger
        Logger instance.
    label
        Block label.
    """

    class _Timer:
        def __enter__(self):
            self.start = time.time()
            logger.info("%s: started", label)
            return self

        def __exit__(self, exc_type, exc, tb):
            elapsed = time.time() - self.start
            if exc_type is None:
                logger.info("%s: completed in %.1f seconds", label, elapsed)
            else:
                logger.error("%s: failed after %.1f seconds", label, elapsed)
            return False

    return _Timer()


def parse_args() -> Config:
    """
    Parse command-line arguments.

    Returns
    -------
    Config
        Parsed configuration object.
    """
    parser = argparse.ArgumentParser(
        description="Map sperm Entrez Gene IDs to HGNC symbols and intersect with GTEx ranked table."
    )
    parser.add_argument(
        "--sperm_tpm_path",
        required=True,
        help="Path to sperm TPM table (TSV or TSV.GZ) with first column 'GeneID'.",
    )
    parser.add_argument(
        "--hgnc_tsv",
        required=True,
        help="Path to HGNC complete set TSV (hgnc_complete_set.txt).",
    )
    parser.add_argument(
        "--gtex_ranked_with_hgnc_tsv",
        default=None,
        help="Optional GTEx ranked TSV containing 'hgnc_symbol' column.",
    )
    parser.add_argument(
        "--output_sperm_mapped_tsv",
        required=True,
        help="Output TSV for sperm TPM table with HGNC symbol mapping.",
    )
    parser.add_argument(
        "--output_gtex_with_sperm_tsv",
        default=None,
        help="Optional output TSV for GTEx ranked table annotated with sperm presence.",
    )
    parser.add_argument(
        "--sperm_tpm_threshold",
        type=float,
        default=0.1,
        help="TPM threshold for calling sperm presence (default: 0.1).",
    )
    parser.add_argument("--log_level", default="INFO", help="Log level (INFO, DEBUG).")
    parser.add_argument("--log_path", default=None, help="Optional log file path.")

    args = parser.parse_args()

    return Config(
        sperm_tpm_path=args.sperm_tpm_path,
        hgnc_tsv=args.hgnc_tsv,
        gtex_ranked_with_hgnc_tsv=args.gtex_ranked_with_hgnc_tsv,
        output_sperm_mapped_tsv=args.output_sperm_mapped_tsv,
        output_gtex_with_sperm_tsv=args.output_gtex_with_sperm_tsv,
        sperm_tpm_threshold=args.sperm_tpm_threshold,
        log_level=args.log_level,
        log_path=args.log_path,
    )


def load_hgnc_entrez_map(*, hgnc_tsv: str, logger: logging.Logger) -> pd.Series:
    """
    Load HGNC table and build an Entrez ID -> HGNC symbol mapping.

    Parameters
    ----------
    hgnc_tsv
        Path to hgnc_complete_set.txt.
    logger
        Logger instance.

    Returns
    -------
    pd.Series
        Series indexed by entrez_id (as string) with values as HGNC symbol.
    """
    with timed(logger=logger, label="Loading HGNC and building Entrez->symbol map"):
        hgnc = pd.read_csv(hgnc_tsv, sep="\t", low_memory=False)
        required = {"entrez_id", "symbol"}
        if not required.issubset(set(hgnc.columns)):
            raise ValueError(
                "HGNC file missing required columns. "
                "Expected 'entrez_id' and 'symbol'. "
                f"Available columns: {list(hgnc.columns)[:30]}"
            )

        sub = hgnc.loc[:, ["entrez_id", "symbol"]].copy()
        sub["entrez_id"] = sub["entrez_id"].map(lambda x: normalise_entrez_id(value=x))

        sub = sub.loc[sub["entrez_id"] != "", :]

        sub["symbol"] = sub["symbol"].astype(str).str.strip()

        sub = sub.loc[(sub["entrez_id"] != "") & (sub["entrez_id"] != "nan") & (sub["symbol"] != ""), :]
        dup_n = int(sub["entrez_id"].duplicated().sum())
        logger.info("HGNC rows with Entrez IDs: %d (duplicate Entrez IDs beyond first: %d)", sub.shape[0], dup_n)

        entrez_to_symbol = sub.groupby("entrez_id", sort=False)["symbol"].first()
        logger.info("Unique Entrez IDs mapped: %d", entrez_to_symbol.shape[0])

        return entrez_to_symbol


def normalise_entrez_id(*, value: object) -> str:
    """
    Normalise Entrez IDs into a comparable string form.

    This removes common formatting issues such as trailing '.0' caused by
    float parsing, and strips whitespace.

    Parameters
    ----------
    value
        Input value representing an Entrez ID.

    Returns
    -------
    str
        Normalised Entrez ID as a string, or empty string if missing/invalid.
    """
    if value is None:
        return ""
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return ""
    # Handle float-like strings, e.g. '653635.0'
    if s.endswith(".0"):
        s = s[:-2]
    # Ensure it is numeric
    if not s.isdigit():
        return ""
    return s



def load_sperm_tpm(*, sperm_tpm_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load sperm TPM table.

    Parameters
    ----------
    sperm_tpm_path
        Path to sperm TPM TSV or TSV.GZ.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by entrez_id (string), columns are sample(s) TPM.
    """
    with timed(logger=logger, label="Loading sperm TPM table"):
        df = pd.read_csv(sperm_tpm_path, sep="\t", low_memory=False)
        logger.info("Loaded sperm table: %d rows x %d columns", df.shape[0], df.shape[1])
        logger.info("First columns: %s", list(df.columns)[:10])

        if df.shape[1] < 2:
            raise ValueError("Sperm TPM table must have at least two columns: GeneID + one sample.")

        gene_col = df.columns[0]
        if str(gene_col).lower() not in {"geneid", "gene_id", "entrez", "entrez_id"}:
            logger.warning("First column is '%s' (expected 'GeneID'). Proceeding anyway.", gene_col)

        df = df.rename(columns={gene_col: "entrez_id"}).copy()

        df["entrez_id"] = df["entrez_id"].map(lambda x: normalise_entrez_id(value=x))
        df = df.loc[df["entrez_id"] != "", :].copy()
        logger.info("Sperm genes with valid Entrez IDs: %d", df.shape[0])

        tpm_cols = [c for c in df.columns if c != "entrez_id"]
        logger.info("TPM columns detected: %d", len(tpm_cols))

        for c in tpm_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

        df = df.set_index("entrez_id", drop=True)
        df = df.loc[df.index != "", :]
        logger.info("Sperm TPM matrix: %d genes x %d samples", df.shape[0], df.shape[1])

        return df


def map_sperm_to_hgnc(
    *,
    sperm_tpm: pd.DataFrame,
    entrez_to_symbol: pd.Series,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Map sperm Entrez IDs to HGNC symbols and return a mapped table.

    Parameters
    ----------
    sperm_tpm
        Sperm TPM matrix indexed by entrez_id.
    entrez_to_symbol
        Series mapping entrez_id to HGNC symbol.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Mapped sperm table with columns: entrez_id, hgnc_symbol, TPM columns...
    """
    with timed(logger=logger, label="Mapping sperm Entrez IDs to HGNC symbols"):
        mapped = sperm_tpm.copy()
        mapped["hgnc_symbol"] = mapped.index.map(lambda x: entrez_to_symbol.get(str(x), ""))
        mapped["hgnc_symbol"] = mapped["hgnc_symbol"].fillna("").astype(str).str.strip()
        n_mapped = int((mapped["hgnc_symbol"].str.len() > 0).sum())



        logger.info("Mapped HGNC symbols for sperm genes: %d / %d", n_mapped, mapped.shape[0])

        mapped = mapped.reset_index().rename(columns={"index": "entrez_id"})
        cols = ["entrez_id", "hgnc_symbol"] + [c for c in mapped.columns if c not in {"entrez_id", "hgnc_symbol"}]
        mapped = mapped.loc[:, cols]

        return mapped


def summarise_sperm_expression(
    *,
    sperm_mapped: pd.DataFrame,
    tpm_threshold: float,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Summarise sperm expression per HGNC symbol across all sperm samples.

    Parameters
    ----------
    sperm_mapped
        Sperm table with 'hgnc_symbol' and TPM columns.
    tpm_threshold
        Threshold for presence calling.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        Per-symbol summary indexed by hgnc_symbol with columns:
          - sperm_present_any
          - sperm_present_frac
          - sperm_tpm_mean
          - sperm_tpm_median
    """
    with timed(logger=logger, label="Summarising sperm expression by HGNC symbol"):
        if "hgnc_symbol" not in sperm_mapped.columns:
            raise ValueError("Expected 'hgnc_symbol' column in sperm_mapped table.")

        tpm_cols = [c for c in sperm_mapped.columns if c not in {"entrez_id", "hgnc_symbol"}]
        if not tpm_cols:
            raise ValueError("No TPM columns found in sperm_mapped table.")

        df = sperm_mapped.loc[sperm_mapped["hgnc_symbol"].astype(str).str.len() > 0, :].copy()
        for c in tpm_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

        # If multiple Entrez IDs map to the same symbol, aggregate by mean TPM per sample.
        agg = df.groupby("hgnc_symbol", sort=False)[tpm_cols].mean()

        present_any = (agg > tpm_threshold).any(axis=1)
        present_frac = (agg > tpm_threshold).sum(axis=1) / float(len(tpm_cols))
        tpm_mean = agg.mean(axis=1)
        tpm_median = agg.median(axis=1)

        out = pd.DataFrame(
            {
                "sperm_present_any": present_any.astype(bool),
                "sperm_present_frac": present_frac.astype(float),
                "sperm_tpm_mean": tpm_mean.astype(float),
                "sperm_tpm_median": tpm_median.astype(float),
            }
        )

        logger.info("Per-symbol sperm summary rows: %d", out.shape[0])
        logger.info(
            "Symbols present in sperm (any sample, TPM>%.3g): %d",
            tpm_threshold,
            int(out["sperm_present_any"].sum()),
        )

        return out


def annotate_gtex_with_sperm(
    *,
    gtex_ranked_path: str,
    sperm_summary: pd.DataFrame,
    output_path: str,
    logger: logging.Logger,
) -> None:
    """
    Add sperm summary columns onto the GTEx ranked table using hgnc_symbol.

    Parameters
    ----------
    gtex_ranked_path
        Path to GTEx ranked table containing 'hgnc_symbol'.
    sperm_summary
        DataFrame indexed by hgnc_symbol with sperm summary columns.
    output_path
        Output TSV path.
    logger
        Logger instance.
    """
    with timed(logger=logger, label="Annotating GTEx ranked table with sperm summary"):
        gtex = pd.read_csv(gtex_ranked_path, sep="\t", index_col=0, low_memory=False)
        if "hgnc_symbol" not in gtex.columns:
            raise ValueError("GTEx ranked table must contain 'hgnc_symbol' column.")

        gtex = gtex.copy()
        gtex["hgnc_symbol_norm"] = gtex["hgnc_symbol"].astype(str).str.upper().str.strip()

        ss = sperm_summary.copy()
        ss.index = ss.index.astype(str).str.upper().str.strip()

        out = gtex.join(ss, on="hgnc_symbol_norm")

        # Fill missing sperm stats for genes not in sperm list
        if "sperm_present_any" in out.columns:
            out["sperm_present_any"] = out["sperm_present_any"].fillna(False).astype(bool)
        for col in ["sperm_present_frac", "sperm_tpm_mean", "sperm_tpm_median"]:
            if col in out.columns:
                out[col] = out[col].fillna(0.0).astype(float)

        n_present = int(out["sperm_present_any"].sum()) if "sperm_present_any" in out.columns else 0
        logger.info("GTEx genes annotated as present in sperm: %d / %d", n_present, out.shape[0])

        out.to_csv(output_path, sep="\t", index=True)
        logger.info("Wrote: %s", output_path)


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run mapping and optional intersection.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    entrez_to_symbol = load_hgnc_entrez_map(hgnc_tsv=cfg.hgnc_tsv, logger=logger)
    sperm_tpm = load_sperm_tpm(sperm_tpm_path=cfg.sperm_tpm_path, logger=logger)
    sperm_mapped = map_sperm_to_hgnc(sperm_tpm=sperm_tpm, entrez_to_symbol=entrez_to_symbol, logger=logger)

    with timed(logger=logger, label="Writing mapped sperm table"):
        sperm_mapped.to_csv(cfg.output_sperm_mapped_tsv, sep="\t", index=False)
        logger.info("Wrote: %s", cfg.output_sperm_mapped_tsv)

    # Optional: annotate GTEx ranked table
    if cfg.gtex_ranked_with_hgnc_tsv and cfg.output_gtex_with_sperm_tsv:
        sperm_summary = summarise_sperm_expression(
            sperm_mapped=sperm_mapped,
            tpm_threshold=cfg.sperm_tpm_threshold,
            logger=logger,
        )
        annotate_gtex_with_sperm(
            gtex_ranked_path=cfg.gtex_ranked_with_hgnc_tsv,
            sperm_summary=sperm_summary,
            output_path=cfg.output_gtex_with_sperm_tsv,
            logger=logger,
        )
    else:
        logger.info(
            "Skipping GTEx annotation (provide both --gtex_ranked_with_hgnc_tsv and --output_gtex_with_sperm_tsv)."
        )


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()
