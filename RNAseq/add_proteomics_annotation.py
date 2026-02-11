#!/usr/bin/env python3
"""
Add proteomics evidence to a gene-level TSV by joining on HGNC gene symbols.

This script reads an Excel proteomics table (e.g. a Proteome Discoverer export)
that contains a 'Gene Symbol' column and one or more per-sample presence columns
(e.g. 'Found in Sample: ...'). It summarises evidence per gene and merges it
onto an existing gene table keyed by HGNC-like symbols (e.g. 'hgnc_symbol_norm').

Outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from dataclasses import dataclass
from typing import List, Optional

import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for proteomics annotation.

    Attributes
    ----------
    input_gene_tsv
        Input TSV to annotate (your master table).
    proteomics_xlsx
        Proteomics Excel file.
    output_tsv
        Output TSV with proteomics columns added.
    gene_symbol_col
        Gene symbol column in the input TSV to join on.
    proteomics_gene_symbol_col
        Gene symbol column in the proteomics Excel file.
    proteomics_accession_col
        UniProt accession column in the proteomics file.
    proteomics_qvalue_col
        Protein-level q-value column in the proteomics file.
    proteomics_fdr_conf_col
        FDR confidence column in the proteomics file.
    require_fdr_confidence
        If set, keep only rows with this confidence (e.g. 'High').
        Set to empty string to disable.
    restrict_species_name
        If set, keep only rows where 'Species Names' matches exactly.
        Set to empty string to disable.
    found_in_sample_prefix
        Prefix used to detect "Found in Sample" columns.
    presence_value_regex
        Regex used to interpret presence columns. If a cell matches this,
        it is counted as present.
    min_unique_peptides
        If > 0, filter out proteins with fewer unique peptides.
    max_accessions_per_gene
        Truncate list of accessions per gene.
    log_level
        Logging level.
    log_path
        Optional log file path.
    """

    input_gene_tsv: str
    proteomics_xlsx: str
    output_tsv: str
    gene_symbol_col: str
    proteomics_gene_symbol_col: str
    proteomics_accession_col: str
    proteomics_qvalue_col: str
    proteomics_fdr_conf_col: str
    require_fdr_confidence: str
    restrict_species_name: str
    found_in_sample_prefix: str
    presence_value_regex: str
    min_unique_peptides: int
    max_accessions_per_gene: int
    log_level: str
    log_path: Optional[str]


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging.

    Parameters
    ----------
    log_level
        Logging level (INFO, DEBUG).
    log_path
        Optional log file path.

    Returns
    -------
    logging.Logger
        Logger.
    """
    logger = logging.getLogger("add_proteomics_annotation")
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
        Parsed configuration.
    """
    p = argparse.ArgumentParser(description="Add proteomics evidence to a gene TSV.")
    p.add_argument("--input_gene_tsv", required=True, help="Input gene TSV (master table).")
    p.add_argument("--proteomics_xlsx", required=True, help="Proteomics Excel file.")
    p.add_argument("--output_tsv", required=True, help="Output TSV with proteomics columns added.")

    p.add_argument(
        "--gene_symbol_col",
        default="hgnc_symbol_norm",
        help="Gene symbol column in input TSV (default: hgnc_symbol_norm).",
    )
    p.add_argument(
        "--proteomics_gene_symbol_col",
        default="Gene Symbol",
        help="Gene symbol column in proteomics file (default: Gene Symbol).",
    )
    p.add_argument(
        "--proteomics_accession_col",
        default="Accession",
        help="Accession column in proteomics file (default: Accession).",
    )
    p.add_argument(
        "--proteomics_qvalue_col",
        default="Exp. q-value: Combined",
        help="Protein q-value column (default: Exp. q-value: Combined).",
    )
    p.add_argument(
        "--proteomics_fdr_conf_col",
        default="Protein FDR Confidence: Combined",
        help="FDR confidence column (default: Protein FDR Confidence: Combined).",
    )
    p.add_argument(
        "--require_fdr_confidence",
        default="High",
        help="Keep only this FDR confidence (default: High). Use empty string to disable.",
    )
    p.add_argument(
        "--restrict_species_name",
        default="Homo sapiens",
        help="Keep only this species (default: Homo sapiens). Use empty string to disable.",
    )
    p.add_argument(
        "--found_in_sample_prefix",
        default="Found in Sample:",
        help="Prefix for sample presence columns (default: Found in Sample:).",
    )
    p.add_argument(
        "--presence_value_regex",
        default=r"^(?i:high|true|yes|1)$",
        help="Regex for values indicating presence (default matches High/TRUE/YES/1).",
    )
    p.add_argument(
        "--min_unique_peptides",
        type=int,
        default=0,
        help="Filter proteins with < this many unique peptides (default: 0).",
    )
    p.add_argument(
        "--max_accessions_per_gene",
        type=int,
        default=20,
        help="Max accessions to keep per gene (default: 20).",
    )
    p.add_argument("--log_level", default="INFO", help="Log level.")
    p.add_argument("--log_path", default=None, help="Optional log file path.")
    a = p.parse_args()

    return Config(
        input_gene_tsv=a.input_gene_tsv,
        proteomics_xlsx=a.proteomics_xlsx,
        output_tsv=a.output_tsv,
        gene_symbol_col=a.gene_symbol_col,
        proteomics_gene_symbol_col=a.proteomics_gene_symbol_col,
        proteomics_accession_col=a.proteomics_accession_col,
        proteomics_qvalue_col=a.proteomics_qvalue_col,
        proteomics_fdr_conf_col=a.proteomics_fdr_conf_col,
        require_fdr_confidence=a.require_fdr_confidence,
        restrict_species_name=a.restrict_species_name,
        found_in_sample_prefix=a.found_in_sample_prefix,
        presence_value_regex=a.presence_value_regex,
        min_unique_peptides=a.min_unique_peptides,
        max_accessions_per_gene=a.max_accessions_per_gene,
        log_level=a.log_level,
        log_path=a.log_path,
    )


def normalise_symbol(*, value: object) -> str:
    """
    Normalise a gene symbol.

    Parameters
    ----------
    value
        Input symbol.

    Returns
    -------
    str
        Uppercased, stripped symbol, or empty string.
    """
    if value is None:
        return ""
    s = str(value).strip().upper()
    if s in {"", "NAN"}:
        return ""
    return s


def dedupe_preserve_order(*, items: List[str], max_items: int) -> List[str]:
    """
    Deduplicate while preserving order, truncated to max_items.

    Parameters
    ----------
    items
        Items.
    max_items
        Max number to keep.

    Returns
    -------
    list
        Deduplicated list
    """
    seen = set()
    out: List[str] = []
    for it in items:
        if it and it not in seen:
            seen.add(it)
            out.append(it)
        if len(out) >= max_items:
            break
    return out


def detect_presence_columns(*, columns: List[str], prefix: str) -> List[str]:
    """
    Detect per-sample presence columns by prefix.

    Parameters
    ----------
    columns
        All column names.
    prefix
        Prefix indicating sample presence columns.

    Returns
    -------
    list
        Presence column names.
    """
    return [c for c in columns if isinstance(c, str) and c.startswith(prefix)]


def parse_numeric_series(*, s: pd.Series) -> pd.Series:
    """
    Convert a series to numeric, coercing errors to NaN.

    Parameters
    ----------
    s
        Series.

    Returns
    -------
    pd.Series
        Numeric series.
    """
    return pd.to_numeric(s, errors="coerce")


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run proteomics annotation and write output.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    t0 = time.time()
    logger.info("Loading input gene TSV: %s", cfg.input_gene_tsv)
    genes = pd.read_csv(cfg.input_gene_tsv, sep="\t", index_col=0, low_memory=False)
    logger.info("Input TSV: %d rows x %d columns", genes.shape[0], genes.shape[1])

    if cfg.gene_symbol_col not in genes.columns:
        raise ValueError(
            f"Missing gene symbol column in input TSV: {cfg.gene_symbol_col}. "
            f"Available columns: {list(genes.columns)[:30]}"
        )

    genes = genes.copy()
    genes["gene_symbol_norm"] = genes[cfg.gene_symbol_col].map(lambda x: normalise_symbol(value=x))

    logger.info("Loading proteomics Excel: %s", cfg.proteomics_xlsx)
    prot = pd.read_excel(cfg.proteomics_xlsx, sheet_name=0)
    logger.info("Proteomics table: %d rows x %d columns", prot.shape[0], prot.shape[1])

    required = [
        cfg.proteomics_gene_symbol_col,
        cfg.proteomics_accession_col,
        cfg.proteomics_qvalue_col,
        cfg.proteomics_fdr_conf_col,
    ]
    for col in required:
        if col not in prot.columns:
            raise ValueError(
                f"Missing required column in proteomics file: {col}. "
                f"Available columns: {list(prot.columns)[:40]}"
            )

    prot = prot.copy()
    prot["gene_symbol_norm"] = prot[cfg.proteomics_gene_symbol_col].map(lambda x: normalise_symbol(value=x))
    prot["accession_str"] = prot[cfg.proteomics_accession_col].astype(str).str.strip()

    # Optional filters
    if cfg.require_fdr_confidence.strip():
        before = prot.shape[0]
        prot = prot.loc[prot[cfg.proteomics_fdr_conf_col].astype(str).str.strip() == cfg.require_fdr_confidence, :].copy()
        logger.info("Filter FDR confidence '%s': %d -> %d", cfg.require_fdr_confidence, before, prot.shape[0])

    if cfg.restrict_species_name.strip() and "Species Names" in prot.columns:
        before = prot.shape[0]
        prot = prot.loc[prot["Species Names"].astype(str).str.strip() == cfg.restrict_species_name, :].copy()
        logger.info("Filter species '%s': %d -> %d", cfg.restrict_species_name, before, prot.shape[0])

    if cfg.min_unique_peptides > 0 and "# Unique Peptides" in prot.columns:
        before = prot.shape[0]
        up = parse_numeric_series(s=prot["# Unique Peptides"])
        prot = prot.loc[up.fillna(0) >= cfg.min_unique_peptides, :].copy()
        logger.info("Filter min unique peptides >= %d: %d -> %d", cfg.min_unique_peptides, before, prot.shape[0])

    prot = prot.loc[prot["gene_symbol_norm"] != "", :].copy()
    logger.info("Proteins with valid gene symbol: %d", prot.shape[0])

    presence_cols = detect_presence_columns(columns=list(prot.columns), prefix=cfg.found_in_sample_prefix)
    logger.info("Detected presence columns: %d", len(presence_cols))

    present_re = re.compile(cfg.presence_value_regex)

    def summarise_gene(df: pd.DataFrame) -> pd.Series:
        """
        Summarise proteomics evidence for one gene.

        Parameters
        ----------
        df
            Per-gene proteomics rows.

        Returns
        -------
        pd.Series
            Per-gene summary fields.
        """
        acc = dedupe_preserve_order(
            items=[str(x).strip() for x in df["accession_str"].tolist() if str(x).strip()],
            max_items=cfg.max_accessions_per_gene,
        )

        qv = parse_numeric_series(s=df[cfg.proteomics_qvalue_col])
        best_q = float(qv.min()) if qv.notna().any() else float("nan")

        pep = parse_numeric_series(s=df["# Peptides"]) if "# Peptides" in df.columns else pd.Series([float("nan")] * df.shape[0])
        uniq = parse_numeric_series(s=df["# Unique Peptides"]) if "# Unique Peptides" in df.columns else pd.Series([float("nan")] * df.shape[0])
        cov = parse_numeric_series(s=df["Coverage [%]"]) if "Coverage [%]" in df.columns else pd.Series([float("nan")] * df.shape[0])

        pep_max = float(pep.max()) if pep.notna().any() else float("nan")
        uniq_max = float(uniq.max()) if uniq.notna().any() else float("nan")
        cov_max = float(cov.max()) if cov.notna().any() else float("nan")

        # Per-sample presence summary
        n_samples = len(presence_cols)
        n_present = 0
        if n_samples > 0:
            any_present_by_sample: List[bool] = []
            for c in presence_cols:
                vals = df[c].astype(str).str.strip()
                is_present = vals.map(lambda v: bool(present_re.match(v)))
                any_present_by_sample.append(bool(is_present.any()))
            n_present = int(sum(any_present_by_sample))

        return pd.Series(
            {
                "prot_accessions": ";".join(acc),
                "prot_n_accessions": len(acc),
                "prot_best_qvalue": best_q,
                "prot_peptides_max": pep_max,
                "prot_unique_peptides_max": uniq_max,
                "prot_coverage_pct_max": cov_max,
                "prot_samples_detected": n_present,
                "prot_samples_total": n_samples,
                "prot_present_any": (n_present > 0) if n_samples > 0 else True,
                "prot_present_fraction": (n_present / n_samples) if n_samples > 0 else float("nan"),
            }
        )

    logger.info("Summarising proteomics per gene")
    prot_ann = prot.groupby("gene_symbol_norm", sort=False).apply(summarise_gene)
    prot_ann.index.name = "gene_symbol_norm"
    logger.info("Proteomics per-gene rows: %d", prot_ann.shape[0])

    logger.info("Merging proteomics onto input table")
    out = genes.join(prot_ann, on="gene_symbol_norm")

    # Fill missing values
    for col in ["prot_n_accessions", "prot_samples_detected", "prot_samples_total"]:
        if col in out.columns:
            out[col] = out[col].fillna(0).astype(int)

    for col in ["prot_accessions"]:
        if col in out.columns:
            out[col] = out[col].fillna("").astype(str)

    for col in ["prot_present_any"]:
        if col in out.columns:
            out[col] = out[col].fillna(False).astype(bool)

    for col in ["prot_best_qvalue", "prot_peptides_max", "prot_unique_peptides_max", "prot_coverage_pct_max", "prot_present_fraction"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    n_present = int(out["prot_present_any"].sum()) if "prot_present_any" in out.columns else 0
    logger.info("Genes with proteomics presence: %d / %d", n_present, out.shape[0])

    logger.info("Writing output TSV: %s", cfg.output_tsv)
    out.to_csv(cfg.output_tsv, sep="\t", index=True)

    logger.info("Done in %.1f seconds", time.time() - t0)


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()
