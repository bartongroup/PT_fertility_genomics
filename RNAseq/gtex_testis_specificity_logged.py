#!/usr/bin/env python3
"""
GTEx v11 testis specificity scoring from gene TPM GCT with detailed logging.

This script:
  1) Loads GTEx gene TPM in GCT format.
  2) Loads GTEx SampleAttributesDS to map sample IDs to tissues (SMTSD by default).
  3) Computes gene × tissue median TPM.
  4) Scores target-tissue (default: Testis) specificity:
       - Tau tissue-specificity score (0..1)
       - Log2 enrichment of target median vs maximum non-target tissue median
       - Target presence fraction among target samples (TPM threshold)
       (filter to decide whether a gene is genuinely expressed in the target tissue, rather than appearing due to noise or a few outlier samples)
       (In how many target samples does this gene pass the TPM threshold?)
       - Top-K tissues by median TPM for quick inspection
  5) Writes tab-separated outputs (TSV).

Outputs
-------
1) Ranked per-gene summary TSV
2) Gene × tissue median TPM matrix TSV

All outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
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
    Configuration parameters for specificity scoring.

    Attributes
    ----------
    gct_path
        Path to GTEx gene TPM GCT (.gct or .gct.gz).
    sample_attributes_path
        Path to GTEx SampleAttributesDS.txt.
    target_tissue
        Tissue label to score as target (default: 'Testis').
    sample_id_col
        Column name containing sample IDs (default: 'SAMPID').
    tissue_col
        Column name containing detailed tissue labels (default: 'SMTSD').
    min_tpm_present
        TPM threshold for calling a gene 'present' in a sample.
    min_target_present_fraction
        Minimum fraction of target samples where TPM >= min_tpm_present to keep gene.
        Use None to disable filtering.
    top_k_tissues
        Number of top tissues to store per gene (by median TPM).
    out_ranked_tsv
        Output TSV path for ranked gene summary.
    out_tissue_medians_tsv
        Output TSV path for gene × tissue median TPM matrix.
    log_path
        Optional path to a log file.
    log_level
        Logging level string.
    epsilon
        Small constant to avoid division by zero.
    """

    gct_path: str
    sample_attributes_path: str
    target_tissue: str
    sample_id_col: str
    tissue_col: str
    min_tpm_present: float
    min_target_present_fraction: Optional[float]
    top_k_tissues: int
    out_ranked_tsv: str
    out_tissue_medians_tsv: str
    log_path: Optional[str]
    log_level: str
    epsilon: float = 1e-8


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging to stdout and optionally a file.

    Parameters
    ----------
    log_level
        Logging level (e.g. 'INFO', 'DEBUG').
    log_path
        Optional log file path.

    Returns
    -------
    logging.Logger
        Configured logger.
    """
    logger = logging.getLogger("gtex_testis_specificity")
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


def parse_args() -> Config:
    """
    Parse command line arguments.

    Returns
    -------
    Config
        Parsed configuration object.
    """
    parser = argparse.ArgumentParser(
        description="Compute GTEx v11 testis tissue-specificity scores from gene TPM GCT (with logging)."
    )
    parser.add_argument("--gct_path", required=True, 
                        help="Path to GTEx gene TPM GCT (.gct or .gct.gz).")
    parser.add_argument(
        "--sample_attributes_path",
        required=True,
        help="Path to GTEx SampleAttributesDS.txt (tab-delimited).",
    )
    parser.add_argument("--target_tissue", default="Testis", 
                        help="Target tissue name to score (default: Testis).")
    parser.add_argument("--sample_id_col", default="SAMPID", 
                        help="Sample ID column (default: SAMPID).")
    parser.add_argument("--tissue_col", default="SMTSD", 
                        help="Tissue column (default: SMTSD).")
    parser.add_argument("--min_tpm_present", type=float, 
                        default=1.0, help="TPM threshold (default: 1.0).")
    parser.add_argument("--min_target_present_fraction",
                        type=float,
                        default=0.5,
                        help="Minimum fraction of target samples with TPM >= threshold (default: 0.5). Use <=0 to disable.",
                        )
    parser.add_argument("--top_k_tissues", type=int, 
                        default=5, help="Top tissues to store (default: 5).")
    parser.add_argument("--out_ranked_tsv", 
                        required=True, help="Output TSV path for ranked gene summary.")
    parser.add_argument("--out_tissue_medians_tsv",
                        required=True,
                        help="Output TSV path for gene × tissue median TPM matrix.",
                        )
    parser.add_argument("--log_path", default=None, help="Optional log file path.")
    parser.add_argument("--log_level",
                        default="INFO",
                        help="Logging level (DEBUG, INFO, WARNING). Default: INFO.",
                        )

    args = parser.parse_args()

    min_frac: Optional[float] = args.min_target_present_fraction
    if min_frac is not None and min_frac <= 0:
        min_frac = None

    return Config(
                gct_path=args.gct_path,
                sample_attributes_path=args.sample_attributes_path,
                target_tissue=args.target_tissue,
                sample_id_col=args.sample_id_col,
                tissue_col=args.tissue_col,
                min_tpm_present=args.min_tpm_present,
                min_target_present_fraction=min_frac,
                top_k_tissues=args.top_k_tissues,
                out_ranked_tsv=args.out_ranked_tsv,
                out_tissue_medians_tsv=args.out_tissue_medians_tsv,
                log_path=args.log_path,
                log_level=args.log_level,
            )


def timed(*, logger: logging.Logger, label: str):
    """
    Context manager for timing blocks.

    Parameters
    ----------
    logger
        Logger instance.
    label
        Label for the timed block.
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


def load_gct(*, gct_path: str, logger: logging.Logger) -> Tuple[pd.Series, pd.DataFrame]:
    """
    Load a GTEx GCT file into (description, expression) with logging.

    Parameters
    ----------
    gct_path
        Path to GCT file.
    logger
        Logger instance.

    Returns
    -------
    Tuple[pd.Series, pd.DataFrame]
        description: Series indexed by gene ID
        expression: DataFrame indexed by gene ID, columns are sample IDs
    """
    with timed(logger=logger, label="Loading GCT"):
        logger.info("Reading GCT: %s", gct_path)
        logger.info("GCT parsing uses header=2 (skips first two lines).")
        df = pd.read_csv(gct_path, sep="\t", header=2, dtype={0: str}, low_memory=False)

        logger.info("GCT loaded: %d rows x %d columns", df.shape[0], df.shape[1])
        logger.info("First 5 columns: %s", list(df.columns)[:5])

        if "Name" not in df.columns:
            raise ValueError("GCT parsing failed: expected 'Name' column at header=2.")

        df = df.set_index("Name", drop=True)
        logger.info("Set index to 'Name' (gene ID). Index size: %d", df.shape[0])

        if "Description" in df.columns:
            description = df["Description"].astype(str).fillna("")
            expr = df.drop(columns=["Description"])
            logger.info("Found 'Description' column; expression columns now: %d", expr.shape[1])
        else:
            description = pd.Series([""] * df.shape[0], index=df.index, name="Description")
            expr = df
            logger.info("No 'Description' column found; proceeding without gene symbols.")

        logger.info("Converting expression values to numeric (non-numeric -> 0).")
        expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        expr.columns = expr.columns.astype(str)

        logger.info("Expression matrix: %d genes x %d samples", expr.shape[0], expr.shape[1])

    return description, expr


def load_sample_attributes(
    *, sample_attributes_path: str, sample_id_col: str, tissue_col: str, logger: logging.Logger
) -> pd.DataFrame:
    """
    Load GTEx SampleAttributesDS with logging.

    Parameters
    ----------
    sample_attributes_path
        Path to SampleAttributesDS.txt.
    sample_id_col
        Sample ID column name.
    tissue_col
        Tissue label column name.
    logger
        Logger instance.

    Returns
    -------
    pd.DataFrame
        DataFrame with [sample_id_col, tissue_col].
    """
    with timed(logger=logger, label="Loading SampleAttributesDS"):
        logger.info("Reading sample attributes: %s", sample_attributes_path)
        df = pd.read_csv(sample_attributes_path, sep="\t", low_memory=False)
        logger.info("SampleAttributesDS loaded: %d rows x %d columns", df.shape[0], df.shape[1])

        for col in [sample_id_col, tissue_col]:
            if col not in df.columns:
                raise ValueError(
                    f"Column '{col}' not found. Example columns: {list(df.columns)[:30]}"
                )

        out = df[[sample_id_col, tissue_col]].copy()
        out[sample_id_col] = out[sample_id_col].astype(str)
        out[tissue_col] = out[tissue_col].astype(str)

        tissue_counts = out[tissue_col].value_counts(dropna=False).head(10)
        logger.info("Top 10 tissues by sample count:\n%s", tissue_counts.to_string())

    return out


def compute_tau(*, tissue_medians: pd.DataFrame, epsilon: float, logger: logging.Logger) -> pd.Series:
    """
    Compute Tau tissue-specificity score per gene from tissue medians, with logging.

    Tau is defined as:
        tau = sum(1 - x_i / x_max) / (n - 1)
    where x_i are tissue expression values (e.g. median TPM per tissue),
    x_max is max across tissues, n is number of tissues.

    Tau ranges:
      - 0: broadly expressed
      - 1: highly tissue-specific

    Parameters
    ----------
    tissue_medians
        DataFrame indexed by gene, columns are tissues.
    epsilon
        Small constant to avoid division by zero.
    logger
        Logger instance.

    Returns
    -------
    pd.Series
        Tau score per gene.
    """
    with timed(logger=logger, label="Computing tau scores"):
        logger.info("Tissue medians matrix: %d genes x %d tissues", tissue_medians.shape[0], tissue_medians.shape[1])
        x = tissue_medians.to_numpy(dtype=float)
        x_max = np.max(x, axis=1) + epsilon
        n_tissues = x.shape[1]
        denom = max(n_tissues - 1, 1)
        tau = np.sum(1.0 - (x / x_max[:, None]), axis=1) / float(denom)
        out = pd.Series(tau, index=tissue_medians.index, name="tau")
        logger.info("Tau computed. Summary: min=%.3g median=%.3g max=%.3g",
                    float(out.min()), float(out.median()), float(out.max()))
    return out


def format_top_tissues(*, row: pd.Series, top_k: int) -> str:
    """
    Format the top tissues for one gene as a compact string.

    Parameters
    ----------
    row
        Tissue median TPM values for a gene.
    top_k
        Number of top tissues to include.

    Returns
    -------
    str
        Semicolon-separated 'Tissue=TPM' entries.
    """
    top = row.sort_values(ascending=False).head(top_k)
    return ";".join([f"{t}={v:.3g}" for t, v in top.items()])


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run tissue median computation and target tissue specificity scoring.

    Parameters
    ----------
    cfg
        Configuration object.
    logger
        Logger instance.
    """
    logger.info("Starting GTEx specificity pipeline")
    logger.info("Config: target_tissue=%s min_tpm_present=%.3g min_target_present_fraction=%s top_k=%d",
                cfg.target_tissue,
                cfg.min_tpm_present,
                "None" if cfg.min_target_present_fraction is None else f"{cfg.min_target_present_fraction:.3g}",
                cfg.top_k_tissues)

    description, expr = load_gct(gct_path=cfg.gct_path, logger=logger)
    sample_attr = load_sample_attributes(
        sample_attributes_path=cfg.sample_attributes_path,
        sample_id_col=cfg.sample_id_col,
        tissue_col=cfg.tissue_col,
        logger=logger,
    )

    with timed(logger=logger, label="Intersecting samples"):
        expr_samples = pd.Index(expr.columns)
        attr_samples = pd.Index(sample_attr[cfg.sample_id_col])
        common_samples = expr_samples.intersection(attr_samples)

        logger.info("Expression samples: %d", len(expr_samples))
        logger.info("Attribute samples: %d", len(attr_samples))
        logger.info("Common samples: %d", len(common_samples))

        if common_samples.empty:
            raise ValueError("No overlapping sample IDs between GCT columns and SampleAttributesDS.")

        expr = expr.loc[:, common_samples]
        sample_attr = sample_attr.loc[sample_attr[cfg.sample_id_col].isin(common_samples), :]

        logger.info("Filtered expression matrix: %d genes x %d samples", expr.shape[0], expr.shape[1])
        logger.info("Filtered sample attributes: %d rows", sample_attr.shape[0])

    with timed(logger=logger, label="Mapping samples to tissues"):
        sample_to_tissue = sample_attr.set_index(cfg.sample_id_col)[cfg.tissue_col]
        tissues = sample_to_tissue.loc[common_samples].astype(str)

        logger.info("Unique tissues in filtered set: %d", tissues.nunique(dropna=False))
        logger.info("Top 10 tissues in filtered set:\n%s",
                    tissues.value_counts(dropna=False).head(10).to_string())

        target_n = int((tissues.values == cfg.target_tissue).sum())
        logger.info("Target tissue '%s' sample count: %d", cfg.target_tissue, target_n)

    with timed(logger=logger, label="Computing gene × tissue medians"):
        logger.info("Transposing expression matrix for groupby median by tissue.")
        tissue_medians = (expr.T.assign(_tissue=tissues.values)
                         .groupby("_tissue", dropna=False)
                         .median(numeric_only=True)
                         .T
                       )

        logger.info("Computed tissue medians: %d genes x %d tissues", tissue_medians.shape[0], tissue_medians.shape[1])

        if cfg.target_tissue not in tissue_medians.columns:
            examples = sorted(tissue_medians.columns.astype(str).tolist())[:25]
            raise ValueError(
                f"Target tissue '{cfg.target_tissue}' not found. Example tissues: {examples}"
            )

        logger.info("Writing tissue medians TSV: %s", cfg.out_tissue_medians_tsv)
        tissue_medians.to_csv(cfg.out_tissue_medians_tsv, sep="\t", index=True)

    tau = compute_tau(tissue_medians=tissue_medians, epsilon=cfg.epsilon, logger=logger)

    with timed(logger=logger, label="Scoring target tissue specificity"):
        target_median = tissue_medians[cfg.target_tissue]
        non_target = tissue_medians.drop(columns=[cfg.target_tissue])

        if non_target.shape[1] > 0:
            max_non_target_tissue = non_target.idxmax(axis=1)
            max_non_target_median = non_target.max(axis=1)
        else:
            max_non_target_tissue = pd.Series([""] * tissue_medians.shape[0], index=tissue_medians.index)
            max_non_target_median = pd.Series([0.0] * tissue_medians.shape[0], index=tissue_medians.index)

        log2_fc = np.log2((target_median + cfg.epsilon) / (max_non_target_median + cfg.epsilon))
        log2_fc.name = "log2_fc_target_vs_max_non_target"

        target_samples = common_samples[tissues.values == cfg.target_tissue]
        if len(target_samples) == 0:
            raise ValueError(f"No samples labelled as '{cfg.target_tissue}' in SampleAttributesDS.")

        logger.info("Computing target presence fraction using TPM >= %.3g", cfg.min_tpm_present)
        target_present_fraction = (
            (expr.loc[:, target_samples] >= cfg.min_tpm_present).sum(axis=1) / float(len(target_samples))
        )
        target_present_fraction.name = "target_present_fraction"

        is_target_max = (tissue_medians.idxmax(axis=1) == cfg.target_tissue).astype(int)
        is_target_max.name = "is_target_max_tissue"

        logger.info("Formatting top tissues strings (top_k=%d). This can take a little while.", cfg.top_k_tissues)
        top_tissues = tissue_medians.apply(lambda r: format_top_tissues(row=r, top_k=cfg.top_k_tissues), axis=1)
        top_tissues.name = "top_tissues_by_median_tpm"

        ranked = pd.DataFrame(
            {
                "Description": description.reindex(tissue_medians.index).fillna(""),
                "target_median_tpm": target_median,
                "max_non_target_tissue": max_non_target_tissue,
                "max_non_target_median_tpm": max_non_target_median,
                "log2_fc_target_vs_max_non_target": log2_fc,
                "tau": tau,
                "target_present_fraction": target_present_fraction,
                "is_target_max_tissue": is_target_max,
                "top_tissues_by_median_tpm": top_tissues,
            }
        )

        logger.info("Ranked table initial size: %d genes", ranked.shape[0])

        if cfg.min_target_present_fraction is not None:
            before = ranked.shape[0]
            ranked = ranked.loc[ranked["target_present_fraction"] >= cfg.min_target_present_fraction, :]
            after = ranked.shape[0]
            logger.info(
                "Applied target presence filter (>= %.3g): %d -> %d genes",
                cfg.min_target_present_fraction,
                before,
                after,
            )

        logger.info("Sorting ranked table by tau, enrichment, and target median TPM.")
        ranked = ranked.sort_values(
            by=["tau", "log2_fc_target_vs_max_non_target", "target_median_tpm"],
            ascending=[False, False, False],
        )

        logger.info("Writing ranked TSV: %s", cfg.out_ranked_tsv)
        ranked.to_csv(cfg.out_ranked_tsv, sep="\t", index=True)

        logger.info("Top 5 ranked genes (index + key metrics):\n%s",
                    ranked[["target_median_tpm", "log2_fc_target_vs_max_non_target", "tau",
                            "max_non_target_tissue", "max_non_target_median_tpm"]].head(5).to_string())

    logger.info("Pipeline finished successfully.")


def main() -> None:
    """
    Command line entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    try:
        run(cfg=cfg, logger=logger)
    except Exception as exc:
        logger.exception("Fatal error: %s", str(exc))
        raise


if __name__ == "__main__":
    main()
