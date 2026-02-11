#!/usr/bin/env python3
"""
Add HPO (Human Phenotype Ontology) gene-to-phenotype annotations to a TSV.

This script joins gene-level phenotype information from HPO's
'genes_to_phenotype.txt' onto an existing table keyed by HGNC-style gene
symbols (for example, 'hgnc_symbol_norm').

It produces per-gene summary columns:
- hpo_term_count: number of unique HPO terms associated with the gene
- hpo_terms: semicolon-separated unique HPO term names
- hpo_ids: semicolon-separated unique HPO IDs
- hpo_disease_ids: semicolon-separated unique disease IDs (OMIM/ORPHA etc.)
- hpo_reproductive_term_count: number of HPO terms matching reproductive keywords
- hpo_reproductive_terms: matching reproductive phenotype names
- hpo_reproductive_ids: matching reproductive HPO IDs

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
    Configuration for HPO annotation.

    Attributes
    ----------
    input_tsv
        Path to the input TSV to annotate.
    output_tsv
        Path to the output TSV with added HPO columns.
    gene_symbol_col
        Column in input_tsv containing gene symbols to join on.
    hpo_genes_to_phenotype_tsv
        Path to HPO genes_to_phenotype.txt file.
    hpo_gene_symbol_col
        Column name in HPO file containing gene symbols.
    hpo_id_col
        Column name in HPO file containing HPO IDs.
    hpo_name_col
        Column name in HPO file containing HPO term names.
    hpo_disease_id_col
        Column name in HPO file containing disease IDs.
    max_terms_per_gene
        Maximum number of terms to keep per gene (highest-level truncation).
    reproductive_regex
        Regular expression used to select reproductive-related phenotype terms.
    log_level
        Logging level.
    log_path
        Optional log file path.
    """

    input_tsv: str
    output_tsv: str
    gene_symbol_col: str
    hpo_genes_to_phenotype_tsv: str
    hpo_gene_symbol_col: str
    hpo_id_col: str
    hpo_name_col: str
    hpo_disease_id_col: str
    max_terms_per_gene: int
    reproductive_regex: str
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
        Configured logger.
    """
    logger = logging.getLogger("add_hpo_annotation")
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
    p = argparse.ArgumentParser(description="Add HPO gene-to-phenotype annotations to a TSV.")
    p.add_argument("--input_tsv", required=True, help="Input TSV to annotate.")
    p.add_argument("--output_tsv", required=True, help="Output TSV with HPO columns added.")
    p.add_argument(
        "--gene_symbol_col",
        default="hgnc_symbol_norm",
        help="Gene symbol column in input TSV (default: hgnc_symbol_norm).",
    )
    p.add_argument(
        "--hpo_genes_to_phenotype_tsv",
        required=True,
        help="Path to HPO genes_to_phenotype.txt (tab-separated).",
    )
    p.add_argument(
        "--hpo_gene_symbol_col",
        default="gene_symbol",
        help="Gene symbol column in HPO file (default: gene_symbol).",
    )
    p.add_argument(
        "--hpo_id_col",
        default="hpo_id",
        help="HPO ID column in HPO file (default: hpo_id).",
    )
    p.add_argument(
        "--hpo_name_col",
        default="hpo_name",
        help="HPO name column in HPO file (default: hpo_name).",
    )
    p.add_argument(
        "--hpo_disease_id_col",
        default="disease_id",
        help="Disease ID column in HPO file (default: disease_id).",
    )
    p.add_argument(
        "--max_terms_per_gene",
        type=int,
        default=50,
        help="Maximum unique HPO terms to keep per gene (default: 50).",
    )
    p.add_argument(
        "--reproductive_regex",
        default=r"(sperm|spermat|testis|testicular|infertil|azoosperm|oligozoosperm|"
                r"asthenozoosperm|teratozoosperm|hypogonad|gonad|cryptorchid|"
                r"androgen|seminifer|epididym|prostat|ejacul|erect|libido|"
                r"ovary|ovarian|uterus|uterine|endometr|oocyte|follic|"
                r"amenorr|menorr|menarch|menopaus|pregnan|miscarriage|"
                r"placent|fetal|foetal)",
        help="Regex to flag reproductive phenotype terms (case-insensitive).",
    )
    p.add_argument("--log_level", default="INFO", help="Log level.")
    p.add_argument("--log_path", default=None, help="Optional log file path.")
    a = p.parse_args()

    return Config(
        input_tsv=a.input_tsv,
        output_tsv=a.output_tsv,
        gene_symbol_col=a.gene_symbol_col,
        hpo_genes_to_phenotype_tsv=a.hpo_genes_to_phenotype_tsv,
        hpo_gene_symbol_col=a.hpo_gene_symbol_col,
        hpo_id_col=a.hpo_id_col,
        hpo_name_col=a.hpo_name_col,
        hpo_disease_id_col=a.hpo_disease_id_col,
        max_terms_per_gene=a.max_terms_per_gene,
        reproductive_regex=a.reproductive_regex,
        log_level=a.log_level,
        log_path=a.log_path,
    )


def normalise_symbol(*, value: object) -> str:
    """
    Normalise a gene symbol to uppercase without surrounding whitespace.

    Parameters
    ----------
    value
        Gene symbol.

    Returns
    -------
    str
        Normalised symbol, or empty string if missing.
    """
    if value is None:
        return ""
    s = str(value).strip().upper()
    if s in {"", "NAN"}:
        return ""
    return s


def dedupe_preserve_order(*, items: List[str], max_items: int) -> List[str]:
    """
    Deduplicate strings while preserving order, truncated to max_items.

    Parameters
    ----------
    items
        List of items.
    max_items
        Max number of items to return.

    Returns
    -------
    list
        Deduplicated list (order-preserving).
    """
    seen = set()
    out: List[str] = []
    for it in items:
        if it not in seen and it != "":
            seen.add(it)
            out.append(it)
        if len(out) >= max_items:
            break
    return out


def load_hpo_gene_to_phenotype(
    *,
    path: str,
    gene_symbol_col: str,
    hpo_id_col: str,
    hpo_name_col: str,
    disease_id_col: str,
    max_terms_per_gene: int,
    reproductive_pattern: re.Pattern,
    logger: logging.Logger,
) -> pd.DataFrame:
    """
    Load and summarise HPO genes_to_phenotype into per-gene annotation rows.

    Parameters
    ----------
    path
        Path to genes_to_phenotype.txt.
    gene_symbol_col
        Column name for gene symbol.
    hpo_id_col
        Column name for HPO ID.
    hpo_name_col
        Column name for HPO term name.
    disease_id_col
        Column name for disease ID.
    max_terms_per_gene
        Maximum number of terms retained per gene.
    reproductive_pattern
        Compiled regex pattern for reproductive terms.
    logger
        Logger.

    Returns
    -------
    pd.DataFrame
        Per-gene annotation table indexed by normalised gene symbol.
    """
    t0 = time.time()
    logger.info("Loading HPO genes_to_phenotype: %s", path)
    hpo = pd.read_csv(path, sep="\t", low_memory=False)
    logger.info("HPO loaded: %d rows x %d columns", hpo.shape[0], hpo.shape[1])

    for col in [gene_symbol_col, hpo_id_col, hpo_name_col, disease_id_col]:
        if col not in hpo.columns:
            raise ValueError(f"Missing required column in HPO file: {col}")

    hpo = hpo.copy()
    hpo["gene_symbol_norm"] = hpo[gene_symbol_col].map(lambda x: normalise_symbol(value=x))
    hpo = hpo.loc[hpo["gene_symbol_norm"] != "", :].copy()

    hpo["hpo_id_str"] = hpo[hpo_id_col].astype(str).str.strip()
    hpo["hpo_name_str"] = hpo[hpo_name_col].astype(str).str.strip()
    hpo["disease_id_str"] = hpo[disease_id_col].astype(str).str.strip()

    # Reproductive flag
    hpo["is_reproductive"] = hpo["hpo_name_str"].str.contains(reproductive_pattern)

    def summarise_group(df: pd.DataFrame) -> pd.Series:
        """
        Summarise one gene's HPO rows into aggregated strings.

        Parameters
        ----------
        df
            Group DataFrame.

        Returns
        -------
        pd.Series
            Aggregated values for the gene.
        """
        ids = dedupe_preserve_order(items=df["hpo_id_str"].tolist(), max_items=max_terms_per_gene)
        names = dedupe_preserve_order(items=df["hpo_name_str"].tolist(), max_items=max_terms_per_gene)
        disease_ids = dedupe_preserve_order(items=df["disease_id_str"].tolist(), max_items=max_terms_per_gene)

        rep = df.loc[df["is_reproductive"], :]
        rep_ids = dedupe_preserve_order(items=rep["hpo_id_str"].tolist(), max_items=max_terms_per_gene)
        rep_names = dedupe_preserve_order(items=rep["hpo_name_str"].tolist(), max_items=max_terms_per_gene)

        return pd.Series(
            {
                "hpo_term_count": len(ids),
                "hpo_ids": ";".join(ids),
                "hpo_terms": ";".join(names),
                "hpo_disease_ids": ";".join(disease_ids),
                "hpo_reproductive_term_count": len(rep_ids),
                "hpo_reproductive_ids": ";".join(rep_ids),
                "hpo_reproductive_terms": ";".join(rep_names),
            }
        )

    logger.info("Summarising HPO annotations per gene (max_terms_per_gene=%d)", max_terms_per_gene)
    ann = hpo.groupby("gene_symbol_norm", sort=False).apply(summarise_group)

    ann.index.name = "gene_symbol_norm"
    ann = ann.copy()

    logger.info(
        "HPO per-gene annotation rows: %d (loaded+summarised in %.1f s)",
        ann.shape[0],
        time.time() - t0,
    )
    return ann


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run the HPO annotation join and write output.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    logger.info("Loading input TSV: %s", cfg.input_tsv)
    df = pd.read_csv(cfg.input_tsv, sep="\t", index_col=0, low_memory=False)
    logger.info("Input loaded: %d rows x %d columns", df.shape[0], df.shape[1])

    if cfg.gene_symbol_col not in df.columns:
        raise ValueError(
            f"Could not find --gene_symbol_col '{cfg.gene_symbol_col}' in input TSV. "
            f"Available columns: {list(df.columns)[:30]}"
        )

    df = df.copy()
    df["gene_symbol_norm"] = df[cfg.gene_symbol_col].map(lambda x: normalise_symbol(value=x))

    reproductive_pattern = re.compile(cfg.reproductive_regex, flags=re.IGNORECASE)
    hpo_ann = load_hpo_gene_to_phenotype(
        path=cfg.hpo_genes_to_phenotype_tsv,
        gene_symbol_col=cfg.hpo_gene_symbol_col,
        hpo_id_col=cfg.hpo_id_col,
        hpo_name_col=cfg.hpo_name_col,
        disease_id_col=cfg.hpo_disease_id_col,
        max_terms_per_gene=cfg.max_terms_per_gene,
        reproductive_pattern=reproductive_pattern,
        logger=logger,
    )

    logger.info("Merging HPO annotations onto input table")
    out = df.join(hpo_ann, on="gene_symbol_norm")

    # Fill missing annotations sensibly
    for col in [
        "hpo_term_count",
        "hpo_reproductive_term_count",
    ]:
        if col in out.columns:
            out[col] = out[col].fillna(0).astype(int)

    for col in [
        "hpo_ids",
        "hpo_terms",
        "hpo_disease_ids",
        "hpo_reproductive_ids",
        "hpo_reproductive_terms",
    ]:
        if col in out.columns:
            out[col] = out[col].fillna("").astype(str)

    n_any = int((out["hpo_term_count"] > 0).sum()) if "hpo_term_count" in out.columns else 0
    n_rep = int((out["hpo_reproductive_term_count"] > 0).sum()) if "hpo_reproductive_term_count" in out.columns else 0

    logger.info("Genes with >=1 HPO term: %d / %d", n_any, out.shape[0])
    logger.info("Genes with >=1 reproductive HPO term: %d / %d", n_rep, out.shape[0])

    logger.info("Writing output TSV: %s", cfg.output_tsv)
    out.to_csv(cfg.output_tsv, sep="\t", index=True)


def main() -> None:
    """
    Entry point.
    """
    cfg = parse_args()
    logger = setup_logging(log_level=cfg.log_level, log_path=cfg.log_path)
    run(cfg=cfg, logger=logger)


if __name__ == "__main__":
    main()
