#!/usr/bin/env python3
"""
Annotate GTEx testis-enriched genes with ClinVar disease links.

This script takes a ranked gene table produced from GTEx (Script 1) and merges:
  - ClinVar gene-to-condition associations (gene_condition_source_id)
  - Optionally, ClinVar variant summaries (variant_summary.txt.gz)

It outputs a TSV with added columns, including a compact list of conditions
and per-gene variant counts by clinical significance.

Notes
-----
- Your GTEx ranked file uses Ensembl gene IDs as the index and may include a
  gene symbol in the 'Description' column. ClinVar files are primarily keyed
  by gene symbol, so this script uses 'Description' as the join key by default.
- If 'Description' is not a gene symbol in your GCT, you can provide an
  alternative symbol column (or supply an external mapping later).

All outputs are tab-separated (TSV).
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd


@dataclass(frozen=True)
class Config:
    """
    Configuration for ClinVar annotation.

    Attributes
    ----------
    ranked_genes_tsv
        Ranked GTEx output TSV from Script 1.
    gene_symbol_col
        Column in ranked_genes_tsv containing gene symbols (default: 'Description').
    clinvar_gene_condition_path
        Path to ClinVar gene_condition_source_id file.
    clinvar_variant_summary_gz
        Optional path to ClinVar variant_summary.txt.gz.
    output_tsv
        Output TSV path.
    max_conditions_per_gene
        Maximum number of condition names to store per gene in the output.
    include_variant_counts
        Whether to parse variant_summary.txt.gz and include counts.
    restrict_to_germline
        If True, restrict variant_summary counts to Germline origin where available.
    log_path
        Optional log file path.
    log_level
        Logging level.
    """

    ranked_genes_tsv: str
    gene_symbol_col: str
    clinvar_gene_condition_path: str
    clinvar_variant_summary_gz: Optional[str]
    output_tsv: str
    max_conditions_per_gene: int = 10
    include_variant_counts: bool = True
    restrict_to_germline: bool = False
    log_path: Optional[str] = None
    log_level: str = "INFO"


def setup_logging(*, log_level: str, log_path: Optional[str]) -> logging.Logger:
    """
    Configure logging to stdout and optionally a file.

    Parameters
    ----------
    log_level
        Logging level string.
    log_path
        Optional file path for logs.

    Returns
    -------
    logging.Logger
        Logger instance.
    """
    logger = logging.getLogger("clinvar_annotator")
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
    Context manager for timing a block.

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
    Parse command line arguments.

    Returns
    -------
    Config
        Parsed configuration.
    """
    parser = argparse.ArgumentParser(
        description="Annotate GTEx testis-enriched genes with ClinVar gene-condition and variant summaries."
    )
    parser.add_argument(
        "--ranked_genes_tsv",
        required=True,
        help="Ranked GTEx output TSV from Script 1.",
    )
    parser.add_argument(
        "--gene_symbol_col",
        default="Description",
        help="Column in ranked_genes_tsv containing gene symbols (default: Description).",
    )
    parser.add_argument(
        "--clinvar_gene_condition_path",
        required=True,
        help="Path to ClinVar gene_condition_source_id (uncompressed).",
    )
    parser.add_argument(
        "--clinvar_variant_summary_gz",
        default=None,
        help="Optional path to ClinVar tab_delimited/variant_summary.txt.gz.",
    )
    parser.add_argument(
        "--include_variant_counts",
        action="store_true",
        help="If set, parse variant_summary and include per-gene variant counts.",
    )
    parser.add_argument(
        "--restrict_to_germline",
        action="store_true",
        help="If set, restrict variant_summary counts to Germline origin (when available).",
    )
    parser.add_argument(
        "--max_conditions_per_gene",
        type=int,
        default=10,
        help="Maximum number of condition names stored per gene (default: 10).",
    )
    parser.add_argument(
        "--output_tsv",
        required=True,
        help="Output TSV path.",
    )
    parser.add_argument("--log_path", default=None, help="Optional log file path.")
    parser.add_argument("--log_level", default="INFO", help="Log level (INFO, DEBUG).")

    args = parser.parse_args()

    return Config(
        ranked_genes_tsv=args.ranked_genes_tsv,
        gene_symbol_col=args.gene_symbol_col,
        clinvar_gene_condition_path=args.clinvar_gene_condition_path,
        clinvar_variant_summary_gz=args.clinvar_variant_summary_gz,
        output_tsv=args.output_tsv,
        max_conditions_per_gene=args.max_conditions_per_gene,
        include_variant_counts=bool(args.include_variant_counts),
        restrict_to_germline=bool(args.restrict_to_germline),
        log_path=args.log_path,
        log_level=args.log_level,
    )


def normalise_gene_symbol(*, symbol: str) -> str:
    """
    Normalise gene symbols for joining.

    Parameters
    ----------
    symbol
        Gene symbol string.

    Returns
    -------
    str
        Normalised symbol.
    """
    if symbol is None:
        return ""
    sym = str(symbol).strip()
    if sym in {"", "nan", "None"}:
        return ""
    return sym.upper()


def load_ranked_genes(*, path: str, gene_symbol_col: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load ranked GTEx gene table.

    Parameters
    ----------
    path
        Ranked TSV.
    gene_symbol_col
        Column containing gene symbol.
    logger
        Logger.

    Returns
    -------
    pd.DataFrame
        Ranked genes with an added 'gene_symbol_norm' column.
    """
    with timed(logger=logger, label="Loading ranked GTEx genes"):
        df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
        logger.info("Loaded ranked genes: %d rows x %d columns", df.shape[0], df.shape[1])

        if gene_symbol_col not in df.columns:
            raise ValueError(
                f"Column '{gene_symbol_col}' not found in ranked file. "
                f"Available columns include: {list(df.columns)[:30]}"
            )

        df = df.copy()
        df["gene_symbol_norm"] = df[gene_symbol_col].map(lambda x: normalise_gene_symbol(symbol=x))
        n_missing = int((df["gene_symbol_norm"] == "").sum())
        logger.info("Normalised gene symbols from column '%s'. Missing/blank: %d", gene_symbol_col, n_missing)

    return df



def parse_gene_condition_source_id(
    *, path: str, logger: logging.Logger, max_conditions_per_gene: int
) -> pd.DataFrame:
    """
    Parse ClinVar gene_condition_source_id into per-gene condition summaries.

    This version supports ClinVar releases where the gene field is stored in
    'AssociatedGenes' (and optionally 'RelatedGenes').

    Parameters
    ----------
    path
        Path to gene_condition_source_id.
    logger
        Logger.
    max_conditions_per_gene
        Maximum number of condition names to keep per gene.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by normalised gene symbol with columns:
          - clinvar_condition_count
          - clinvar_conditions (semicolon-separated)
    """
    with timed(logger=logger, label="Parsing ClinVar gene_condition_source_id"):
        df = pd.read_csv(path, sep="\t", low_memory=False)
        logger.info("Loaded gene_condition_source_id: %d rows x %d columns", df.shape[0], df.shape[1])
        logger.info("Columns: %s", list(df.columns))

        # ClinVar commonly uses 'AssociatedGenes' for gene symbols in this file.
        possible_gene_cols = [
            "GeneSymbol",
            "gene_symbol",
            "Symbol",
            "symbol",
            "Gene",
            "AssociatedGenes",
            "RelatedGenes",
        ]
        gene_col = next((c for c in possible_gene_cols if c in df.columns), None)
        if gene_col is None:
            raise ValueError(
                "Could not find a gene symbol column in gene_condition_source_id. "
                "Expected one of: GeneSymbol, AssociatedGenes, RelatedGenes."
            )

        # ClinVar condition name column
        possible_condition_cols = ["ConditionName", "condition_name", "DiseaseName", "disease_name", "Phenotype"]
        condition_col = next((c for c in possible_condition_cols if c in df.columns), None)
        if condition_col is None:
            candidates = [c for c in df.columns if "condition" in c.lower() or "disease" in c.lower()]
            raise ValueError(
                "Could not find a condition name column in gene_condition_source_id. "
                f"Candidates containing 'condition'/'disease': {candidates}"
            )

        logger.info("Using gene column: %s", gene_col)
        logger.info("Using condition column: %s", condition_col)

        df = df[[gene_col, condition_col]].copy()
        df[gene_col] = df[gene_col].astype(str).fillna("")
        df[condition_col] = df[condition_col].astype(str).fillna("")

        def _split_genes(g: str) -> List[str]:
            # Handle multiple genes in one field
            # Seen delimiters in ClinVar-like files: comma, semicolon, pipe
            raw = str(g).strip()
            if not raw:
                return []
            for delim in ["|", ";", ","]:
                raw = raw.replace(delim, ",")
            parts = [p.strip() for p in raw.split(",") if p.strip()]
            return parts

        records: List[Tuple[str, str]] = []
        for _, row in df.iterrows():
            cond = str(row[condition_col]).strip()
            if not cond:
                continue
            genes = _split_genes(str(row[gene_col]))
            for g in genes:
                g_norm = normalise_gene_symbol(symbol=g)
                if g_norm:
                    records.append((g_norm, cond))

        if not records:
            raise ValueError("No (gene, condition) pairs could be parsed from gene_condition_source_id.")

        parsed = pd.DataFrame(records, columns=["gene_symbol_norm", "condition"])
        logger.info("Parsed gene-condition pairs: %d", parsed.shape[0])

        def _collapse_conditions(conds: pd.Series) -> str:
            uniq: List[str] = []
            seen: Set[str] = set()
            for c in conds.astype(str):
                c2 = c.strip()
                if not c2 or c2 in seen:
                    continue
                uniq.append(c2)
                seen.add(c2)
                if len(uniq) >= max_conditions_per_gene:
                    break
            return ";".join(uniq)

        grouped = parsed.groupby("gene_symbol_norm")["condition"].agg(
            clinvar_condition_count="nunique",
            clinvar_conditions=_collapse_conditions,
        )

        logger.info("Generated per-gene condition summaries: %d genes", grouped.shape[0])
        return grouped



def parse_variant_summary_counts(
    *,
    path_gz: str,
    logger: logging.Logger,
    restrict_to_germline: bool,
) -> pd.DataFrame:
    """
    Parse ClinVar variant_summary.txt.gz and compute per-gene counts.

    Parameters
    ----------
    path_gz
        Path to variant_summary.txt.gz.
    logger
        Logger.
    restrict_to_germline
        Whether to restrict to Germline origin where available.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by normalised gene symbol with columns like:
          - clinvar_variants_total
          - clinvar_pathogenic
          - clinvar_likely_pathogenic
          - clinvar_benign
          - clinvar_likely_benign
          - clinvar_uncertain_significance
    """
    with timed(logger=logger, label="Parsing ClinVar variant_summary counts"):
        logger.info("Reading gzipped variant summary: %s", path_gz)

        # Read header only to detect columns.
        with gzip.open(path_gz, "rt") as handle:
            header_line = handle.readline().rstrip("\n")
        cols = header_line.split("\t")
        logger.info("variant_summary header columns: %d", len(cols))
        logger.info("First 25 columns: %s", cols[:25])

        possible_gene_cols = ["GeneSymbol", "GeneSymbolList", "GeneSymbols", "Gene"]
        gene_col = next((c for c in possible_gene_cols if c in cols), None)
        if gene_col is None:
            raise ValueError(
                "Could not find a gene symbol column in variant_summary. "
                "Expected one of: GeneSymbol, GeneSymbolList."
            )

        clin_sig_col = "ClinicalSignificance"
        if clin_sig_col not in cols:
            raise ValueError("ClinicalSignificance column not found in variant_summary.")

        origin_col = "OriginSimple"
        has_origin = origin_col in cols

        usecols = [gene_col, clin_sig_col]
        if restrict_to_germline and has_origin:
            usecols.append(origin_col)

        df = pd.read_csv(path_gz, sep="\t", compression="gzip", usecols=usecols, low_memory=False)
        logger.info("variant_summary loaded with usecols=%s: %d rows", usecols, df.shape[0])

        df["gene_symbol_norm"] = df[gene_col].map(lambda x: normalise_gene_symbol(symbol=x))
        df = df.loc[df["gene_symbol_norm"] != "", :].copy()

        if restrict_to_germline:
            if not has_origin:
                logger.warning(
                    "OriginSimple not present in variant_summary; cannot restrict to germline. Proceeding without filter."
                )
            else:
                before = df.shape[0]
                df = df.loc[df[origin_col].astype(str).str.contains("germline", case=False, na=False), :]
                logger.info("Restricted to germline: %d -> %d rows", before, df.shape[0])

        def _map_sig(sig: str) -> str:
            s = str(sig).strip().lower()
            if "pathogenic" in s and "likely" not in s:
                return "pathogenic"
            if "likely pathogenic" in s:
                return "likely_pathogenic"
            if "benign" in s and "likely" not in s:
                return "benign"
            if "likely benign" in s:
                return "likely_benign"
            if "uncertain significance" in s or "vus" in s:
                return "uncertain_significance"
            return "other"

        df["sig_bucket"] = df[clin_sig_col].map(_map_sig)

        counts_total = df.groupby("gene_symbol_norm").size().rename("clinvar_variants_total")
        counts_bucket = (
            df.groupby(["gene_symbol_norm", "sig_bucket"])
            .size()
            .unstack(fill_value=0)
        )

        # Ensure consistent columns
        for col in [
            "pathogenic",
            "likely_pathogenic",
            "benign",
            "likely_benign",
            "uncertain_significance",
            "other",
        ]:
            if col not in counts_bucket.columns:
                counts_bucket[col] = 0

        counts_bucket = counts_bucket.rename(
            columns={
                "pathogenic": "clinvar_pathogenic",
                "likely_pathogenic": "clinvar_likely_pathogenic",
                "benign": "clinvar_benign",
                "likely_benign": "clinvar_likely_benign",
                "uncertain_significance": "clinvar_uncertain_significance",
                "other": "clinvar_other_significance",
            }
        )

        out = pd.concat([counts_total, counts_bucket], axis=1).fillna(0).astype(int)
        logger.info("Generated per-gene variant counts: %d genes", out.shape[0])

    return out


def run(*, cfg: Config, logger: logging.Logger) -> None:
    """
    Run ClinVar annotation and write output.

    Parameters
    ----------
    cfg
        Configuration.
    logger
        Logger.
    """
    ranked = load_ranked_genes(path=cfg.ranked_genes_tsv, gene_symbol_col=cfg.gene_symbol_col, logger=logger)

    gene_conditions = parse_gene_condition_source_id(
        path=cfg.clinvar_gene_condition_path,
        logger=logger,
        max_conditions_per_gene=cfg.max_conditions_per_gene,
    )

    with timed(logger=logger, label="Merging gene-condition annotations"):
        out = ranked.join(gene_conditions, on="gene_symbol_norm")
        out["clinvar_condition_count"] = out["clinvar_condition_count"].fillna(0).astype(int)
        out["clinvar_conditions"] = out["clinvar_conditions"].fillna("")
        logger.info("After merge: %d rows x %d columns", out.shape[0], out.shape[1])

        n_annot = int((out["clinvar_condition_count"] > 0).sum())
        logger.info(
            "Genes with >=1 ClinVar condition: %d / %d (%.1f%%)",
            n_annot,
            out.shape[0],
            100.0 * n_annot / out.shape[0],
        )


    if cfg.include_variant_counts:
        if not cfg.clinvar_variant_summary_gz:
            raise ValueError("--include_variant_counts was set but --clinvar_variant_summary_gz is missing.")
        variant_counts = parse_variant_summary_counts(
            path_gz=cfg.clinvar_variant_summary_gz,
            logger=logger,
            restrict_to_germline=cfg.restrict_to_germline,
        )
        with timed(logger=logger, label="Merging variant count annotations"):
            out = out.join(variant_counts, on="gene_symbol_norm")
            for col in [
                "clinvar_variants_total",
                "clinvar_pathogenic",
                "clinvar_likely_pathogenic",
                "clinvar_benign",
                "clinvar_likely_benign",
                "clinvar_uncertain_significance",
                "clinvar_other_significance",
            ]:
                if col in out.columns:
                    out[col] = out[col].fillna(0).astype(int)
            logger.info("After merge: %d rows x %d columns", out.shape[0], out.shape[1])
    else:
        logger.info("Variant counts disabled; skipping variant_summary parsing.")

    with timed(logger=logger, label="Writing output"):
        logger.info("Writing annotated TSV: %s", cfg.output_tsv)
        out.to_csv(cfg.output_tsv, sep="\t", index=True)

        preview_cols = [
            cfg.gene_symbol_col,
            "target_median_tpm",
            "tau",
            "clinvar_condition_count",
        ]
        preview_cols = [c for c in preview_cols if c in out.columns]
        logger.info("Preview (top 10 rows):\n%s", out[preview_cols].head(10).to_string())

    logger.info("Done.")


def main() -> None:
    """
    Entry point.
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
