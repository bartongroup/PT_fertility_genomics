#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summerise_number_in_RNA.py

Summarise overlaps between evidence layers (RNA and proteomics) using TSV inputs.

This script is designed for the sperm/testis prioritisation workflow and aims to
produce publication-friendly summaries and reproducible gene lists.

Inputs
------
1) GTEx testis specificity TSV (required)
   Expected columns (some may be optional depending on your upstream tables):
   - ensembl_gene_id or gene_key
   - tau
   - log2_fc_target_vs_max_non_target
   - is_target_max_tissue

2) Sperm evidence TSV (required)
   Expected columns (some may be optional depending on your upstream tables):
   - ensembl_gene_id or gene_key
   - sperm_tpm_mean and/or sperm_tpm_median and/or sperm_present_any
   - prot_present_any (for internal proteomics support if present)

3) Optional public proteomics TSV (optional)
   Must contain ensembl_gene_id or gene_key. Optionally may include a presence
   column, but by default any row with a gene id is treated as present.

Outputs (written to --out_dir, default: --base_dir)
---------------------------------------------------
- layer_counts.tsv
- pairwise_overlaps.tsv (includes percentage overlaps)
- combination_overlaps.tsv (all combinations, includes percentage overlaps)
- gene_lists/
  - <LAYER>.tsv (all genes passing criteria for that layer)
  - <LAYER>__ONLY.tsv (genes unique to that layer)
  - <LAYER_A>__AND__<LAYER_B>.tsv (pairwise intersections)
  - <LAYER_A>__AND__<LAYER_B>__AND__<LAYER_C>.tsv (triple intersections)
  - ... and the full intersection if applicable
  - sperm_only_genes.tsv (alias for RNA_sperm__ONLY.tsv)

Notes
-----
- All outputs are tab-separated (TSV).
- Gene keys are normalised to strip Ensembl version suffixes (e.g. ENSG... .12).

Example
-------
python summerise_number_in_RNA.py \
  --base_dir /home/pthorpe001/data/2026_sperm_Gates/results/testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --gtex_tsv /home/pthorpe001/data/2026_sperm_Gates/results/testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1/02_hgnc_mapping/gtex_v11_testis_specificity_ranked_with_hgnc.tsv \
  --sperm_tsv /home/pthorpe001/data/2026_sperm_Gates/results/testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1/09_high_confidence_final/high_confidence_testis_sperm_genes_function_hpo_prot.tsv \
  --tau_threshold 0.95 \
  --log2fc_threshold 0 \
  --require_testis_dominance \
  --sperm_tpm_threshold 0.1 \
  --verbose
"""

from __future__ import annotations

import argparse
import itertools
import logging
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd


LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class LayerSpec:
    """Definition of a layer to summarise and export."""

    name: str
    genes: Set[str]


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise overlaps between RNA and proteomics evidence layers."
    )

    parser.add_argument(
        "--base_dir",
        required=True,
        help="Base directory for the run (used for default output location).",
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory (default: --base_dir).",
    )

    parser.add_argument(
        "--gtex_tsv",
        required=True,
        help="GTEx testis specificity TSV.",
    )
    parser.add_argument(
        "--sperm_tsv",
        required=True,
        help="Sperm evidence TSV (RNA + may include internal proteomics columns).",
    )
    parser.add_argument(
        "--public_proteomics_tsv",
        default=None,
        help=(
            "Optional public proteomics TSV. If provided, a 'Prot_public' layer is created. "
            "Any row with a gene id is treated as present."
        ),
    )
    parser.add_argument(
        "--public_proteomics_present_column",
        default=None,
        help=(
            "Optional boolean column in the public proteomics TSV to filter present proteins "
            "(e.g. 'present_any'). If not provided, all rows with a gene id are used."
        ),
    )

    parser.add_argument(
        "--tau_threshold",
        type=float,
        default=0.95,
        help="Tau threshold for testis-specific genes (default: 0.95).",
    )
    parser.add_argument(
        "--log2fc_threshold",
        type=float,
        default=0.0,
        help="Minimum log2 FC target vs max non-target for testis-specific genes (default: 0).",
    )
    parser.add_argument(
        "--require_testis_dominance",
        action="store_true",
        help="Require is_target_max_tissue == True in GTEx table.",
    )

    parser.add_argument(
        "--sperm_tpm_threshold",
        type=float,
        default=0.1,
        help="Minimum sperm TPM threshold for sperm-present genes (default: 0.1).",
    )
    parser.add_argument(
        "--sperm_tpm_column",
        default=None,
        help=(
            "Optional TPM column to use for sperm-present filtering. If not set, the script tries "
            "sperm_tpm_mean, then sperm_tpm_median, then sperm_present_any."
        ),
    )

    parser.add_argument(
        "--internal_proteomics_present_column",
        default="prot_present_any",
        help=(
            "Column in the sperm TSV indicating internal proteomics support (default: prot_present_any). "
            "If the column is absent, the internal proteomics layer is skipped."
        ),
    )

    parser.add_argument(
        "--internal_proteomics_tsv",
        default=None,
        help=(
            "Optional internal proteomics gene-level TSV. If provided, a "
            "'Prot_internal_full' layer is created from this file."
        ),
    )
    parser.add_argument(
        "--internal_proteomics_full_present_column",
        default=None,
        help=(
            "Optional boolean-like column in the internal proteomics TSV to filter "
            "present proteins (e.g. 'prot_present_any'). If not provided, all rows "
            "with a gene id are used."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )

    return parser.parse_args()


def setup_logging(*, verbose: bool) -> None:
    """
    Configure logging.

    Parameters
    ----------
    verbose : bool
        If True, enable DEBUG logging.

    Returns
    -------
    None
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def read_tsv(*, path: str) -> pd.DataFrame:
    """
    Read a TSV file.

    Parameters
    ----------
    path : str
        Path to TSV.

    Returns
    -------
    pandas.DataFrame
        Loaded data.
    """
    return pd.read_csv(path, sep="\t", dtype=str)


def normalise_gene_key(*, series: pd.Series) -> pd.Series:
    """
    Normalise a gene key series by stripping whitespace and Ensembl version suffixes.

    Parameters
    ----------
    series : pandas.Series
        Series containing gene identifiers.

    Returns
    -------
    pandas.Series
        Normalised gene identifiers.
    """
    s = series.astype(str).str.strip()
    # Strip Ensembl version suffix if present (e.g. ENSG0000... .12)
    s = s.str.replace(r"\.\d+$", "", regex=True)
    # Treat common "nan" strings as missing
    s = s.replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})
    return s


def ensure_gene_key_column(*, df: pd.DataFrame, context: str) -> pd.DataFrame:
    """
    Ensure a DataFrame contains a 'gene_key' column by renaming common columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.
    context : str
        Label for error messages.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing a 'gene_key' column.

    Raises
    ------
    KeyError
        If no suitable gene identifier column is found.
    """
    if "gene_key" in df.columns:
        return df

    candidates = [
        "ensembl_gene_id",
        "gene_id",
        "Gene stable ID",
        "Gene stable ID version",
        "hgnc_symbol",
        "hgnc_symbol_norm",
        "symbol",
        "gene_symbol",
        "Name",
    ]

    for col in candidates:
        if col in df.columns:
            LOGGER.debug("Renaming %s -> gene_key for %s", col, context)
            return df.rename(columns={col: "gene_key"})

    raise KeyError(
        f"Could not find a gene identifier column to use as 'gene_key' in {context}. "
        f"Columns present: {list(df.columns)}"
    )


def coerce_numeric(*, series: pd.Series, context: str) -> pd.Series:
    """
    Coerce a series to numeric, safely.

    Parameters
    ----------
    series : pandas.Series
        Series to convert.
    context : str
        Description used for logging.

    Returns
    -------
    pandas.Series
        Numeric series (float), with non-parsable values as NaN.
    """
    out = pd.to_numeric(series, errors="coerce")
    if out.isna().all():
        LOGGER.warning("All values became NaN when converting to numeric for: %s", context)
    return out


def filter_gtex_testis_specific(
    *,
    gtex_df: pd.DataFrame,
    tau_threshold: float,
    log2fc_threshold: float,
    require_testis_dominance: bool,
) -> pd.DataFrame:
    """
    Filter GTEx table to define testis-specific genes.

    Parameters
    ----------
    gtex_df : pandas.DataFrame
        GTEx table with gene identifiers and specificity metrics.
    tau_threshold : float
        Minimum tau value.
    log2fc_threshold : float
        Minimum log2 FC (target vs max non-target).
    require_testis_dominance : bool
        If True, require is_target_max_tissue == True.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame.
    """
    df = gtex_df.copy()

    if "tau" not in df.columns:
        raise KeyError("GTEx TSV is missing required column: tau")
    if "log2_fc_target_vs_max_non_target" not in df.columns:
        raise KeyError("GTEx TSV is missing required column: log2_fc_target_vs_max_non_target")

    df["tau_num"] = coerce_numeric(series=df["tau"], context="GTEx tau")
    df["log2fc_num"] = coerce_numeric(
        series=df["log2_fc_target_vs_max_non_target"],
        context="GTEx log2_fc_target_vs_max_non_target",
    )

    keep = (df["tau_num"] >= tau_threshold) & (df["log2fc_num"] >= log2fc_threshold)

    if require_testis_dominance:
        if "is_target_max_tissue" not in df.columns:
            raise KeyError("GTEx TSV is missing required column for --require_testis_dominance: is_target_max_tissue")
        dom = df["is_target_max_tissue"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
        keep = keep & dom

    df = df.loc[keep].copy()
    return df


def filter_sperm_present(
    *,
    sperm_df: pd.DataFrame,
    sperm_tpm_threshold: float,
    sperm_tpm_column: Optional[str],
) -> pd.DataFrame:
    """
    Filter sperm table to define sperm-present genes.

    Parameters
    ----------
    sperm_df : pandas.DataFrame
        Sperm evidence table.
    sperm_tpm_threshold : float
        Minimum TPM threshold.
    sperm_tpm_column : str or None
        Specific TPM column to use. If None, the script tries reasonable defaults.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing sperm-present genes.
    """
    df = sperm_df.copy()

    chosen_col = sperm_tpm_column
    if chosen_col is None:
        for candidate in ["sperm_tpm_mean", "sperm_tpm_median"]:
            if candidate in df.columns:
                chosen_col = candidate
                break

    if chosen_col is not None and chosen_col in df.columns:
        df["sperm_tpm_num"] = coerce_numeric(series=df[chosen_col], context=f"Sperm {chosen_col}")
        keep = df["sperm_tpm_num"] >= sperm_tpm_threshold
        return df.loc[keep].copy()

    if "sperm_present_any" in df.columns:
        keep = df["sperm_present_any"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
        return df.loc[keep].copy()

    raise KeyError(
        "Could not determine sperm-present genes: no TPM column found and sperm_present_any is absent. "
        f"Columns present: {list(df.columns)}"
    )


def derive_proteomics_layer_from_sperm(
    *,
    sperm_df: pd.DataFrame,
    present_column: str,
) -> Optional[pd.DataFrame]:
    """
    Derive an internal proteomics layer from columns in the sperm TSV.

    Parameters
    ----------
    sperm_df : pandas.DataFrame
        Sperm evidence table (may include internal proteomics annotations).
    present_column : str
        Boolean-like column indicating proteomics presence.

    Returns
    -------
    pandas.DataFrame or None
        Filtered DataFrame for proteomics-supported genes, or None if the column is absent.
    """
    if present_column not in sperm_df.columns:
        LOGGER.info(
            "Internal proteomics column '%s' not found in sperm TSV; skipping Prot_internal layer.",
            present_column,
        )
        return None

    df = sperm_df.copy()
    keep = df[present_column].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    return df.loc[keep].copy()


def filter_present_rows(
    *,
    df: pd.DataFrame,
    present_column: Optional[str],
    context: str,
) -> pd.DataFrame:
    """
    Optionally filter a DataFrame using a boolean-like present column.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.
    present_column : str or None
        Column used to filter for present rows. If None, no filtering is applied.
    context : str
        Label used in error messages.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame.
    """
    if present_column is None:
        return df

    if present_column not in df.columns:
        raise KeyError(
            f"{context} TSV does not contain the specified present column: {present_column}. "
            f"Columns present: {list(df.columns)}"
        )

    keep = df[present_column].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    return df.loc[keep].copy()


def derive_public_proteomics_layer(
    *,
    public_df: pd.DataFrame,
    present_column: Optional[str],
) -> pd.DataFrame:
    """
    Derive a public proteomics layer from a TSV.

    Parameters
    ----------
    public_df : pandas.DataFrame
        Public proteomics table.
    present_column : str or None
        Optional boolean-like column indicating presence. If None, all rows are used.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame for public proteomics-supported genes.
    """
    df = public_df.copy()

    if present_column is None:
        return df

    if present_column not in df.columns:
        raise KeyError(
            f"Public proteomics TSV does not contain the specified present column: {present_column}. "
            f"Columns present: {list(df.columns)}"
        )

    keep = df[present_column].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    return df.loc[keep].copy()


def as_gene_set(*, df: pd.DataFrame) -> Set[str]:
    """
    Convert a DataFrame with a 'gene_key' column to a set of non-missing gene keys.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data.

    Returns
    -------
    set
        Set of gene keys.
    """
    s = normalise_gene_key(series=df["gene_key"])
    return set(s.dropna().unique().tolist())


def safe_makedirs(*, path: str) -> None:
    """
    Create a directory if it does not exist.

    Parameters
    ----------
    path : str
        Directory path.

    Returns
    -------
    None
    """
    os.makedirs(path, exist_ok=True)


def write_gene_list(*, genes: Set[str], out_path: str) -> None:
    """
    Write a list of gene identifiers to a TSV file (one column: gene_key).

    Parameters
    ----------
    genes : set
        Set of gene identifiers.
    out_path : str
        Output file path.

    Returns
    -------
    None
    """
    df = pd.DataFrame({"gene_key": sorted(genes)})
    df.to_csv(out_path, sep="\t", index=False)


def format_combo_name(*, layer_names: Iterable[str], joiner: str = "__AND__") -> str:
    """
    Format a filename-safe combination name.

    Parameters
    ----------
    layer_names : Iterable[str]
        Layer names to combine.
    joiner : str
        Join string.

    Returns
    -------
    str
        Combined name.
    """
    return joiner.join(layer_names)


def compute_pairwise_overlaps(*, layers: List[LayerSpec]) -> pd.DataFrame:
    """
    Compute pairwise overlap counts and percentages.

    Percentages are reported relative to each layer (A and B).

    Parameters
    ----------
    layers : list of LayerSpec
        Layers to compare.

    Returns
    -------
    pandas.DataFrame
        Pairwise overlap summary.
    """
    rows: List[dict] = []
    for a, b in itertools.combinations(layers, 2):
        inter = a.genes & b.genes
        a_n = len(a.genes)
        b_n = len(b.genes)
        inter_n = len(inter)

        rows.append(
            {
                "Layer_A": a.name,
                "Layer_B": b.name,
                "A_count": a_n,
                "B_count": b_n,
                "Overlap_count": inter_n,
                "Overlap_pct_of_A": (inter_n / a_n * 100.0) if a_n else 0.0,
                "Overlap_pct_of_B": (inter_n / b_n * 100.0) if b_n else 0.0,
            }
        )

    return pd.DataFrame(rows)


def compute_all_combinations(*, layers: List[LayerSpec]) -> pd.DataFrame:
    """
    Compute overlap counts and percentages for all combinations (size 1..N).

    Percentages are reported relative to each member layer in the combination.

    Parameters
    ----------
    layers : list of LayerSpec
        Layers to combine.

    Returns
    -------
    pandas.DataFrame
        Combination overlap summary.
    """
    rows: List[dict] = []
    name_to_set: Dict[str, Set[str]] = {l.name: l.genes for l in layers}

    layer_names = [l.name for l in layers]
    for r in range(1, len(layer_names) + 1):
        for combo in itertools.combinations(layer_names, r):
            combo_sets = [name_to_set[name] for name in combo]
            inter = set.intersection(*combo_sets) if combo_sets else set()
            inter_n = len(inter)

            row = {
                "Combination": format_combo_name(layer_names=combo),
                "Combination_size": r,
                "Overlap_count": inter_n,
            }

            for name in combo:
                denom = len(name_to_set[name])
                row[f"Overlap_pct_of_{name}"] = (inter_n / denom * 100.0) if denom else 0.0

            rows.append(row)

    df = pd.DataFrame(rows)
    df = df.sort_values(by=["Combination_size", "Overlap_count"], ascending=[True, False]).reset_index(drop=True)
    return df


def compute_unique_sets(*, layers: List[LayerSpec]) -> Dict[str, Set[str]]:
    """
    Compute genes unique to each layer (layer-only sets).

    Parameters
    ----------
    layers : list of LayerSpec
        Layers.

    Returns
    -------
    dict
        Mapping of layer name to its unique gene set.
    """
    name_to_set: Dict[str, Set[str]] = {l.name: l.genes for l in layers}
    unique: Dict[str, Set[str]] = {}

    for layer in layers:
        others = set()
        for other_name, other_set in name_to_set.items():
            if other_name == layer.name:
                continue
            others |= other_set
        unique[layer.name] = layer.genes - others

    return unique


def main() -> int:
    """
    Run the overlap summary workflow.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    args = parse_args()
    setup_logging(verbose=args.verbose)

    out_dir = args.out_dir if args.out_dir is not None else args.base_dir
    gene_lists_dir = os.path.join(out_dir, "gene_lists")
    safe_makedirs(path=out_dir)
    safe_makedirs(path=gene_lists_dir)

    LOGGER.info("GTEx table: %s", args.gtex_tsv)
    LOGGER.info("Sperm table: %s", args.sperm_tsv)
    if args.public_proteomics_tsv:
        LOGGER.info("Public proteomics table: %s", args.public_proteomics_tsv)

    # Load and standardise GTEx
    gtex_df = read_tsv(path=args.gtex_tsv)
    gtex_df = ensure_gene_key_column(df=gtex_df, context="GTEx TSV")
    gtex_df["gene_key"] = normalise_gene_key(series=gtex_df["gene_key"])

    # Load and standardise sperm
    sperm_df = read_tsv(path=args.sperm_tsv)
    sperm_df = ensure_gene_key_column(df=sperm_df, context="Sperm TSV")
    sperm_df["gene_key"] = normalise_gene_key(series=sperm_df["gene_key"])

    # Filter layers
    gtex_filtered = filter_gtex_testis_specific(
        gtex_df=gtex_df,
        tau_threshold=args.tau_threshold,
        log2fc_threshold=args.log2fc_threshold,
        require_testis_dominance=args.require_testis_dominance,
    )
    sperm_filtered = filter_sperm_present(
        sperm_df=sperm_df,
        sperm_tpm_threshold=args.sperm_tpm_threshold,
        sperm_tpm_column=args.sperm_tpm_column,
    )

    rna_testis_set = as_gene_set(df=gtex_filtered)
    rna_sperm_set = as_gene_set(df=sperm_filtered)

    layers: List[LayerSpec] = [
        LayerSpec(name="RNA_testis", genes=rna_testis_set),
        LayerSpec(name="RNA_sperm", genes=rna_sperm_set),
    ]

    # Internal proteomics (from sperm TSV)
    internal_df = derive_proteomics_layer_from_sperm(
        sperm_df=sperm_df,
        present_column=args.internal_proteomics_present_column,
    )
    if internal_df is not None:
        internal_df = ensure_gene_key_column(df=internal_df, context="Internal proteomics (from sperm TSV)")
        internal_df["gene_key"] = normalise_gene_key(series=internal_df["gene_key"])
        prot_internal_set = as_gene_set(df=internal_df)
        layers.append(LayerSpec(name="Prot_internal_supported_in_sperm_tsv", genes=prot_internal_set))
    # Public proteomics (separate TSV)
    # Full internal proteomics (separate gene-level TSV)
    if args.internal_proteomics_tsv:
        internal_full_df = read_tsv(path=args.internal_proteomics_tsv)
        internal_full_df = ensure_gene_key_column(
            df=internal_full_df,
            context="Internal proteomics TSV",
        )
        internal_full_df["gene_key"] = normalise_gene_key(series=internal_full_df["gene_key"])

        internal_full_df = filter_present_rows(
            df=internal_full_df,
            present_column=args.internal_proteomics_full_present_column,
            context="Internal proteomics",
        )

        prot_internal_full_set = as_gene_set(df=internal_full_df)
        layers.append(LayerSpec(name="Prot_internal_full", genes=prot_internal_full_set))
        
        
    if args.public_proteomics_tsv:
        public_df = read_tsv(path=args.public_proteomics_tsv)
        public_df = ensure_gene_key_column(df=public_df, context="Public proteomics TSV")
        public_df["gene_key"] = normalise_gene_key(series=public_df["gene_key"])
        public_filtered = derive_public_proteomics_layer(
            public_df=public_df,
            present_column=args.public_proteomics_present_column,
        )
        prot_public_set = as_gene_set(df=public_filtered)
        layers.append(LayerSpec(name="Prot_public", genes=prot_public_set))

    # Print headline counts (human-readable)
    counts = {l.name: len(l.genes) for l in layers}
    for name in sorted(counts.keys()):
        LOGGER.info("%s genes: %d", name, counts[name])

    # RNA overlap headline (matches your current outputs)
    overlap_rna = len(rna_testis_set & rna_sperm_set)
    testis_only = len(rna_testis_set - rna_sperm_set)
    sperm_only = len(rna_sperm_set - rna_testis_set)

    print(f"Testis-specific genes: {len(rna_testis_set)}")
    print(f"Sperm-present genes: {len(rna_sperm_set)}")
    print(f"Overlap (both): {overlap_rna}")
    print(f"Testis-only: {testis_only}")
    print(f"Sperm-only: {sperm_only}")

    # Write per-layer gene lists
    for layer in layers:
        out_path = os.path.join(gene_lists_dir, f"{layer.name}.tsv")
        write_gene_list(genes=layer.genes, out_path=out_path)

    # Unique (only) sets
    unique_sets = compute_unique_sets(layers=layers)
    for layer_name, gene_set in unique_sets.items():
        out_path = os.path.join(gene_lists_dir, f"{layer_name}__ONLY.tsv")
        write_gene_list(genes=gene_set, out_path=out_path)

    # Convenience alias requested (134 sperm-only genes)
    sperm_only_path = os.path.join(gene_lists_dir, "sperm_only_genes.tsv")
    write_gene_list(genes=unique_sets["RNA_sperm"], out_path=sperm_only_path)

    # Pairwise overlaps
    pairwise_df = compute_pairwise_overlaps(layers=layers)
    pairwise_out = os.path.join(out_dir, "pairwise_overlaps.tsv")
    pairwise_df.to_csv(pairwise_out, sep="\t", index=False)

    # Write pairwise overlap gene lists
    for a, b in itertools.combinations(layers, 2):
        name = format_combo_name(layer_names=[a.name, b.name])
        out_path = os.path.join(gene_lists_dir, f"{name}.tsv")
        write_gene_list(genes=a.genes & b.genes, out_path=out_path)

    # All combinations (1..N)
    combos_df = compute_all_combinations(layers=layers)
    combos_out = os.path.join(out_dir, "combination_overlaps.tsv")
    combos_df.to_csv(combos_out, sep="\t", index=False)

    # Write triple+ overlap gene lists (and also full intersection)
    # (Pairwise already written above; unique sets already written.)
    layer_names = [l.name for l in layers]
    name_to_set = {l.name: l.genes for l in layers}

    for r in range(3, len(layer_names) + 1):
        for combo in itertools.combinations(layer_names, r):
            inter = set.intersection(*[name_to_set[n] for n in combo])
            name = format_combo_name(layer_names=combo)
            out_path = os.path.join(gene_lists_dir, f"{name}.tsv")
            write_gene_list(genes=inter, out_path=out_path)

    # Layer counts table
    layer_counts_df = pd.DataFrame(
        [{"Layer": l.name, "Gene_count": len(l.genes)} for l in layers]
    ).sort_values(by="Gene_count", ascending=False)
    layer_counts_out = os.path.join(out_dir, "layer_counts.tsv")
    layer_counts_df.to_csv(layer_counts_out, sep="\t", index=False)

    LOGGER.info("Wrote: %s", layer_counts_out)
    LOGGER.info("Wrote: %s", pairwise_out)
    LOGGER.info("Wrote: %s", combos_out)
    LOGGER.info("Gene lists directory: %s", gene_lists_dir)
    LOGGER.info("Sperm-only genes (GO-ready): %s", sperm_only_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
