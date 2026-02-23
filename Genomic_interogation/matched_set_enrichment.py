#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
matched_set_enrichment.py

Permutation testing using matched random gene sets to assess whether a target
gene set is enriched for particular genomic-context features.

This script consumes the per-gene feature table produced by
build_gene_context_table.py and performs the following:

- Defines a target set using a boolean column (e.g., is_sperm_gene) or by a list
  of gene_ids / gene_keys (optional modes).
- Generates N matched random control sets (default 10,000) using matching strata
  (by default chromosome + gene length bins; optional GC bins).
- Computes an observed effect size for each feature as the difference in medians
  between target and one matched control sample.
- Builds null distributions of effect sizes across permutations and calculates
  empirical two-sided p-values.
- Applies Benjamini–Hochberg FDR correction using statsmodels.

This script does not create plots.

Example
python matched_set_enrichment.py \
  --features_tsv gene_context_features.tsv \
  --target_gene_list_tsv sperm_genes.tsv \
  --gene_key_column gene_key \
  --n_permutations 10000 \
  --n_length_bins 10 \
  --match_on_gc \
  --n_gc_bins 10 \
  --out_prefix sperm_vs_null
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def setup_logger(*, verbose: bool) -> None:
    """Configure logging for the script."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_tsv(*, path: str) -> pd.DataFrame:
    """Read a TSV file into a DataFrame."""
    return pd.read_csv(path, sep="\t", low_memory=False)


def write_tsv(*, df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame as TSV, creating parent directories if needed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def read_gene_list(*, path: str, gene_key_column: str) -> pd.Series:
    """Read a gene list TSV and return normalised gene keys."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if gene_key_column not in df.columns:
        raise ValueError(f"Gene list TSV must contain column '{gene_key_column}'.")
    keys = df[gene_key_column].astype(str).str.strip()
    keys = keys[keys != ""].drop_duplicates()
    return keys


def normalise_gene_key(*, s: pd.Series) -> pd.Series:
    """Normalise gene keys for matching (upper-case, trimmed)."""
    return s.astype(str).str.strip().str.upper()


def add_target_flag(
    *, features: pd.DataFrame, gene_list_tsv: str, gene_key_column: str
) -> pd.DataFrame:
    """
    Add a boolean 'is_target' column based on presence in a gene list TSV.

    Matching is performed using 'gene_key' in the features table by default;
    if not present, we fall back to 'gene_name' if available.
    """
    keys = normalise_gene_key(s=read_gene_list(path=gene_list_tsv, gene_key_column=gene_key_column))

    df = features.copy()
    if "gene_key" in df.columns:
        feature_key = normalise_gene_key(s=df["gene_key"])
    elif "gene_name" in df.columns:
        feature_key = normalise_gene_key(s=df["gene_name"])
    else:
        raise ValueError("Features table must contain 'gene_key' or 'gene_name' to match targets.")

    df["is_target"] = feature_key.isin(set(keys.to_list()))
    return df


def make_bins(
    *, values: pd.Series, n_bins: int, label_prefix: str
) -> pd.Series:
    """
    Create quantile bins for a numeric Series.

    Falls back to fewer bins if there are not enough distinct values.
    """
    v = pd.to_numeric(values, errors="coerce")
    v = v.replace([np.inf, -np.inf], np.nan)
    if v.notna().sum() == 0:
        return pd.Series(["bin_na"] * len(values), index=values.index)

    n_unique = int(v.dropna().nunique())
    use_bins = int(min(max(2, n_bins), max(2, n_unique)))
    try:
        b = pd.qcut(v, q=use_bins, duplicates="drop")
    except ValueError:
        return pd.Series(["bin_na"] * len(values), index=values.index)

    return b.astype(str).radd(f"{label_prefix}:")


def prepare_matching_strata(
    *,
    df: pd.DataFrame,
    n_length_bins: int,
    match_on_gc: bool,
    n_gc_bins: int,
) -> pd.DataFrame:
    """Add matching strata columns used for control sampling."""
    out = df.copy()

    if "Chromosome" not in out.columns:
        raise ValueError("Features table must contain a 'Chromosome' column for matching.")

    if "gene_length_bp" not in out.columns:
        raise ValueError("Features table must contain 'gene_length_bp'.")

    # Length bins on log scale to stabilise long-tailed distributions
    gene_len = pd.to_numeric(out["gene_length_bp"], errors="coerce").clip(lower=1)
    out["gene_length_log10"] = np.log10(gene_len.astype(float))
    out["len_bin"] = make_bins(values=out["gene_length_log10"], n_bins=n_length_bins, label_prefix="len")

    if match_on_gc:
        if "gc_gene_body_pct" not in out.columns:
            raise ValueError("match_on_gc requested but 'gc_gene_body_pct' is missing.")
        out["gc_bin"] = make_bins(values=out["gc_gene_body_pct"], n_bins=n_gc_bins, label_prefix="gc")
    else:
        out["gc_bin"] = "gc:all"

    out["match_key"] = (
        out["Chromosome"].astype(str) + "|" + out["len_bin"].astype(str) + "|" + out["gc_bin"].astype(str)
    )
    return out


def choose_controls_for_targets(
    *,
    df: pd.DataFrame,
    rng: np.random.Generator,
    allow_replacement: bool,
) -> pd.Index:
    """
    Sample a matched control gene for each target gene using match_key strata.

    Controls are drawn from non-target genes in the same stratum.
    If a stratum is empty, the gene is skipped and will contribute NaNs.

    Returns the index (row labels) of selected controls, aligned to target rows.
    """
    targets = df[df["is_target"]].copy()
    background = df[~df["is_target"]].copy()

    if targets.empty:
        raise ValueError("No target genes found (is_target is empty).")

    # Build pools by match_key
    pools: Dict[str, np.ndarray] = {}
    for key, sub in background.groupby("match_key"):
        pools[str(key)] = sub.index.to_numpy()

    chosen: List[int] = []
    used_by_key: Dict[str, set] = {}

    for row in targets.itertuples(index=True):
        key = str(getattr(row, "match_key"))
        pool = pools.get(key, np.array([], dtype=int))

        if pool.size == 0:
            chosen.append(-1)
            continue

        if allow_replacement:
            chosen_idx = int(rng.choice(pool, size=1)[0])
            chosen.append(chosen_idx)
            continue

        used = used_by_key.setdefault(key, set())
        available = np.array([i for i in pool if int(i) not in used], dtype=int)

        if available.size == 0:
            chosen.append(int(rng.choice(pool, size=1)[0]))
        else:
            pick = int(rng.choice(available, size=1)[0])
            used.add(pick)
            chosen.append(pick)

    return pd.Index(chosen, name="control_index")


def get_numeric_feature_columns(*, df: pd.DataFrame) -> List[str]:
    """
    Identify numeric feature columns to test.

    Excludes obvious identifier columns and matching helper columns.
    """
    exclude = {
        "Start",
        "End",
        "tss_0based",
        "gene_length_bp",
        "gene_length_log10",
    }
    exclude_prefixes = (
        "is_within_",
    )
    exclude_cols = {
        "gene_id",
        "gene_name",
        "gene_key",
        "Chromosome",
        "Strand",
        "is_target",
        "len_bin",
        "gc_bin",
        "match_key",
    }
    cols: List[str] = []
    for c in df.columns:
        if c in exclude_cols or c in exclude:
            continue
        if any(c.startswith(p) for p in exclude_prefixes):
            # flags are boolean; treat separately if wanted later
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            cols.append(str(c))
    return cols


def median_effect(*, target: np.ndarray, control: np.ndarray) -> float:
    """Compute median(target) - median(control), ignoring NaNs."""
    t = target[np.isfinite(target)]
    c = control[np.isfinite(control)]
    if t.size == 0 or c.size == 0:
        return float("nan")
    return float(np.median(t) - np.median(c))


def empirical_p_value_two_sided(*, observed: float, null: np.ndarray) -> float:
    """Compute two-sided empirical p-value with +1 smoothing."""
    null_finite = null[np.isfinite(null)]
    if not np.isfinite(observed) or null_finite.size == 0:
        return float("nan")
    ge = float(np.sum(np.abs(null_finite) >= abs(observed)))
    return float((1.0 + ge) / (1.0 + float(null_finite.size)))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Matched random-set permutation testing for gene context features.")
    parser.add_argument("--features_tsv", required=True, help="TSV produced by build_gene_context_table.py.")
    parser.add_argument("--target_gene_list_tsv", required=True, help="TSV with a column listing target genes.")
    parser.add_argument("--gene_key_column", default="gene_key", help="Column name in target list TSV (default: gene_key).")
    parser.add_argument("--n_permutations", type=int, default=10_000, help="Number of matched sets (default: 10000).")
    parser.add_argument("--seed", type=int, default=1, help="Random seed (default: 1).")
    parser.add_argument("--n_length_bins", type=int, default=10, help="Number of gene-length bins (default: 10).")
    parser.add_argument("--match_on_gc", action="store_true", help="Also match on GC% bins.")
    parser.add_argument("--n_gc_bins", type=int, default=10, help="Number of GC bins if --match_on_gc is set (default: 10).")
    parser.add_argument("--allow_replacement", action="store_true", help="Allow control reuse within strata (default: off).")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for TSV files.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv=argv)
    setup_logger(verbose=bool(args.verbose))

    logging.info("Loading features TSV: %s", args.features_tsv)
    features = read_tsv(path=args.features_tsv)

    logging.info("Marking targets using gene list: %s", args.target_gene_list_tsv)
    features = add_target_flag(
        features=features,
        gene_list_tsv=args.target_gene_list_tsv,
        gene_key_column=str(args.gene_key_column),
    )

    n_targets = int(features["is_target"].sum())
    if n_targets == 0:
        logging.error("No target genes found after matching. Check gene keys and mapping.")
        return 2
    logging.info("Target genes found: %s", n_targets)

    logging.info("Preparing matching strata")
    features = prepare_matching_strata(
        df=features,
        n_length_bins=int(args.n_length_bins),
        match_on_gc=bool(args.match_on_gc),
        n_gc_bins=int(args.n_gc_bins),
    )

    numeric_cols = get_numeric_feature_columns(df=features)
    if not numeric_cols:
        logging.error("No numeric feature columns found to test.")
        return 2
    logging.info("Testing %s numeric features", len(numeric_cols))

    rng = np.random.default_rng(int(args.seed))

    # Observed: one matched control set with a fixed seed offset for reproducibility
    obs_rng = np.random.default_rng(int(args.seed) + 999)
    obs_controls = choose_controls_for_targets(
        df=features, rng=obs_rng, allow_replacement=bool(args.allow_replacement)
    )
    target_df = features[features["is_target"]].copy()
    control_df_obs = features.reindex(obs_controls).copy()
    missing_mask = ~pd.Index(ctrl_idx).isin(features.index)
    if missing_mask.any():
        logging.warning(
            "Missing controls for %s/%s targets (no match in background strata).",
            int(missing_mask.sum()),
            int(missing_mask.size),
        )
    ctrl_df.loc[missing_mask, numeric_cols] = np.nan
    control_df_obs.index = target_df.index  # align row-wise (may include -1 rows as NaN)

    # Replace invalid selections (-1) with NaNs
    invalid_mask = obs_controls.to_numpy() == -1
    if invalid_mask.any():
        control_df_obs.loc[invalid_mask, numeric_cols] = np.nan

    observed_effects: Dict[str, float] = {}
    for col in numeric_cols:
        observed_effects[col] = median_effect(
            target=target_df[col].to_numpy(dtype=float),
            control=control_df_obs[col].to_numpy(dtype=float),
        )

    # Null distributions
    logging.info("Generating %s matched permutations", int(args.n_permutations))
    null_effects = {col: np.full(int(args.n_permutations), np.nan, dtype=float) for col in numeric_cols}

    for i in range(int(args.n_permutations)):
        ctrl_idx = choose_controls_for_targets(
            df=features, rng=rng, allow_replacement=bool(args.allow_replacement)
        )
        #ctrl_df = features.loc[ctrl_idx].copy()
        ctrl_df = features.reindex(ctrl_idx).copy()
        ctrl_df.index = target_df.index  # align row-wise (may include -1 rows as NaN)

        invalid = ctrl_idx.to_numpy() == -1
        if invalid.any():
            ctrl_df.loc[invalid, numeric_cols] = np.nan

        for col in numeric_cols:
            null_effects[col][i] = median_effect(
                target=target_df[col].to_numpy(dtype=float),
                control=ctrl_df[col].to_numpy(dtype=float),
            )

        if (i + 1) % 1000 == 0:
            logging.info("Completed %s/%s permutations", i + 1, int(args.n_permutations))

    # Summarise and p-values
    rows: List[Dict[str, object]] = []
    for col in numeric_cols:
        obs = float(observed_effects[col])
        null = null_effects[col]
        p = empirical_p_value_two_sided(observed=obs, null=null)
        rows.append(
            {
                "feature": col,
                "observed_median_diff": obs,
                "null_median_diff_mean": float(np.nanmean(null)),
                "null_median_diff_sd": float(np.nanstd(null)),
                "empirical_p_two_sided": p,
                "n_targets": n_targets,
                "n_permutations": int(args.n_permutations),
                "match_on_gc": bool(args.match_on_gc),
                "n_length_bins": int(args.n_length_bins),
                "n_gc_bins": int(args.n_gc_bins) if bool(args.match_on_gc) else 0,
                "allow_replacement": bool(args.allow_replacement),
            }
        )

    summary = pd.DataFrame(rows).sort_values(by="empirical_p_two_sided", ascending=True).reset_index(drop=True)

    # BH correction
    pvals = summary["empirical_p_two_sided"].to_numpy(dtype=float)
    ok = np.isfinite(pvals)
    qvals = np.full_like(pvals, np.nan, dtype=float)
    if ok.any():
        _, q, _, _ = multipletests(pvals[ok], method="fdr_bh")
        qvals[ok] = q
    summary["bh_q_value"] = qvals

    out_summary = f"{args.out_prefix}.permutation_summary.tsv"
    out_summary_bh = f"{args.out_prefix}.permutation_summary_bh.tsv"

    logging.info("Writing summary TSV: %s", out_summary)
    write_tsv(df=summary.drop(columns=["bh_q_value"]), path=out_summary)

    logging.info("Writing BH-adjusted TSV: %s", out_summary_bh)
    write_tsv(df=summary, path=out_summary_bh)

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
