#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genome_organisation_tests.py

Chromosome enrichment + nearest-neighbour clustering tests for a target gene set,
plus a chromosome ideogram-style plot showing target gene locations.

Inputs
------
1) --features_tsv
   TSV produced by your genomic-context pipeline (e.g. build_gene_context_table.py)
   Must include at least:
     - Chromosome
     - Start
     - End
     - gene_length_bp
   And one of:
     - gene_key (preferred)
     - gene_name (fallback)

2) --target_gene_list_tsv
   TSV containing a column of target genes (e.g. sperm-expressed genes)

Outputs (TSV, not comma-separated)
---------------------------------
- <out_prefix>.chromosome_enrichment.tsv
- <out_prefix>.chromosome_enrichment_bh.tsv
- <out_prefix>.nearest_neighbour_summary.tsv
- <out_prefix>.nearest_neighbour_null_summary.tsv
- <out_prefix>.ideogram.png
- <out_prefix>.ideogram.pdf

What it tests
-------------
A) Chromosome enrichment:
   For each chromosome, tests whether the target set is over/under-represented
   compared with the background gene universe using Fisher's exact test and BH FDR.

B) Nearest-neighbour clustering ("gene islands"):
   Computes, for each target gene, the distance (bp) to the nearest OTHER target
   gene on the same chromosome (based on gene midpoint).
   Then tests whether the target set has unusually small nearest-neighbour
   distances compared with matched random gene sets (matched by chromosome and
   gene length bins; optional GC bins), using permutation testing.

Notes
-----
- This script does not assume anything about pathogens; it is generic genome organisation testing.
- Outputs are designed to drop into your Results/Supplementary materials directly.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
from dataclasses import dataclass
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


def normalise_gene_key(*, s: pd.Series) -> pd.Series:
    """Normalise gene keys for matching (upper-case, trimmed)."""
    return s.astype(str).str.strip().str.upper()


def read_gene_list(*, path: str, gene_key_column: str) -> pd.Series:
    """Read a gene list TSV and return normalised gene keys."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if gene_key_column not in df.columns:
        raise ValueError(f"Target gene list TSV must contain column '{gene_key_column}'.")
    keys = df[gene_key_column].astype(str).str.strip()
    keys = keys[keys != ""].drop_duplicates()
    return keys


def add_target_flag(
    *, features: pd.DataFrame, gene_list_tsv: str, gene_key_column: str
) -> pd.DataFrame:
    """
    Add boolean 'is_target' column based on presence in a gene list TSV.

    Matching uses 'gene_key' in features if present; otherwise falls back to 'gene_name'.
    """
    keys = normalise_gene_key(s=read_gene_list(path=gene_list_tsv, gene_key_column=gene_key_column))
    key_set = set(keys.to_list())

    df = features.copy()

    if "gene_key" in df.columns:
        feature_key = normalise_gene_key(s=df["gene_key"])
    elif "gene_name" in df.columns:
        feature_key = normalise_gene_key(s=df["gene_name"])
    else:
        raise ValueError("Features table must contain 'gene_key' or 'gene_name' to match targets.")

    df["is_target"] = feature_key.isin(key_set)
    return df


def make_bins(*, values: pd.Series, n_bins: int, label_prefix: str) -> pd.Series:
    """
    Create quantile bins for a numeric Series.

    Falls back to fewer bins if there are not enough distinct values.
    """
    v = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
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
    """Add matching strata columns used for permutation sampling."""
    out = df.copy()

    required_cols = ["Chromosome", "gene_length_bp"]
    for col in required_cols:
        if col not in out.columns:
            raise ValueError(f"Features table must contain '{col}'.")

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


def chromosome_sort_key(*, chrom: str) -> Tuple[int, str]:
    """
    Sort key for chromosomes: 1..22, X, Y, MT/M, then others.

    This is cosmetic (tables/plots) and does not affect statistical tests.
    """
    c = str(chrom).replace("chr", "").replace("CHR", "").strip()
    if c.isdigit():
        return (0, f"{int(c):02d}")
    if c in {"X"}:
        return (1, "X")
    if c in {"Y"}:
        return (2, "Y")
    if c in {"MT", "M", "MITO"}:
        return (3, "MT")
    return (4, c)


def fisher_exact_test(*, a: int, b: int, c: int, d: int) -> float:
    """
    Two-sided Fisher's exact test p-value for 2x2 table:
        [[a, b],
         [c, d]]

    Uses scipy if available; otherwise falls back to hypergeometric enumeration.
    """
    try:
        from scipy.stats import fisher_exact  # type: ignore

        _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        return float(p)
    except Exception:
        # Fallback: compute two-sided p by summing probabilities <= observed
        # under the hypergeometric distribution (fixed margins).
        # This is slower but fine per chromosome (few tests).
        n1 = a + b
        n2 = c + d
        m1 = a + c
        m2 = b + d
        n = n1 + n2

        def log_choose(n_: int, k_: int) -> float:
            return math.lgamma(n_ + 1) - math.lgamma(k_ + 1) - math.lgamma(n_ - k_ + 1)

        def pmf(x: int) -> float:
            # Hypergeometric probability of x successes in sample of size n1
            # given m1 successes in population of size n
            return math.exp(log_choose(m1, x) + log_choose(m2, n1 - x) - log_choose(n, n1))

        obs = a
        lo = max(0, n1 - m2)
        hi = min(n1, m1)

        p_obs = pmf(obs)
        p_two_sided = 0.0
        for x in range(lo, hi + 1):
            px = pmf(x)
            if px <= p_obs + 1e-15:
                p_two_sided += px
        return float(min(1.0, p_two_sided))


def run_chromosome_enrichment(*, df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-chromosome enrichment of targets using Fisher's exact test + BH FDR."""
    if "Chromosome" not in df.columns or "is_target" not in df.columns:
        raise ValueError("DataFrame must contain 'Chromosome' and 'is_target'.")

    total_targets = int(df["is_target"].sum())
    total_background = int((~df["is_target"]).sum())

    rows: List[Dict[str, object]] = []
    for chrom, sub in df.groupby("Chromosome"):
        in_chr_targets = int(sub["is_target"].sum())
        in_chr_background = int((~sub["is_target"]).sum())

        a = in_chr_targets
        b = total_targets - in_chr_targets
        c = in_chr_background
        d = total_background - in_chr_background

        p = fisher_exact_test(a=a, b=b, c=c, d=d)
        frac_targets = a / total_targets if total_targets > 0 else float("nan")
        frac_genome = (a + c) / (total_targets + total_background) if (total_targets + total_background) > 0 else float("nan")

        rows.append(
            {
                "Chromosome": str(chrom),
                "targets_on_chromosome": a,
                "targets_total": total_targets,
                "background_on_chromosome": c,
                "background_total": total_background,
                "fraction_targets_on_chromosome": frac_targets,
                "fraction_all_genes_on_chromosome": frac_genome,
                "fisher_p_two_sided": p,
            }
        )

    out = pd.DataFrame(rows)
    out = out.sort_values(by="Chromosome", key=lambda s: s.map(lambda x: chromosome_sort_key(chrom=x))).reset_index(drop=True)

    pvals = out["fisher_p_two_sided"].to_numpy(dtype=float)
    ok = np.isfinite(pvals)
    qvals = np.full_like(pvals, np.nan, dtype=float)
    if ok.any():
        _, q, _, _ = multipletests(pvals[ok], method="fdr_bh")
        qvals[ok] = q
    out["bh_q_value"] = qvals

    return out


def gene_midpoint_bp(*, start: pd.Series, end: pd.Series) -> pd.Series:
    """Compute gene midpoint (bp) from Start/End columns."""
    s = pd.to_numeric(start, errors="coerce")
    e = pd.to_numeric(end, errors="coerce")
    return ((s + e) / 2.0).round().astype("Int64")


def nearest_neighbour_distances_for_targets(*, df: pd.DataFrame) -> pd.Series:
    """
    Compute nearest-neighbour distances (bp) among target genes only.

    For each target gene, distance is to the nearest OTHER target gene on the same chromosome,
    based on gene midpoint. If a chromosome has <2 targets, distances are NaN for those genes.
    """
    required = ["Chromosome", "Start", "End", "is_target"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Features table must contain '{col}' for nearest-neighbour analysis.")

    work = df.copy()
    work["mid_bp"] = gene_midpoint_bp(start=work["Start"], end=work["End"])
    work = work[work["is_target"]].copy()

    dists: List[pd.Series] = []
    for chrom, sub in work.groupby("Chromosome"):
        sub = sub.dropna(subset=["mid_bp"]).sort_values(by="mid_bp").copy()
        if sub.shape[0] < 2:
            dists.append(pd.Series([np.nan] * sub.shape[0], index=sub.index))
            continue

        pos = sub["mid_bp"].astype(float).to_numpy()
        left = np.r_[np.nan, np.diff(pos)]
        right = np.r_[np.diff(pos), np.nan]
        nn = np.nanmin(np.vstack([left, right]), axis=0)
        dists.append(pd.Series(nn, index=sub.index))

    if not dists:
        return pd.Series(dtype=float)

    return pd.concat(dists).sort_index()


def sample_matched_gene_set_indices(
    *,
    df: pd.DataFrame,
    rng: np.random.Generator,
    allow_replacement: bool,
) -> pd.Index:
    """
    Sample a matched random gene set with the same size and stratum composition as targets.

    Matching is done via df['match_key'] and draws from non-target genes.
    """
    targets = df[df["is_target"]].copy()
    background = df[~df["is_target"]].copy()

    if targets.empty:
        raise ValueError("No target genes found (is_target is empty).")

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
            chosen.append(int(rng.choice(pool, size=1)[0]))
            continue

        used = used_by_key.setdefault(key, set())
        available = np.array([i for i in pool if int(i) not in used], dtype=int)

        if available.size == 0:
            chosen.append(int(rng.choice(pool, size=1)[0]))
        else:
            pick = int(rng.choice(available, size=1)[0])
            used.add(pick)
            chosen.append(pick)

    return pd.Index(chosen, name="sample_index")


@dataclass
class NearestNeighbourTestResult:
    """Container for nearest-neighbour permutation test outputs."""
    observed_median_nn_bp: float
    observed_mean_nn_bp: float
    empirical_p_two_sided_median: float
    empirical_p_two_sided_mean: float


def empirical_p_value_two_sided(*, observed: float, null: np.ndarray) -> float:
    """Compute two-sided empirical p-value with +1 smoothing."""
    null_finite = null[np.isfinite(null)]
    if not np.isfinite(observed) or null_finite.size == 0:
        return float("nan")
    ge = float(np.sum(np.abs(null_finite - np.nanmedian(null_finite)) >= abs(observed - np.nanmedian(null_finite))))
    return float((1.0 + ge) / (1.0 + float(null_finite.size)))


def run_nearest_neighbour_permutation_test(
    *,
    df: pd.DataFrame,
    n_permutations: int,
    seed: int,
    allow_replacement: bool,
) -> Tuple[NearestNeighbourTestResult, pd.DataFrame]:
    """
    Permutation test for target nearest-neighbour distances using matched random gene sets.

    Returns:
      - summary result object
      - null summary DataFrame (one row per permutation with median/mean NN distance)
    """
    rng = np.random.default_rng(int(seed))

    obs_nn = nearest_neighbour_distances_for_targets(df=df)
    obs_median = float(np.nanmedian(obs_nn.to_numpy(dtype=float)))
    obs_mean = float(np.nanmean(obs_nn.to_numpy(dtype=float)))

    null_medians = np.full(int(n_permutations), np.nan, dtype=float)
    null_means = np.full(int(n_permutations), np.nan, dtype=float)

    # Precompute positions for the whole genome once
    base = df.copy()
    base["mid_bp"] = gene_midpoint_bp(start=base["Start"], end=base["End"])

    target_size = int(base["is_target"].sum())
    if target_size < 2:
        raise ValueError("Nearest-neighbour test requires at least 2 target genes.")

    for i in range(int(n_permutations)):
        idx = sample_matched_gene_set_indices(df=base, rng=rng, allow_replacement=allow_replacement)
        sampled = base.reindex(idx).copy()

        # -1 indicates no available match in that stratum
        valid_mask = idx.to_numpy() != -1
        sampled = sampled.loc[valid_mask].copy()

        # Make this sampled set behave like targets for the NN computation
        sampled["is_target"] = True

        nn = nearest_neighbour_distances_for_targets(df=sampled)
        null_medians[i] = float(np.nanmedian(nn.to_numpy(dtype=float)))
        null_means[i] = float(np.nanmean(nn.to_numpy(dtype=float)))

        if (i + 1) % 1000 == 0:
            logging.info("Nearest-neighbour permutations completed %s/%s", i + 1, int(n_permutations))

    # Two-sided empirical p-values relative to the null distributions
    # For clustering we mainly care about unusually SMALL distances. Two-sided is safer for publication.
    p_median = float((1.0 + np.sum(np.abs(null_medians - np.nanmedian(null_medians)) >= abs(obs_median - np.nanmedian(null_medians))))
                     / (1.0 + np.sum(np.isfinite(null_medians))))
    p_mean = float((1.0 + np.sum(np.abs(null_means - np.nanmedian(null_means)) >= abs(obs_mean - np.nanmedian(null_means))))
                   / (1.0 + np.sum(np.isfinite(null_means))))

    result = NearestNeighbourTestResult(
        observed_median_nn_bp=obs_median,
        observed_mean_nn_bp=obs_mean,
        empirical_p_two_sided_median=p_median,
        empirical_p_two_sided_mean=p_mean,
    )

    null_df = pd.DataFrame(
        {
            "permutation": np.arange(1, int(n_permutations) + 1, dtype=int),
            "null_median_nn_bp": null_medians,
            "null_mean_nn_bp": null_means,
        }
    )

    return result, null_df


def plot_ideogram(
    *,
    df: pd.DataFrame,
    out_png: str,
    out_pdf: str,
    title: str,
) -> None:
    """
    Create a simple chromosome ideogram-style plot with target gene positions.

    Chromosomes are shown as horizontal lines (length = max End).
    Target genes are plotted as tick marks at their midpoints.
    """
    import matplotlib.pyplot as plt  # local import to keep CLI lightweight

    required = ["Chromosome", "Start", "End", "is_target"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Features table must contain '{col}' for ideogram plotting.")

    work = df.copy()
    work["mid_bp"] = gene_midpoint_bp(start=work["Start"], end=work["End"])
    work = work.dropna(subset=["mid_bp"]).copy()

    chrom_lengths = (
        work.groupby("Chromosome")["End"]
        .apply(lambda s: pd.to_numeric(s, errors="coerce").max())
        .dropna()
        .to_dict()
    )

    chroms = sorted(chrom_lengths.keys(), key=lambda x: chromosome_sort_key(chrom=str(x)))
    chrom_to_y = {c: i for i, c in enumerate(chroms)}

    fig = plt.figure(figsize=(12, max(4, int(len(chroms) * 0.35))))
    ax = fig.add_subplot(111)

    # Draw chromosome baselines
    for c in chroms:
        y = chrom_to_y[c]
        ax.hlines(y=y, xmin=0, xmax=float(chrom_lengths[c]), linewidth=2)

    # Plot target ticks
    targets = work[work["is_target"]].copy()
    for c, sub in targets.groupby("Chromosome"):
        if c not in chrom_to_y:
            continue
        y = chrom_to_y[c]
        x = sub["mid_bp"].astype(float).to_numpy()
        ax.plot(x, np.full_like(x, y, dtype=float), marker="|", linestyle="None", markersize=6)

    ax.set_yticks(list(chrom_to_y.values()))
    ax.set_yticklabels([str(c) for c in chroms])
    ax.set_xlabel("Genomic position (bp)")
    ax.set_title(title)
    ax.grid(True, axis="x", linewidth=0.5, alpha=0.4)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_png)) or ".", exist_ok=True)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Chromosome enrichment + nearest-neighbour clustering for target gene sets."
    )
    parser.add_argument("--features_tsv", required=True, help="TSV produced by build_gene_context_table.py.")
    parser.add_argument("--target_gene_list_tsv", required=True, help="TSV listing target genes.")
    parser.add_argument(
        "--gene_key_column",
        default="gene_key",
        help="Column name in target list TSV (default: gene_key).",
    )
    parser.add_argument("--out_prefix", required=True, help="Output prefix for all results.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed (default: 1).")
    parser.add_argument(
        "--n_permutations",
        type=int,
        default=10_000,
        help="Number of permutations for nearest-neighbour test (default: 10000).",
    )
    parser.add_argument(
        "--n_length_bins",
        type=int,
        default=10,
        help="Number of gene-length bins for matching (default: 10).",
    )
    parser.add_argument("--match_on_gc", action="store_true", help="Also match on GC bins.")
    parser.add_argument(
        "--n_gc_bins",
        type=int,
        default=10,
        help="Number of GC bins if --match_on_gc is set (default: 10).",
    )
    parser.add_argument(
        "--allow_replacement",
        action="store_true",
        help="Allow control reuse within strata (default: off).",
    )
    parser.add_argument("--plot_title", default="", help="Optional title override for ideogram.")
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

    # Chromosome enrichment (does not match on chromosome; it tests chromosome distribution)
    logging.info("Running chromosome enrichment tests")
    chr_enrich = run_chromosome_enrichment(df=features)

    out_chr = f"{args.out_prefix}.chromosome_enrichment.tsv"
    out_chr_bh = f"{args.out_prefix}.chromosome_enrichment_bh.tsv"
    write_tsv(df=chr_enrich.drop(columns=["bh_q_value"]), path=out_chr)
    write_tsv(df=chr_enrich, path=out_chr_bh)
    logging.info("Wrote chromosome enrichment TSVs: %s ; %s", out_chr, out_chr_bh)

    # Prepare matching strata for nearest-neighbour permutations
    logging.info("Preparing matching strata for nearest-neighbour permutation test")
    features_match = prepare_matching_strata(
        df=features,
        n_length_bins=int(args.n_length_bins),
        match_on_gc=bool(args.match_on_gc),
        n_gc_bins=int(args.n_gc_bins),
    )

    logging.info("Running nearest-neighbour clustering permutation test (%s permutations)", int(args.n_permutations))
    nn_result, nn_null_df = run_nearest_neighbour_permutation_test(
        df=features_match,
        n_permutations=int(args.n_permutations),
        seed=int(args.seed),
        allow_replacement=bool(args.allow_replacement),
    )

    out_nn = f"{args.out_prefix}.nearest_neighbour_summary.tsv"
    out_nn_null = f"{args.out_prefix}.nearest_neighbour_null_summary.tsv"

    nn_summary_df = pd.DataFrame(
        [
            {
                "n_targets": n_targets,
                "n_permutations": int(args.n_permutations),
                "n_length_bins": int(args.n_length_bins),
                "match_on_gc": bool(args.match_on_gc),
                "n_gc_bins": int(args.n_gc_bins) if bool(args.match_on_gc) else 0,
                "allow_replacement": bool(args.allow_replacement),
                "observed_median_nearest_neighbour_bp": nn_result.observed_median_nn_bp,
                "observed_mean_nearest_neighbour_bp": nn_result.observed_mean_nn_bp,
                "empirical_p_two_sided_median": nn_result.empirical_p_two_sided_median,
                "empirical_p_two_sided_mean": nn_result.empirical_p_two_sided_mean,
            }
        ]
    )
    write_tsv(df=nn_summary_df, path=out_nn)
    write_tsv(df=nn_null_df, path=out_nn_null)
    logging.info("Wrote nearest-neighbour TSVs: %s ; %s", out_nn, out_nn_null)

    # Ideogram plot
    out_png = f"{args.out_prefix}.ideogram.png"
    out_pdf = f"{args.out_prefix}.ideogram.pdf"
    title = args.plot_title.strip()
    if not title:
        title = f"Target gene distribution across chromosomes (n={n_targets})"

    logging.info("Creating ideogram plot")
    plot_ideogram(df=features, out_png=out_png, out_pdf=out_pdf, title=title)
    logging.info("Wrote ideogram: %s ; %s", out_png, out_pdf)

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
