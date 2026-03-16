#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarise_overlap_results.py

Create a concise, reproducible summary of multi-omic overlap outputs produced by
summerise_number_in_RNA.py.

This script reads:
- layer_counts.tsv
- pairwise_overlaps.tsv
- combination_overlaps.tsv

and writes:
- overlap_summary_report.md (human-readable report you can paste into Results)
- overlap_summary_key_numbers.tsv (key overlaps and percentages)
- overlap_core_sets.tsv (core 4-way/5-way intersection rows, if present)

All outputs are tab-separated where applicable.

Example
-------
python summarise_overlap_results.py \
  --layer_counts_tsv layer_counts.tsv \
  --pairwise_overlaps_tsv pairwise_overlaps.tsv \
  --combination_overlaps_tsv combination_overlaps.tsv \
  --out_dir .

Notes
-----
- Percentages are assumed to already be in the input TSVs (as produced by
  summerise_number_in_RNA.py).
- Jaccard indices are computed here from A_count, B_count, and Overlap_count.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Optional, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise overlap outputs from summerise_number_in_RNA.py."
    )
    parser.add_argument(
        "--layer_counts_tsv",
        required=True,
        help="Path to layer_counts.tsv.",
    )
    parser.add_argument(
        "--pairwise_overlaps_tsv",
        required=True,
        help="Path to pairwise_overlaps.tsv.",
    )
    parser.add_argument(
        "--combination_overlaps_tsv",
        required=True,
        help="Path to combination_overlaps.tsv.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for report and summary TSVs.",
    )
    return parser.parse_args()


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


def read_tsv(*, path: str) -> pd.DataFrame:
    """
    Read a TSV file.

    Parameters
    ----------
    path : str
        TSV path.

    Returns
    -------
    pandas.DataFrame
        Loaded DataFrame.
    """
    return pd.read_csv(path, sep="\t")


def compute_jaccard(*, df: pd.DataFrame) -> pd.DataFrame:
    """
    Add union and Jaccard index columns to a pairwise overlap table.

    Parameters
    ----------
    df : pandas.DataFrame
        Pairwise overlap DataFrame.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with Union_count and Jaccard columns.
    """
    out = df.copy()
    out["Union_count"] = out["A_count"] + out["B_count"] - out["Overlap_count"]
    out["Jaccard"] = (out["Overlap_count"] / out["Union_count"]).fillna(0.0)
    return out


def get_pair_row(
    *,
    pair_df: pd.DataFrame,
    a: str,
    b: str,
) -> Optional[pd.Series]:
    """
    Retrieve a pairwise overlap row regardless of ordering.

    Parameters
    ----------
    pair_df : pandas.DataFrame
        Pairwise overlap table.
    a : str
        Layer name A.
    b : str
        Layer name B.

    Returns
    -------
    pandas.Series or None
        Matching row, or None if not found.
    """
    mask = ((pair_df["Layer_A"] == a) & (pair_df["Layer_B"] == b)) | (
        (pair_df["Layer_A"] == b) & (pair_df["Layer_B"] == a)
    )
    if not mask.any():
        return None
    return pair_df.loc[mask].iloc[0]


def fmt_int(*, x: object) -> str:
    """
    Format an integer-like value with commas.

    Parameters
    ----------
    x : object
        Value to format.

    Returns
    -------
    str
        Formatted integer string.
    """
    try:
        return f"{int(x):,}"
    except Exception:
        return str(x)


def fmt_pct(*, x: object, ndp: int = 1) -> str:
    """
    Format a percentage value.

    Parameters
    ----------
    x : object
        Value to format.
    ndp : int
        Number of decimal places.

    Returns
    -------
    str
        Formatted percentage string.
    """
    try:
        return f"{float(x):.{ndp}f}%"
    except Exception:
        return str(x)


def fmt_float(*, x: object, ndp: int = 3) -> str:
    """
    Format a float.

    Parameters
    ----------
    x : object
        Value to format.
    ndp : int
        Number of decimal places.

    Returns
    -------
    str
        Formatted float string.
    """
    try:
        return f"{float(x):.{ndp}f}"
    except Exception:
        return str(x)


def build_key_numbers_table(*, pair_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a compact table of key overlaps with percentages and Jaccard.

    Parameters
    ----------
    pair_df : pandas.DataFrame
        Pairwise overlap table (with Jaccard computed).

    Returns
    -------
    pandas.DataFrame
        Key overlaps summary.
    """
    desired_pairs = [
        ("RNA_testis", "RNA_sperm"),
        ("RNA_sperm", "Prot_internal_supported_in_sperm_tsv"),
        ("RNA_testis", "Prot_internal_supported_in_sperm_tsv"),
        ("RNA_sperm", "Prot_internal_full"),
        ("RNA_testis", "Prot_internal_full"),
        ("RNA_sperm", "Prot_public"),
        ("RNA_testis", "Prot_public"),
        ("Prot_internal_full", "Prot_public"),
        ("Prot_internal_supported_in_sperm_tsv", "Prot_public"),
        ("Prot_internal_supported_in_sperm_tsv", "Prot_internal_full"),
    ]

    rows = []
    for a, b in desired_pairs:
        r = get_pair_row(pair_df=pair_df, a=a, b=b)
        if r is None:
            continue
        rows.append(
            {
                "Layer_A": a,
                "Layer_B": b,
                "A_count": int(r["A_count"]),
                "B_count": int(r["B_count"]),
                "Overlap_count": int(r["Overlap_count"]),
                "Overlap_pct_of_A": float(r["Overlap_pct_of_A"]),
                "Overlap_pct_of_B": float(r["Overlap_pct_of_B"]),
                "Jaccard": float(r["Jaccard"]),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = out.sort_values("Overlap_count", ascending=False).reset_index(drop=True)
    return out


def find_core_rows(*, combo_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pull the largest-combination rows (e.g. 4-way/5-way core sets).

    Parameters
    ----------
    combo_df : pandas.DataFrame
        Combination overlap table.

    Returns
    -------
    pandas.DataFrame
        Subset of combination rows for the maximum combination size.
    """
    if combo_df.empty or "Combination_size" not in combo_df.columns:
        return pd.DataFrame()

    max_size = int(combo_df["Combination_size"].max())
    core = combo_df.loc[combo_df["Combination_size"] == max_size].copy()
    if core.empty:
        return core
    core = core.sort_values("Overlap_count", ascending=False).reset_index(drop=True)
    return core


def write_markdown_report(
    *,
    out_path: str,
    layer_df: pd.DataFrame,
    key_df: pd.DataFrame,
    combo_df: pd.DataFrame,
) -> None:
    """
    Write a markdown report summarising layer sizes and key overlaps.

    Parameters
    ----------
    out_path : str
        Output markdown path.
    layer_df : pandas.DataFrame
        Layer counts table.
    key_df : pandas.DataFrame
        Key overlaps table.
    combo_df : pandas.DataFrame
        Combination overlaps table.

    Returns
    -------
    None
    """
    layer_df_sorted = layer_df.sort_values("Gene_count", ascending=False).copy()

    def md_table(df: pd.DataFrame) -> str:
        return df.to_markdown(index=False)

    # Extract headline RNA numbers if available
    rna_row = None
    if not key_df.empty:
        mask = (key_df["Layer_A"] == "RNA_testis") & (key_df["Layer_B"] == "RNA_sperm")
        if mask.any():
            rna_row = key_df.loc[mask].iloc[0]

    # Find a clean "core" 4-way set if present
    core_4_name = "RNA_testis__AND__RNA_sperm__AND__Prot_internal_full__AND__Prot_public"
    core_4 = None
    if "Combination" in combo_df.columns:
        m4 = combo_df["Combination"] == core_4_name
        if m4.any():
            core_4 = combo_df.loc[m4].iloc[0].to_dict()

    core_5 = None
    core_max = find_core_rows(combo_df=combo_df)
    if not core_max.empty:
        core_5 = core_max.iloc[0].to_dict()

    lines = []
    lines.append("# Multi-omic overlap summary")
    lines.append("")
    lines.append("## Layer sizes")
    lines.append("")
    lines.append(md_table(layer_df_sorted[["Layer", "Gene_count"]]))
    lines.append("")

    lines.append("## Key pairwise overlaps")
    lines.append("")
    if key_df.empty:
        lines.append("No key overlap rows were found in the pairwise table.")
    else:
        key_show = key_df.copy()
        key_show["A_count"] = key_show["A_count"].map(lambda x: fmt_int(x=x))
        key_show["B_count"] = key_show["B_count"].map(lambda x: fmt_int(x=x))
        key_show["Overlap_count"] = key_show["Overlap_count"].map(lambda x: fmt_int(x=x))
        key_show["Overlap_pct_of_A"] = key_show["Overlap_pct_of_A"].map(lambda x: fmt_pct(x=x))
        key_show["Overlap_pct_of_B"] = key_show["Overlap_pct_of_B"].map(lambda x: fmt_pct(x=x))
        key_show["Jaccard"] = key_show["Jaccard"].map(lambda x: fmt_float(x=x))
        lines.append(md_table(key_show))
    lines.append("")

    lines.append("## Results-ready statements")
    lines.append("")
    if rna_row is not None:
        lines.append(
            (
                f"- Testis-specific genes (RNA_testis): {fmt_int(x=rna_row['A_count'])}\n"
                f"- Sperm-present genes (RNA_sperm): {fmt_int(x=rna_row['B_count'])}\n"
                f"- Overlap (RNA_testis ∩ RNA_sperm): {fmt_int(x=rna_row['Overlap_count'])} "
                f"({fmt_pct(x=rna_row['Overlap_pct_of_A'])} of RNA_testis; "
                f"{fmt_pct(x=rna_row['Overlap_pct_of_B'])} of RNA_sperm)"
            )
        )
    else:
        lines.append("- RNA overlap statement not available (RNA_testis/RNA_sperm row missing).")

    if core_4 is not None:
        lines.append("")
        lines.append(
            (
                f"- Core 4-way overlap (RNA_testis ∩ RNA_sperm ∩ Prot_internal_full ∩ Prot_public): "
                f"{fmt_int(x=core_4['Overlap_count'])} genes "
                f"({fmt_pct(x=core_4.get('Overlap_pct_of_RNA_testis', 'NA'))} of RNA_testis; "
                f"{fmt_pct(x=core_4.get('Overlap_pct_of_RNA_sperm', 'NA'))} of RNA_sperm; "
                f"{fmt_pct(x=core_4.get('Overlap_pct_of_Prot_internal_full', 'NA'))} of Prot_internal_full; "
                f"{fmt_pct(x=core_4.get('Overlap_pct_of_Prot_public', 'NA'))} of Prot_public)"
            )
        )

    if core_5 is not None and core_5.get("Combination_size", 0) >= 5:
        lines.append("")
        lines.append(
            (
                f"- Maximum-overlap set (size {int(core_5['Combination_size'])}): "
                f"{fmt_int(x=core_5['Overlap_count'])} genes "
                f"(Combination: {core_5['Combination']})"
            )
        )

    lines.append("")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def main() -> int:
    """
    Run summary generation.

    Returns
    -------
    int
        Exit code.
    """
    args = parse_args()
    safe_makedirs(path=args.out_dir)

    layer_df = read_tsv(path=args.layer_counts_tsv)
    pair_df = read_tsv(path=args.pairwise_overlaps_tsv)
    combo_df = read_tsv(path=args.combination_overlaps_tsv)

    # Compute Jaccard for pairwise table
    pair_df = compute_jaccard(df=pair_df)

    key_df = build_key_numbers_table(pair_df=pair_df)
    core_df = find_core_rows(combo_df=combo_df)

    key_out = os.path.join(args.out_dir, "overlap_summary_key_numbers.tsv")
    core_out = os.path.join(args.out_dir, "overlap_core_sets.tsv")
    md_out = os.path.join(args.out_dir, "overlap_summary_report.md")

    key_df.to_csv(key_out, sep="\t", index=False)
    core_df.to_csv(core_out, sep="\t", index=False)

    write_markdown_report(
        out_path=md_out,
        layer_df=layer_df,
        key_df=key_df,
        combo_df=combo_df,
    )

    print(f"Wrote: {key_out}")
    print(f"Wrote: {core_out}")
    print(f"Wrote: {md_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
