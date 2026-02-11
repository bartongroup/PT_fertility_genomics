#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fertility_gene_set_plots.py

Generate summary visualisations (UpSet plots, Venn diagrams, bar charts, and
pairwise Jaccard heatmap) from the master Excel workbook.

This script is designed to read the tier/flag summary sheet produced by the
master workbook builder, and then visualise the overlap between evidence sets.

Outputs
-------
- UpSet plots as PNG:
  - upset_all_flags.png
  - upset_core_flags.png
- Venn diagrams as PNG (3-way, where possible):
  - venn_hpo_clinvarhcpath_testishc.png
  - venn_testishc_proteomics_sperm.png
- Bar chart of set sizes:
  - set_sizes.png
- Pairwise Jaccard heatmap:
  - jaccard_heatmap.png
- TSV summaries:
  - set_sizes.tsv
  - intersection_counts_top.tsv

All arguments are named. The input workbook is assumed to contain a sheet with
a 'gene_key' column and multiple boolean-like flag columns.

Notes
-----
- Venn diagrams are only sensible up to 3 sets.
- UpSet plots can handle many sets and is preferred for the full overview.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd
import matplotlib.pyplot as plt

try:
    from upsetplot import UpSet, from_indicators
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "Missing dependency 'upsetplot'. Install with: python -m pip install upsetplot"
    ) from exc

try:
    from matplotlib_venn import venn3
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "Missing dependency 'matplotlib-venn'. Install with: python -m pip install matplotlib-venn"
    ) from exc


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create UpSet/Venn/summary plots from the master Excel workbook."
    )
    parser.add_argument(
        "--in_xlsx",
        required=True,
        type=Path,
        help="Input Excel workbook (.xlsx), e.g. results/master_workbook/master_fertility_genes.xlsx",
    )
    parser.add_argument(
        "--sheet_name",
        required=False,
        default="Tier_Summary_With_Omics",
        type=str,
        help="Sheet to read from the workbook (default: Tier_Summary_With_Omics).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        type=Path,
        help="Output directory for plots and summaries.",
    )
    parser.add_argument(
        "--gene_col",
        required=False,
        default="gene_key",
        type=str,
        help="Gene identifier column (default: gene_key).",
    )
    parser.add_argument(
        "--max_intersections",
        required=False,
        default=40,
        type=int,
        help="Max number of intersections to display in UpSet (default: 40).",
    )
    return parser.parse_args()


def series_to_bool(series: pd.Series) -> pd.Series:
    """
    Convert a pandas Series of mixed values into a boolean Series.

    Parameters
    ----------
    series : pd.Series
        Input series containing True/False, 0/1, yes/no, strings, or blanks.

    Returns
    -------
    pd.Series
        Boolean series.
    """
    s = series.fillna("").astype(str).str.strip().str.lower()
    return s.isin({"true", "t", "1", "yes", "y"})


def load_flags_table(
    in_xlsx: Path,
    sheet_name: str,
    gene_col: str,
) -> pd.DataFrame:
    """
    Load the tier/flag sheet and coerce flag columns to booleans.

    Parameters
    ----------
    in_xlsx : Path
        Input Excel workbook.
    sheet_name : str
        Sheet name to read.
    gene_col : str
        Column containing gene identifiers.

    Returns
    -------
    pd.DataFrame
        DataFrame with gene_col plus boolean flag columns.
    """
    df = pd.read_excel(in_xlsx, sheet_name=sheet_name, dtype=str)
    if gene_col not in df.columns:
        raise ValueError(
            f"Expected a '{gene_col}' column in sheet '{sheet_name}'. "
            f"Columns: {list(df.columns)}"
        )

    df[gene_col] = df[gene_col].fillna("").astype(str).str.strip()
    df = df[df[gene_col] != ""].copy()

    # Identify candidate flag columns (everything except gene_col).
    flag_cols = [c for c in df.columns if c != gene_col]

    # Coerce to booleans.
    for c in flag_cols:
        df[c] = series_to_bool(df[c])

    # Collapse duplicates if they exist (any True wins).
    df = df.groupby(gene_col, as_index=False)[flag_cols].any()

    return df


def write_set_sizes(df: pd.DataFrame, gene_col: str, out_path: Path) -> pd.DataFrame:
    """
    Write per-flag set sizes to TSV.

    Parameters
    ----------
    df : pd.DataFrame
        Flags table.
    gene_col : str
        Gene identifier column.
    out_path : Path
        Output TSV path.

    Returns
    -------
    pd.DataFrame
        Set size summary DataFrame.
    """
    flag_cols = [c for c in df.columns if c != gene_col]
    rows = []
    for c in flag_cols:
        rows.append({"set": c, "n_genes": int(df[c].sum())})
    out = pd.DataFrame(rows).sort_values("n_genes", ascending=False)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def plot_set_sizes(set_sizes: pd.DataFrame, out_png: Path) -> None:
    """
    Plot bar chart of set sizes.

    Parameters
    ----------
    set_sizes : pd.DataFrame
        DataFrame with columns: set, n_genes.
    out_png : Path
        Output PNG path.
    """
    plt.figure()
    plt.bar(set_sizes["set"], set_sizes["n_genes"])
    plt.xticks(rotation=60, ha="right")
    plt.ylabel("Number of genes")
    plt.title("Evidence set sizes")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def upset_plot(
    df: pd.DataFrame,
    gene_col: str,
    flag_cols: Sequence[str],
    out_png: Path,
    title: str,
    max_intersections: int,
) -> pd.DataFrame:
    """
    Generate an UpSet plot and return the intersection size table.

    Parameters
    ----------
    df : pd.DataFrame
        Flags table.
    gene_col : str
        Gene identifier column.
    flag_cols : Sequence[str]
        Columns to include in the UpSet plot.
    out_png : Path
        Output PNG path.
    title : str
        Plot title.
    max_intersections : int
        Maximum intersections to display.

    Returns
    -------
    pd.DataFrame
        Intersection counts as a DataFrame.
    """
    work = df[[gene_col] + list(flag_cols)].copy()
    indicators = work[list(flag_cols)].astype(bool)

    upset_data = from_indicators(indicators=indicators, data=work[gene_col])

    plt.figure()
    upset = UpSet(
        upset_data,
        show_counts=True,
        sort_by="cardinality",
        intersection_plot_elements=max_intersections,
    )
    upset.plot()
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    # Build intersection count table
    counts = upset_data.groupby(level=list(range(len(flag_cols)))).size()
    counts = counts.sort_values(ascending=False).head(max_intersections)
    out_rows: List[Dict[str, object]] = []
    for idx, n in counts.items():
        # idx is a tuple of booleans in the same order as flag_cols
        included = [flag_cols[i] for i, v in enumerate(idx) if bool(v)]
        out_rows.append({"intersection": "&".join(included) if included else "(none)", "n_genes": int(n)})

    return pd.DataFrame(out_rows)


def jaccard_similarity_matrix(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    """
    Compute pairwise Jaccard similarity between all boolean flag columns.

    Parameters
    ----------
    df : pd.DataFrame
        Flags table.
    gene_col : str
        Gene identifier column.

    Returns
    -------
    pd.DataFrame
        Square DataFrame of Jaccard similarities.
    """
    flag_cols = [c for c in df.columns if c != gene_col]
    mat = pd.DataFrame(index=flag_cols, columns=flag_cols, dtype=float)

    for a in flag_cols:
        A = df[a].astype(bool)
        for b in flag_cols:
            B = df[b].astype(bool)
            inter = int((A & B).sum())
            union = int((A | B).sum())
            mat.loc[a, b] = (inter / union) if union > 0 else 0.0
    return mat


def plot_jaccard_heatmap(mat: pd.DataFrame, out_png: Path) -> None:
    """
    Plot a pairwise Jaccard similarity heatmap.

    Parameters
    ----------
    mat : pd.DataFrame
        Jaccard similarity matrix.
    out_png : Path
        Output PNG path.
    """
    plt.figure()
    plt.imshow(mat.values, aspect="auto")
    plt.colorbar(label="Jaccard similarity")
    plt.xticks(range(len(mat.columns)), mat.columns, rotation=60, ha="right")
    plt.yticks(range(len(mat.index)), mat.index)
    plt.title("Pairwise overlap (Jaccard)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_venn3_from_flags(
    df: pd.DataFrame,
    gene_col: str,
    a: str,
    b: str,
    c: str,
    out_png: Path,
    title: str,
) -> None:
    """
    Plot a 3-way Venn diagram for three boolean columns.

    Parameters
    ----------
    df : pd.DataFrame
        Flags table.
    gene_col : str
        Gene identifier column.
    a : str
        First set column.
    b : str
        Second set column.
    c : str
        Third set column.
    out_png : Path
        Output PNG path.
    title : str
        Plot title.
    """
    for col in (a, b, c):
        if col not in df.columns:
            raise ValueError(f"Missing column for Venn: {col}")

    A = set(df.loc[df[a], gene_col])
    B = set(df.loc[df[b], gene_col])
    C = set(df.loc[df[c], gene_col])

    plt.figure()
    venn3([A, B, C], set_labels=(a, b, c))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = load_flags_table(
        in_xlsx=args.in_xlsx,
        sheet_name=args.sheet_name,
        gene_col=args.gene_col,
    )

    # Set sizes
    set_sizes = write_set_sizes(
        df=df,
        gene_col=args.gene_col,
        out_path=args.out_dir / "set_sizes.tsv",
    )
    plot_set_sizes(set_sizes=set_sizes, out_png=args.out_dir / "set_sizes.png")

    # UpSet: all flags
    all_flag_cols = [c for c in df.columns if c != args.gene_col]
    intersections_all = upset_plot(
        df=df,
        gene_col=args.gene_col,
        flag_cols=all_flag_cols,
        out_png=args.out_dir / "upset_all_flags.png",
        title="UpSet: all evidence flags",
        max_intersections=args.max_intersections,
    )
    intersections_all.to_csv(args.out_dir / "intersection_counts_top.tsv", sep="\t", index=False)

    # UpSet: core subset (more presentation-friendly)
    core_cols = [
        "in_hpo_gene_set",
        "clinvar_hc_pathogenic_present",
        "in_testis_high_conf_final",
        "proteomics_present",
        "sperm_rnaseq_present",
    ]
    core_cols = [c for c in core_cols if c in df.columns]
    if len(core_cols) >= 2:
        _ = upset_plot(
            df=df,
            gene_col=args.gene_col,
            flag_cols=core_cols,
            out_png=args.out_dir / "upset_core_flags.png",
            title="UpSet: core evidence flags",
            max_intersections=args.max_intersections,
        )

    # Venn diagrams (only if all three columns exist)
    if all(c in df.columns for c in ["in_hpo_gene_set", "clinvar_hc_pathogenic_present", "in_testis_high_conf_final"]):
        plot_venn3_from_flags(
            df=df,
            gene_col=args.gene_col,
            a="in_hpo_gene_set",
            b="clinvar_hc_pathogenic_present",
            c="in_testis_high_conf_final",
            out_png=args.out_dir / "venn_hpo_clinvarhcpath_testishc.png",
            title="Venn: HPO vs ClinVar HC pathogenic vs Testis HC final",
        )

    if all(c in df.columns for c in ["in_testis_high_conf_final", "proteomics_present", "sperm_rnaseq_present"]):
        plot_venn3_from_flags(
            df=df,
            gene_col=args.gene_col,
            a="in_testis_high_conf_final",
            b="proteomics_present",
            c="sperm_rnaseq_present",
            out_png=args.out_dir / "venn_testishc_proteomics_sperm.png",
            title="Venn: Testis HC final vs Proteomics vs Sperm RNA-seq",
        )

    # Pairwise overlap heatmap
    jac = jaccard_similarity_matrix(df=df, gene_col=args.gene_col)
    plot_jaccard_heatmap(mat=jac, out_png=args.out_dir / "jaccard_heatmap.png")
    jac.to_csv(args.out_dir / "jaccard_matrix.tsv", sep="\t")

    # Quick log to stdout (handy in cluster logs)
    print(f"Loaded genes: {df.shape[0]}")
    print(f"Flags: {len([c for c in df.columns if c != args.gene_col])}")
    print(f"Wrote outputs to: {args.out_dir}")


if __name__ == "__main__":
    main()
