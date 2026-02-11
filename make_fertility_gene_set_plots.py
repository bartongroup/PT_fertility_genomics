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
- UpSet plots as pdf:
  - upset_all_flags.pdf
  - upset_core_flags.pdf
- Venn diagrams as pdf (3-way, where possible):
  - venn_hpo_clinvarhcpath_testishc.pdf
  - venn_testishc_proteomics_sperm.pdf
- Bar chart of set sizes:
  - set_sizes.pdf
- Pairwise Jaccard heatmap:
  - jaccard_heatmap.pdf
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
import matplotlib.p
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import re
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

import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


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


def plot_set_sizes(set_sizes: pd.DataFrame, out_pdf: Path) -> None:
    """
    Plot bar chart of set sizes.

    Parameters
    ----------
    set_sizes : pd.DataFrame
        DataFrame with columns: set, n_genes.
    out_pdf : Path
        Output pdf path.
    """
    plt.figure()
    plt.bar(set_sizes["set"], set_sizes["n_genes"])
    plt.xticks(rotation=60, ha="right")
    plt.ylabel("Number of genes")
    plt.title("Evidence set sizes")
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=200)
    plt.close()



def upset_plot(
    df: pd.DataFrame,
    gene_col: str,
    flag_cols: Sequence[str],
    out_pdf: Path,
    title: str,
    max_intersections: int,
) -> pd.DataFrame:
    """
    Generate an UpSet plot and return the intersection size table.

    This implementation uses subset_size="count" to safely handle non-unique
    membership patterns (many genes share the same combination of flags).

    Parameters
    ----------
    df : pd.DataFrame
        Flags table (one row per gene, boolean flag columns).
    gene_col : str
        Gene identifier column.
    flag_cols : Sequence[str]
        Columns to include in the UpSet plot.
    out_pdf : Path
        Output pdf path.
    title : str
        Plot title.
    max_intersections : int
        Maximum intersections to display.

    Returns
    -------
    pd.DataFrame
        Intersection counts table (top intersections).
    """
    work = df[[gene_col] + list(flag_cols)].copy()

    # Ensure booleans (defensive)
    for c in flag_cols:
        work[c] = work[c].fillna(False).astype(bool)

    # Build a Series of 1s indexed by the membership pattern.
    # Duplicated membership patterns are fine when subset_size="count".
    upset_series = (
        work.assign(_count=1)
        .set_index(list(flag_cols))["_count"]
    )

    plt.figure()
    upset = UpSet(
        upset_series,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        intersection_plot_elements=max_intersections,
    )
    upset.plot()
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=200)
    plt.close()

    # Intersection counts table (top N)
    counts = upset_series.groupby(level=list(range(len(flag_cols)))).sum()
    counts = counts.sort_values(ascending=False).head(max_intersections)

    out_rows: List[Dict[str, object]] = []
    for idx, n in counts.items():
        included = [flag_cols[i] for i, v in enumerate(idx) if bool(v)]
        out_rows.append(
            {
                "intersection": "&".join(included) if included else "(none)",
                "n_genes": int(n),
            }
        )

    return pd.DataFrame(out_rows)



def _series_to_bool(series: pd.Series) -> pd.Series:
    """
    Convert a pandas Series of mixed values (True/False, 0/1, yes/no, strings)
    into a boolean Series.

    Parameters
    ----------
    series : pd.Series
        Input Series.

    Returns
    -------
    pd.Series
        Boolean Series.
    """
    s = series.fillna("").astype(str).str.strip().str.lower()
    return s.isin({"true", "t", "1", "yes", "y"})


def _safe_float(series: pd.Series) -> pd.Series:
    """
    Coerce a Series to float safely (invalid parses become NaN).

    Parameters
    ----------
    series : pd.Series
        Input Series.

    Returns
    -------
    pd.Series
        Float Series.
    """
    return pd.to_numeric(series, errors="coerce")


def make_evidence_score_table(
    tier_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    gene_col: str,
    out_tsv: Path,
    top_n: int = 200,
) -> pd.DataFrame:
    """
    Create a ranked 'top genes' table using a simple evidence score.

    Scoring (default, simple and explainable)
    ----------------------------------------
    - in_hpo_gene_set: 1
    - sperm_rnaseq_present: 1
    - proteomics_present: 1
    - clinvar_best_present: 1
    - clinvar_best_pathogenic_present: 2
    - clinvar_hc_present: 2
    - clinvar_hc_pathogenic_present: 3
    - in_testis_high_conf_final: 2

    Parameters
    ----------
    tier_df : pd.DataFrame
        Tier sheet with boolean flags.
    genes_df : pd.DataFrame
        Genes_Master sheet with tau/testis TPM etc.
    gene_col : str
        Column name for gene key in both tables (usually 'gene_key').
    out_tsv : Path
        Output TSV path.
    top_n : int, optional
        Keep top N genes in the output table, by default 200.

    Returns
    -------
    pd.DataFrame
        Ranked evidence table.
    """
    flags = [
        "in_hpo_gene_set",
        "sperm_rnaseq_present",
        "proteomics_present",
        "clinvar_best_present",
        "clinvar_best_pathogenic_present",
        "clinvar_hc_present",
        "clinvar_hc_pathogenic_present",
        "in_testis_high_conf_final",
    ]
    weights: Dict[str, int] = {
        "in_hpo_gene_set": 1,
        "sperm_rnaseq_present": 1,
        "proteomics_present": 1,
        "clinvar_best_present": 1,
        "clinvar_best_pathogenic_present": 2,
        "clinvar_hc_present": 2,
        "clinvar_hc_pathogenic_present": 3,
        "in_testis_high_conf_final": 2,
    }

    work = tier_df[[gene_col] + [c for c in flags if c in tier_df.columns]].copy()
    work[gene_col] = work[gene_col].fillna("").astype(str).str.strip()
    work = work[work[gene_col] != ""].drop_duplicates(subset=[gene_col])

    for c in flags:
        if c not in work.columns:
            work[c] = False
        work[c] = _series_to_bool(work[c])

    score = np.zeros(len(work), dtype=int)
    for c in flags:
        score += work[c].astype(int).values * int(weights[c])
    work["evidence_score"] = score

    # Bring in tau/testis TPM for convenience in the ranked table
    cols_from_genes = [gene_col]
    for c in ["tau", "target_median_tpm", "sperm_tpm_median", "sperm_tpm_mean"]:
        if c in genes_df.columns:
            cols_from_genes.append(c)

    merged = work.merge(genes_df[cols_from_genes].drop_duplicates(subset=[gene_col]),
                        on=gene_col, how="left")

    if "tau" in merged.columns:
        merged["tau"] = _safe_float(merged["tau"])
    if "target_median_tpm" in merged.columns:
        merged["target_median_tpm"] = _safe_float(merged["target_median_tpm"])

    merged = merged.sort_values(
        by=["evidence_score", "tau", "target_median_tpm"],
        ascending=[False, False, False],
        na_position="last",
    )

    merged_top = merged.head(int(top_n)).copy()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    merged_top.to_csv(out_tsv, sep="\t", index=False)

    return merged_top


def assign_clinvar_tier_label(df: pd.DataFrame) -> pd.Series:
    """
    Assign a single ClinVar tier label per gene from boolean flags.

    Priority (highest to lowest)
    ----------------------------
    1) HC pathogenic
    2) HC any
    3) Best pathogenic
    4) Best any
    5) None

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing ClinVar boolean flags.

    Returns
    -------
    pd.Series
        Tier label per row.
    """
    for col in [
        "clinvar_hc_pathogenic_present",
        "clinvar_hc_present",
        "clinvar_best_pathogenic_present",
        "clinvar_best_present",
    ]:
        if col not in df.columns:
            df[col] = False
        df[col] = _series_to_bool(df[col])

    labels: List[str] = []
    for _, row in df.iterrows():
        if bool(row["clinvar_hc_pathogenic_present"]):
            labels.append("ClinVar HC pathogenic")
        elif bool(row["clinvar_hc_present"]):
            labels.append("ClinVar HC")
        elif bool(row["clinvar_best_pathogenic_present"]):
            labels.append("ClinVar best pathogenic")
        elif bool(row["clinvar_best_present"]):
            labels.append("ClinVar best")
        else:
            labels.append("No ClinVar tier")
    return pd.Series(labels, index=df.index)


def plot_tau_vs_testis_tpm(
    genes_df: pd.DataFrame,
    gene_col: str,
    out_pdf: Path,
    title: str,
) -> None:
    """
    Scatter plot of tau vs testis median TPM, coloured by ClinVar tier.

    Parameters
    ----------
    genes_df : pd.DataFrame
        Genes_Master table.
    gene_col : str
        Gene key column.
    out_pdf : Path
        Output pdf path.
    title : str
        Plot title.
    """
    required = {"tau", "target_median_tpm", gene_col}
    missing = [c for c in required if c not in genes_df.columns]
    if missing:
        raise ValueError(f"Genes_Master missing required columns: {missing}")

    work = genes_df[[gene_col, "tau", "target_median_tpm",
                     "clinvar_best_present", "clinvar_best_pathogenic_present",
                     "clinvar_hc_present", "clinvar_hc_pathogenic_present"]].copy()

    work[gene_col] = work[gene_col].fillna("").astype(str).str.strip()
    work = work[work[gene_col] != ""].drop_duplicates(subset=[gene_col])

    work["tau"] = _safe_float(work["tau"])
    work["target_median_tpm"] = _safe_float(work["target_median_tpm"])
    work = work.dropna(subset=["tau", "target_median_tpm"])

    work["clinvar_tier"] = assign_clinvar_tier_label(work)

    # Plot
    plt.figure()
    categories = [
        "ClinVar HC pathogenic",
        "ClinVar HC",
        "ClinVar best pathogenic",
        "ClinVar best",
        "No ClinVar tier",
    ]

    for cat in categories:
        sub = work[work["clinvar_tier"] == cat]
        if sub.empty:
            continue
        plt.scatter(sub["tau"], sub["target_median_tpm"], label=cat, alpha=0.7)

    plt.xlabel("Tau (testis specificity)")
    plt.ylabel("Testis median TPM")
    plt.yscale("log")
    plt.title(title)
    plt.legend(loc="best", frameon=False)
    plt.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, dpi=200)
    plt.close()


def plot_tissue_heatmap_for_shortlist(
    tissue_df: pd.DataFrame,
    shortlist_gene_ids: List[str],
    id_col: str,
    out_pdf: Path,
    title: str,
    max_genes: int = 50,
    top_tissues: int = 15,
) -> None:
    """
    Plot a heatmap of GTEx tissue median TPM for a shortlist of genes.

    Parameters
    ----------
    tissue_df : pd.DataFrame
        GTEx tissue median table from the workbook.
    shortlist_gene_ids : List[str]
        List of gene identifiers to match (e.g. Ensembl gene IDs).
    id_col : str
        Column in tissue_df containing those identifiers (commonly 'Name').
    out_pdf : Path
        Output pdf path.
    title : str
        Plot title.
    max_genes : int
        Max genes to plot.
    top_tissues : int
        Max tissues to plot (highest variance across shortlist).
    """
    if id_col not in tissue_df.columns:
        raise ValueError(f"GTEx_Tissue_Medians is missing expected id column: {id_col}")

    df = tissue_df.copy()
    df[id_col] = df[id_col].fillna("").astype(str).str.strip()
    df["_id_norm"] = df["Name"].str.replace(r"\.\d+$", "", regex=True)

    shortlist_norm = (
        pd.Series(shortlist_gene_ids, dtype="string")
        .fillna("")
        .astype(str)
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
        .tolist()
    )

    df_sub = df[df["_id_norm"].isin(set(shortlist_norm))].copy()
    if df_sub.empty:
        print(f"Shortlist size: {len(shortlist_norm)}")
        print(f"Tissue IDs size: {df['_id_norm'].nunique()}")
        print("Example shortlist IDs:", shortlist_norm[:5])
        print("Example tissue IDs:", df["_id_norm"].head(5).tolist())
        raise ValueError("No shortlist genes matched rows in GTEx_Tissue_Medians.")

    # Keep only numeric tissue columns (everything except identifiers/description)
    drop_cols = {id_col, "_id_norm", "Description"}
    tissue_cols = [c for c in df_sub.columns if c not in drop_cols]

    # Convert to numeric TPM
    for c in tissue_cols:
        df_sub[c] = pd.to_numeric(df_sub[c], errors="coerce").fillna(0.0)

    # Pick top tissues by variance across shortlist
    variances = df_sub[tissue_cols].var(axis=0).sort_values(ascending=False)
    keep_tissues = variances.head(top_tissues).index.tolist()

    # Limit genes
    df_sub = df_sub.head(max_genes)

    mat = df_sub[keep_tissues].to_numpy()
    plt.figure(figsize=(max(10, len(keep_tissues) * 0.6), max(8, df_sub.shape[0] * 0.25)))
    plt.imshow(mat, aspect="auto")
    plt.yticks(range(df_sub.shape[0]), df_sub[id_col].tolist())
    plt.xticks(range(len(keep_tissues)), keep_tissues, rotation=90)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=200)
    plt.close()


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


def plot_jaccard_heatmap(mat: pd.DataFrame, out_pdf: Path) -> None:
    """
    Plot a pairwise Jaccard similarity heatmap.

    Parameters
    ----------
    mat : pd.DataFrame
        Jaccard similarity matrix.
    out_pdf : Path
        Output pdf path.
    """
    plt.figure()
    plt.imshow(mat.values, aspect="auto")
    plt.colorbar(label="Jaccard similarity")
    plt.xticks(range(len(mat.columns)), mat.columns, rotation=60, ha="right")
    plt.yticks(range(len(mat.index)), mat.index)
    plt.title("Pairwise overlap (Jaccard)")
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=200)
    plt.close()


def plot_venn3_from_flags(
    df: pd.DataFrame,
    gene_col: str,
    a: str,
    b: str,
    c: str,
    out_pdf: Path,
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
    out_pdf : Path
        Output pdf path.
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
    plt.savefig(out_pdf, dpi=200)
    plt.close()



def main() -> None:
    """
    Entry point.

    This main function:
    1) Loads the Tier_Summary flags table (boolean membership per gene).
    2) Produces set size summaries, UpSet plots, Venn diagrams, and Jaccard heatmap.
    3) Additionally reads Genes_Master (and optionally GTEx_Tissue_Medians) to make:
       - top genes table ranked by a simple evidence score
       - tau vs testis TPM scatter coloured by ClinVar tier
       - tissue median TPM heatmap for the top-ranked shortlist (if available)
    """
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1) Load gene-set flags (Tier summary sheet)
    # ------------------------------------------------------------------
    df = load_flags_table(
        in_xlsx=args.in_xlsx,
        sheet_name=args.sheet_name,
        gene_col=args.gene_col,
    )


    genes_master = pd.read_excel(
    io=args.in_xlsx,
    sheet_name="Genes_Master",
    dtype=str,
    )

    genes_master["gene_key"] = genes_master["gene_key"].fillna("").astype(str).str.strip()
    genes_master["Name"] = genes_master["Name"].fillna("").astype(str).str.strip()

    shortlist_symbols = set(
        df.loc[df["in_testis_high_conf_final"].astype(bool), args.gene_col]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    shortlist_ensembl_versioned = (
        genes_master.loc[genes_master["gene_key"].isin(shortlist_symbols), "Name"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )


    # ------------------------------------------------------------------
    # 2) Set sizes + set size plot
    # ------------------------------------------------------------------
    set_sizes = write_set_sizes(
        df=df,
        gene_col=args.gene_col,
        out_path=args.out_dir / "set_sizes.tsv",
    )
    plot_set_sizes(set_sizes=set_sizes, out_pdf=args.out_dir / "set_sizes.pdf")

    # ------------------------------------------------------------------
    # 3) UpSet plots
    # ------------------------------------------------------------------
    all_flag_cols = [c for c in df.columns if c != args.gene_col]
    intersections_all = upset_plot(
        df=df,
        gene_col=args.gene_col,
        flag_cols=all_flag_cols,
        out_pdf=args.out_dir / "upset_all_flags.pdf",
        title="UpSet: all evidence flags",
        max_intersections=args.max_intersections,
    )
    intersections_all.to_csv(
        args.out_dir / "intersection_counts_top.tsv",
        sep="\t",
        index=False,
    )

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
            out_pdf=args.out_dir / "upset_core_flags.pdf",
            title="UpSet: core evidence flags",
            max_intersections=args.max_intersections,
        )

    # ------------------------------------------------------------------
    # 4) Venn diagrams
    # ------------------------------------------------------------------
    if all(
        c in df.columns
        for c in [
            "in_hpo_gene_set",
            "clinvar_hc_pathogenic_present",
            "in_testis_high_conf_final",
        ]
    ):
        plot_venn3_from_flags(
            df=df,
            gene_col=args.gene_col,
            a="in_hpo_gene_set",
            b="clinvar_hc_pathogenic_present",
            c="in_testis_high_conf_final",
            out_pdf=args.out_dir / "venn_hpo_clinvarhcpath_testishc.pdf",
            title="Venn: HPO vs ClinVar HC pathogenic vs Testis HC final",
        )

    if all(
        c in df.columns
        for c in [
            "in_testis_high_conf_final",
            "proteomics_present",
            "sperm_rnaseq_present",
        ]
    ):
        plot_venn3_from_flags(
            df=df,
            gene_col=args.gene_col,
            a="in_testis_high_conf_final",
            b="proteomics_present",
            c="sperm_rnaseq_present",
            out_pdf=args.out_dir / "venn_testishc_proteomics_sperm.pdf",
            title="Venn: Testis HC final vs Proteomics vs Sperm RNA-seq",
        )

    # ------------------------------------------------------------------
    # 5) Pairwise overlap heatmap
    # ------------------------------------------------------------------
    jac = jaccard_similarity_matrix(df=df, gene_col=args.gene_col)
    plot_jaccard_heatmap(mat=jac, out_pdf=args.out_dir / "jaccard_heatmap.pdf")
    jac.to_csv(args.out_dir / "jaccard_matrix.tsv", sep="\t")

    # ------------------------------------------------------------------
    # 6) Extra visualisations from Genes_Master / GTEx_Tissue_Medians
    # ------------------------------------------------------------------

    ranked_tsv = args.out_dir / "top_genes_by_evidence_score.tsv"
    ranked = make_evidence_score_table(
        tier_df=df,
        genes_df=genes_master,
        gene_col=args.gene_col,
        out_tsv=ranked_tsv,
        top_n=200,
    )

    plot_tau_vs_testis_tpm(
        genes_df=genes_master,
        gene_col=args.gene_col,
        out_pdf=args.out_dir / "tau_vs_testis_median_tpm_by_clinvar_tier.pdf",
        title="Tau vs testis median TPM (colour by ClinVar tier)",
    )

    try:
        tissue_medians = pd.read_excel(
            args.in_xlsx,
            sheet_name="GTEx_Tissue_Medians",
            dtype=str,
        )
    except ValueError:
        tissue_medians = None

    if tissue_medians is not None:
        shortlist_symbols = ranked[args.gene_col].head(50).astype(str).tolist()

        top_symbols = ranked[args.gene_col].head(50).astype(str).str.strip().tolist()

        shortlist_ensembl_versioned = (
            genes_master.loc[genes_master[args.gene_col].isin(top_symbols), "Name"]
            .dropna()
            .astype(str)
            .str.strip()
            .tolist()
        )


        gm = genes_master.copy()
        if args.gene_col not in gm.columns and "gene_key" in gm.columns:
            gm[args.gene_col] = gm["gene_key"]

        gm[args.gene_col] = gm[args.gene_col].fillna("").astype(str).str.strip()
        gm["ensembl_gene_id"] = gm.get("ensembl_gene_id", "").fillna("").astype(str).str.strip()

        shortlist_map = gm.loc[
            gm[args.gene_col].isin(shortlist_symbols),
            ["ensembl_gene_id", args.gene_col],
        ].drop_duplicates()

        shortlist_ensembl = shortlist_map["ensembl_gene_id"].tolist()

        shortlist_ensembl = [
                            re.sub(r"\.\d+$", "", str(x).strip())
                            for x in shortlist_ensembl
                            if str(x).strip() != ""
                        ]


        plot_tissue_heatmap_for_shortlist(
            tissue_df=tissue_medians,
            shortlist_gene_ids=shortlist_ensembl_versioned,
            id_col="Name",
            out_pdf=args.out_dir / "heatmap_tissue_medians_top50.pdf",
            title="GTEx tissue median TPM (top 50 by evidence score)",
            max_genes=50,
            top_tissues=15,
        )

        shortlist_map.to_csv(
            args.out_dir / "top50_symbol_to_ensembl.tsv",
            sep="\t",
            index=False,
        )


    # ------------------------------------------------------------------
    # 7) log to stdout (handy in cluster logs)
    # ------------------------------------------------------------------
    print(f"Loaded genes: {df.shape[0]}")
    print(f"Flags: {len([c for c in df.columns if c != args.gene_col])}")
    print(f"Wrote outputs to: {args.out_dir}")


if __name__ == "__main__":
    main()



