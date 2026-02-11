#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_master_fertility_workbook.py

Create a single Excel workbook that consolidates:
- GTEx testis-specificity metrics (tau, TPMs, tissue medians pointers)
- Sperm RNA-seq presence/TPM summary
- Functional annotation (e.g. MyGene fields if present)
- HPO infertility gene associations (gene-level summary)
- ClinVar evidence tiers (best, high-confidence, high-confidence pathogenic)
- Proteomics presence and summary fields (if present)

Inputs are expected to live under the user's existing layout:
  ~/data/2026_sperm_Gates/results/<testis_run_id>/
  ~/data/2026_sperm_Gates/results/<male_infertility_run_id>/   (or similar)

The script writes a multi-sheet .xlsx file, with a gene-centric "Genes_Master"
sheet plus optional variant-level sheets.

All arguments are named; no positional arguments are used.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Build a master fertility genetics Excel workbook from results folders."
    )
    parser.add_argument(
        "--base_dir",
        required=True,
        type=Path,
        help="Base project directory, e.g. /home/pthorpe001/data/2026_sperm_Gates",
    )
    parser.add_argument(
        "--results_dir",
        required=False,
        type=Path,
        default=None,
        help=(
            "Results directory root. Defaults to <base_dir>/results if not provided."
        ),
    )
    parser.add_argument(
        "--testis_run_id",
        required=True,
        type=str,
        help=(
            "Run folder name for the testis pipeline under results/, "
            "e.g. testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1"
        ),
    )
    parser.add_argument(
        "--male_infertility_run_relpath",
        required=True,
        type=Path,
        help=(
            "Relative path (from base_dir) to male-infertility results folder, "
            "e.g. hpo_data/male_infertility_out_phrase_seeds OR results/01_hpo style folder."
        ),
    )
    parser.add_argument(
        "--clinvar_results_relpath",
        required=True,
        type=Path,
        help=(
            "Relative path (from base_dir) to the ClinVar filtered outputs folder, "
            "e.g. clinvar/clinvar_filtered OR results/..../05_clinvar_high_confidence."
        ),
    )
    parser.add_argument(
        "--out_xlsx",
        required=True,
        type=Path,
        help="Output Excel workbook path (.xlsx).",
    )
    parser.add_argument(
        "--include_variant_sheets",
        required=False,
        action="store_true",
        help="If set, include variant-level sheets (can make workbook larger).",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    """
    Read a tab-separated file safely.

    Parameters
    ----------
    path : Path
        Path to TSV.

    Returns
    -------
    pd.DataFrame
        DataFrame with strings preserved where possible.
    """
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def normalise_gene_symbol(value: str) -> str:
    """
    Normalise a gene symbol for consistent joining.

    Parameters
    ----------
    value : str
        Input symbol.

    Returns
    -------
    str
        Normalised symbol.
    """
    if value is None:
        return ""
    v = str(value).strip()
    v = re.sub(r"\s+", "", v)
    return v.upper()


def ensure_gene_key(df: pd.DataFrame, candidate_cols: Iterable[str]) -> pd.DataFrame:
    """
    Ensure a DataFrame contains a 'gene_key' column built from the first
    existing column in candidate_cols.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    candidate_cols : Iterable[str]
        Column names to try in order.

    Returns
    -------
    pd.DataFrame
        DataFrame with a 'gene_key' column.
    """
    for col in candidate_cols:
        if col in df.columns:
            df = df.copy()
            df["gene_key"] = df[col].fillna("").map(normalise_gene_symbol)
            return df
    raise ValueError(
        f"Could not find any gene symbol column in: {list(candidate_cols)}. "
        f"Available columns: {list(df.columns)}"
    )


def summarise_clinvar_by_gene(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    """
    Summarise ClinVar variant rows to gene-level counts by key categories.

    Parameters
    ----------
    df : pd.DataFrame
        Variant-level ClinVar DataFrame.
    gene_col : str
        Column containing gene symbols.

    Returns
    -------
    pd.DataFrame
        Gene-level summary DataFrame.
    """
    work = df.copy()
    if gene_col not in work.columns:
        raise ValueError(f"ClinVar gene column not found: {gene_col}")

    work = ensure_gene_key(work, candidate_cols=[gene_col])

    def _count(series: pd.Series) -> int:
        return int(series.dropna().shape[0])

    group = work.groupby("gene_key", dropna=False)
    out = pd.DataFrame(
        {
            "clinvar_rows": group.size().astype(int),
            "clinvar_unique_allele_ids": group["AlleleID"].nunique(dropna=True)
            if "AlleleID" in work.columns
            else group.size().astype(int),
            "clinvar_unique_variation_ids": group["VariationID"].nunique(dropna=True)
            if "VariationID" in work.columns
            else group.size().astype(int),
        }
    ).reset_index()

    if "ClinicalSignificance" in work.columns:
        cs = (
            work.groupby(["gene_key", "ClinicalSignificance"])
            .size()
            .reset_index(name="n")
        )
        cs_pivot = cs.pivot_table(
            index="gene_key", columns="ClinicalSignificance", values="n", fill_value=0
        )
        cs_pivot.columns = [f"clinvar_cs__{c}" for c in cs_pivot.columns]
        cs_pivot = cs_pivot.reset_index()
        out = out.merge(cs_pivot, on="gene_key", how="left")

    if "ReviewStatus" in work.columns:
        rs = (
            work.groupby(["gene_key", "ReviewStatus"])
            .size()
            .reset_index(name="n")
        )
        rs_pivot = rs.pivot_table(
            index="gene_key", columns="ReviewStatus", values="n", fill_value=0
        )
        rs_pivot.columns = [f"clinvar_rs__{c}" for c in rs_pivot.columns]
        rs_pivot = rs_pivot.reset_index()
        out = out.merge(rs_pivot, on="gene_key", how="left")

    return out


def write_sheet(
    wb: Workbook,
    title: str,
    df: pd.DataFrame,
    freeze_panes: str = "A2",
    autofilter: bool = True,
) -> None:
    """
    Write a DataFrame to an Excel sheet with basic formatting.

    Parameters
    ----------
    wb : Workbook
        Openpyxl workbook.
    title : str
        Sheet name.
    df : pd.DataFrame
        Data to write.
    freeze_panes : str
        Cell reference to freeze panes.
    autofilter : bool
        Whether to add an auto-filter to the header row.
    """
    ws = wb.create_sheet(title=title)

    header_fill = PatternFill("solid", fgColor="1F4E79")
    header_font = Font(bold=True, color="FFFFFF")
    header_align = Alignment(vertical="center", wrap_text=True)

    # Header
    ws.append(list(df.columns))
    for col_idx in range(1, len(df.columns) + 1):
        cell = ws.cell(row=1, column=col_idx)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = header_align

    # Rows
    for _, row in df.iterrows():
        ws.append([row.get(c, "") for c in df.columns])

    ws.freeze_panes = freeze_panes

    if autofilter:
        ws.auto_filter.ref = ws.dimensions

    # Column widths (simple heuristic)
    for col_idx, col_name in enumerate(df.columns, start=1):
        max_len = max(
            len(str(col_name)),
            *(len(str(v)) for v in df[col_name].astype(str).head(200).values),
        )
        width = min(max(10, max_len + 2), 60)
        ws.column_dimensions[get_column_letter(col_idx)].width = width


def main() -> None:
    """
    Main entry point.
    """
    args = parse_args()

    base_dir: Path = args.base_dir
    results_dir: Path = args.results_dir if args.results_dir else base_dir / "results"

    testis_run_dir = results_dir / args.testis_run_id
    if not testis_run_dir.exists():
        raise FileNotFoundError(f"Testis run folder not found: {testis_run_dir}")

    male_inf_dir = base_dir / args.male_infertility_run_relpath
    if not male_inf_dir.exists():
        raise FileNotFoundError(f"Male infertility folder not found: {male_inf_dir}")

    clinvar_dir = base_dir / args.clinvar_results_relpath
    if not clinvar_dir.exists():
        raise FileNotFoundError(f"ClinVar results folder not found: {clinvar_dir}")

    # Testis pipeline key files (from your run folder)
    testis_annot_path = (
        testis_run_dir
        / "08_proteomics_annotation"
        / "gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv"
    )
    if not testis_annot_path.exists():
        raise FileNotFoundError(
            "Could not find the fully annotated testis table at:\n"
            f"{testis_annot_path}\n"
            "Check your run folder structure and filenames."
        )

    tissue_medians_path = (
        testis_run_dir / "01_gtex_specificity" / "gtex_v11_tissue_medians.tsv"
    )

    # Male infertility (HPO-derived) gene summary
    # Supports either:
    #  - male_infertility_genes_summary.tsv in the male_inf_dir
    #  - or nested under results/01_hpo/
    candidates = [
        male_inf_dir / "male_infertility_genes_summary.tsv",
        male_inf_dir / "01_hpo" / "male_infertility_genes_summary.tsv",
    ]
    hpo_genes_summary_path = next((p for p in candidates if p.exists()), None)
    if hpo_genes_summary_path is None:
        raise FileNotFoundError(
            "Could not find male_infertility_genes_summary.tsv under:\n"
            f"{male_inf_dir}\n"
            "Expected either direct file or under 01_hpo/."
        )

    # ClinVar tier outputs (support either original names or your results-renamed files)
    def pick_first(paths: List[Path], label: str) -> Path:
        for p in paths:
            if p.exists():
                return p
        raise FileNotFoundError(
            f"Could not find {label}. Tried:\n" + "\n".join(str(p) for p in paths)
        )

    clinvar_best_path = pick_first(
        paths=[
            clinvar_dir / "male_infertility_clinvar_variants_best.tsv",
            clinvar_dir / "clinvar_best.tsv",
        ],
        label="ClinVar best TSV",
    )
    clinvar_hc_path = pick_first(
        paths=[
            clinvar_dir / "male_infertility_clinvar_variants_high_confidence.tsv",
            clinvar_dir / "clinvar_high_confidence.tsv",
        ],
        label="ClinVar high-confidence TSV",
    )
    clinvar_hc_pathogenic_path = pick_first(
        paths=[
            clinvar_dir / "male_infertility_clinvar_variants_high_confidence_pathogenic.tsv",
            clinvar_dir / "clinvar_high_confidence_pathogenic.tsv",
        ],
        label="ClinVar high-confidence pathogenic TSV",
    )

    # Load tables
    testis_df = read_tsv(testis_annot_path)
    testis_df = ensure_gene_key(
        testis_df,
        candidate_cols=[
            "gene_symbol_norm",
            "hgnc_symbol_norm",
            "hgnc_symbol",
            "hgnc_symbol_from_ensembl",
        ],
    )

    hpo_df = read_tsv(hpo_genes_summary_path)
    hpo_df = ensure_gene_key(hpo_df, candidate_cols=["gene_symbol", "GeneSymbol", "symbol"])

    clinvar_best_df = read_tsv(clinvar_best_path)
    clinvar_best_df = ensure_gene_key(clinvar_best_df, candidate_cols=["GeneSymbol"])

    clinvar_hc_df = read_tsv(clinvar_hc_path)
    clinvar_hc_df = ensure_gene_key(clinvar_hc_df, candidate_cols=["GeneSymbol"])

    clinvar_hc_path_df = read_tsv(clinvar_hc_pathogenic_path)
    clinvar_hc_path_df = ensure_gene_key(clinvar_hc_path_df, candidate_cols=["GeneSymbol"])

    # Gene-level ClinVar summaries
    clinvar_best_gene = summarise_clinvar_by_gene(clinvar_best_df, gene_col="GeneSymbol")
    clinvar_hc_gene = summarise_clinvar_by_gene(clinvar_hc_df, gene_col="GeneSymbol")
    clinvar_hc_path_gene = summarise_clinvar_by_gene(
        clinvar_hc_path_df, gene_col="GeneSymbol"
    )

    # Flag columns for tiers
    clinvar_best_gene["clinvar_best_present"] = True
    clinvar_hc_gene["clinvar_high_confidence_present"] = True
    clinvar_hc_path_gene["clinvar_high_confidence_pathogenic_present"] = True

    # Combine gene-level master
    genes_master = testis_df.copy()

    # Select a sensible subset of columns first (keep everything else later if needed)
    preferred_cols = [
        "gene_key",
        "hgnc_symbol",
        "hgnc_symbol_norm",
        "gene_symbol_norm",
        "ensembl_gene_id",
        "tau",
        "target_median_tpm",
        "target_present_fraction",
        "max_non_target_median_tpm",
        "log2_fc_target_vs_max_non_target",
        "sperm_present_any",
        "sperm_present_frac",
        "sperm_tpm_mean",
        "sperm_tpm_median",
        "hpo_term_count",
        "hpo_ids",
        "hpo_terms",
        "hpo_disease_ids",
        "prot_present_any",
        "prot_present_fraction",
        "prot_accessions",
        "prot_n_accessions",
        "prot_best_qvalue",
        "prot_unique_peptides_max",
        "prot_coverage_pct_max",
        "gene_type",
        "go_bp",
        "go_mf",
        "go_cc",
        "gene_summary",
    ]
    keep_cols = [c for c in preferred_cols if c in genes_master.columns]
    other_cols = [c for c in genes_master.columns if c not in keep_cols]
    genes_master = genes_master[keep_cols + other_cols]

    # Merge HPO gene summary fields (gene-level)
    hpo_cols = [c for c in hpo_df.columns if c != "gene_key"]
    hpo_rename = {c: f"hpo_genesummary__{c}" for c in hpo_cols}
    genes_master = genes_master.merge(
        hpo_df.rename(columns=hpo_rename),
        on="gene_key",
        how="left",
    )

    # Merge ClinVar summaries
    genes_master = genes_master.merge(clinvar_best_gene, on="gene_key", how="left")
    genes_master = genes_master.merge(clinvar_hc_gene, on="gene_key", how="left", suffixes=("", "_hc"))
    genes_master = genes_master.merge(
        clinvar_hc_path_gene, on="gene_key", how="left", suffixes=("", "_hc_path")
    )

    # Fill tier flags
    for col in [
        "clinvar_best_present",
        "clinvar_high_confidence_present",
        "clinvar_high_confidence_pathogenic_present",
    ]:
        if col in genes_master.columns:
            genes_master[col] = genes_master[col].fillna(False)

    # Workbook
    wb = Workbook()
    # Remove default sheet
    wb.remove(wb.active)

    # README sheet
    readme = pd.DataFrame(
        {
            "field": [
                "created",
                "base_dir",
                "results_dir",
                "testis_run_id",
                "male_infertility_relpath",
                "clinvar_results_relpath",
                "annotated_testis_table",
                "hpo_genes_summary",
                "clinvar_best",
                "clinvar_high_confidence",
                "clinvar_high_confidence_pathogenic",
                "include_variant_sheets",
            ],
            "value": [
                dt.datetime.now().isoformat(timespec="seconds"),
                str(base_dir),
                str(results_dir),
                args.testis_run_id,
                str(args.male_infertility_run_relpath),
                str(args.clinvar_results_relpath),
                str(testis_annot_path),
                str(hpo_genes_summary_path),
                str(clinvar_best_path),
                str(clinvar_hc_path),
                str(clinvar_hc_pathogenic_path),
                str(bool(args.include_variant_sheets)),
            ],
        }
    )
    write_sheet(wb, title="README", df=readme, freeze_panes="A2", autofilter=False)

    # Genes master
    write_sheet(wb, title="Genes_Master", df=genes_master, freeze_panes="A2", autofilter=True)

    # Optional supportive sheets
    if tissue_medians_path.exists():
        tissue_df = read_tsv(tissue_medians_path)
        write_sheet(wb, title="GTEx_Tissue_Medians", df=tissue_df, freeze_panes="A2")

    write_sheet(wb, title="HPO_Genes_Summary", df=hpo_df, freeze_panes="A2")

    # Variant sheets (optional, can be large)
    if args.include_variant_sheets:
        write_sheet(wb, title="ClinVar_Best", df=clinvar_best_df, freeze_panes="A2")
        write_sheet(wb, title="ClinVar_HighConfidence", df=clinvar_hc_df, freeze_panes="A2")
        write_sheet(
            wb,
            title="ClinVar_HC_Pathogenic",
            df=clinvar_hc_path_df,
            freeze_panes="A2",
        )

    # Save
    args.out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb.save(args.out_xlsx)


if __name__ == "__main__":
    main()
