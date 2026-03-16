
#!/usr/bin/env python3
"""
Analyse fertility gene lists from a column-based Excel workbook.

Outputs TSV files for:
1. Recovery of validated infertility genes across selected gene sets.
2. Overlap with a curated sperm flagellar motility gene module.
3. Novel candidates in the strict prioritised set.
4. Novel candidates shared by the strict prioritised and druggable sets.

Usage:
    python fertility_gene_list_validation.py \
        --input SUMMARY_fertility_evidence_gene_lists_only.xlsx \
        --output-dir results
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Iterable, Set

import pandas as pd


LOGGER = logging.getLogger(__name__)


FLAGELLAR_MODULE = {
    "AK7",
    "CCDC39",
    "CCDC40",
    "CCDC42",
    "CCDC63",
    "CCDC65",
    "CCDC103",
    "CFAP43",
    "CFAP44",
    "CFAP45",
    "CFAP46",
    "CFAP52",
    "CFAP54",
    "CFAP58",
    "CFAP61",
    "CFAP69",
    "CFAP70",
    "CFAP91",
    "CFAP157",
    "CFAP221",
    "CFAP298",
    "DNAH1",
    "DNAH2",
    "DNAH5",
    "DNAH6",
    "DNAH7",
    "DNAH8",
    "DNAH9",
    "DNAH10",
    "DNAH11",
    "DNAH12",
    "DNAH17",
    "DNAI1",
    "DNAI2",
    "DNAAF1",
    "DNAAF2",
    "DNAAF3",
    "DNAAF4",
    "DNAAF5",
    "DNAAF6",
    "DNALI1",
    "DYNLRB1",
    "DYNLRB2",
    "LRRC6",
    "ODF1",
    "ODF2",
    "RSPH1",
    "RSPH3",
    "RSPH4A",
    "RSPH6A",
    "RSPH9",
    "SPAG1",
    "SPAG6",
    "SPAG16",
    "SPAG17",
    "ZMYND10",
}


DEFAULT_TARGET_SETS = [
    "Testis_high_confidence_final",
    "Sperm_RNAseq_present",
    "Proteomics_any_source_present",
    "Strict_final_GOI_insperm",
    "druggable_targets",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate fertility gene lists and export TSV summaries."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to the input Excel workbook containing column-based gene lists.",
    )
    parser.add_argument(
        "--sheet",
        default="Gene_Lists",
        help="Workbook sheet containing the gene lists.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for TSV outputs.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def clean_gene_set(values: Iterable[object]) -> Set[str]:
    genes: Set[str] = set()

    for value in values:
        if pd.isna(value):
            continue

        gene = str(value).strip().upper()

        if gene and gene != "NAN":
            genes.add(gene)

    return genes


def load_gene_sets(input_path: Path, sheet_name: str) -> Dict[str, Set[str]]:
    LOGGER.info("Loading workbook: %s", input_path)
    dataframe = pd.read_excel(input_path, sheet_name=sheet_name)

    gene_sets = {
        column: clean_gene_set(dataframe[column])
        for column in dataframe.columns
    }

    LOGGER.info("Loaded %d gene-list columns", len(gene_sets))
    return gene_sets


def make_validated_infertility_set(gene_sets: Dict[str, Set[str]]) -> Set[str]:
    validated = (
        gene_sets["HPO_gene_set"]
        | gene_sets["ClinVar_best_pathogenic"]
        | gene_sets["ClinVar_high_confidence_pathogenic"]
    )
    LOGGER.info(
        "Validated infertility reference set contains %d genes",
        len(validated),
    )
    return validated


def build_recovery_summary(
    gene_sets: Dict[str, Set[str]],
    validated_set: Set[str],
) -> pd.DataFrame:
    rows = []

    for set_name in DEFAULT_TARGET_SETS:
        current_set = gene_sets[set_name]
        overlap = sorted(current_set & validated_set)

        rows.append(
            {
                "gene_set": set_name,
                "total_genes": len(current_set),
                "validated_infertility_genes_recovered": len(overlap),
                "validated_recovery_fraction": (
                    len(overlap) / len(validated_set)
                    if validated_set
                    else 0.0
                ),
                "validated_recovery_percent": (
                    100.0 * len(overlap) / len(validated_set)
                    if validated_set
                    else 0.0
                ),
                "overlap_genes": ";".join(overlap),
            }
        )

    return pd.DataFrame(rows).sort_values(
        by="validated_infertility_genes_recovered",
        ascending=False,
    )


def build_flagellar_summary(gene_sets: Dict[str, Set[str]]) -> pd.DataFrame:
    rows = []

    for set_name in DEFAULT_TARGET_SETS:
        current_set = gene_sets[set_name]
        overlap = sorted(current_set & FLAGELLAR_MODULE)

        rows.append(
            {
                "gene_set": set_name,
                "total_genes": len(current_set),
                "flagellar_module_genes_recovered": len(overlap),
                "flagellar_module_size": len(FLAGELLAR_MODULE),
                "flagellar_recovery_fraction": len(overlap) / len(FLAGELLAR_MODULE),
                "flagellar_recovery_percent": (
                    100.0 * len(overlap) / len(FLAGELLAR_MODULE)
                ),
                "overlap_genes": ";".join(overlap),
            }
        )

    return pd.DataFrame(rows).sort_values(
        by="flagellar_module_genes_recovered",
        ascending=False,
    )


def build_novel_candidate_tables(
    gene_sets: Dict[str, Set[str]],
    validated_set: Set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    literature_set = gene_sets["Literature_gene_set"]
    strict_set = gene_sets["Strict_final_GOI_insperm"]
    druggable_set = gene_sets["druggable_targets"]

    novel_strict = sorted(strict_set - validated_set - literature_set)
    novel_strict_and_druggable = sorted(
        (strict_set & druggable_set) - validated_set - literature_set
    )

    novel_strict_df = pd.DataFrame(
        {
            "gene": novel_strict,
            "in_strict_final_goi_insperm": True,
            "in_druggable_targets": [
                gene in druggable_set for gene in novel_strict
            ],
            "in_validated_infertility_set": False,
            "in_literature_gene_set": False,
        }
    )

    novel_shared_df = pd.DataFrame(
        {
            "gene": novel_strict_and_druggable,
            "in_strict_final_goi_insperm": True,
            "in_druggable_targets": True,
            "in_validated_infertility_set": False,
            "in_literature_gene_set": False,
        }
    )

    LOGGER.info(
        "Novel genes in strict prioritised set: %d",
        len(novel_strict_df),
    )
    LOGGER.info(
        "Novel genes shared by strict prioritised and druggable sets: %d",
        len(novel_shared_df),
    )

    return novel_strict_df, novel_shared_df


def write_tsv(dataframe: pd.DataFrame, output_path: Path) -> None:
    LOGGER.info("Writing TSV: %s", output_path)
    dataframe.to_csv(output_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    gene_sets = load_gene_sets(args.input, args.sheet)
    validated_set = make_validated_infertility_set(gene_sets)

    recovery_summary = build_recovery_summary(gene_sets, validated_set)
    flagellar_summary = build_flagellar_summary(gene_sets)
    novel_strict_df, novel_shared_df = build_novel_candidate_tables(
        gene_sets,
        validated_set,
    )

    write_tsv(
        recovery_summary,
        args.output_dir / "validated_infertility_recovery_summary.tsv",
    )
    write_tsv(
        flagellar_summary,
        args.output_dir / "flagellar_module_recovery_summary.tsv",
    )
    write_tsv(
        novel_strict_df,
        args.output_dir / "novel_candidates_strict_final.tsv",
    )
    write_tsv(
        novel_shared_df,
        args.output_dir / "novel_candidates_strict_and_druggable.tsv",
    )

    LOGGER.info("Analysis complete")


if __name__ == "__main__":
    main()
