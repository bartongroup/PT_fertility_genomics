#!/usr/bin/env python3
"""
Overlap testis-specific genes with ClinVar male infertility genes, tracking evidence tiers.

Inputs:
- testis_specific_gene_symbols.tsv (TSV with 'gene_symbol')
- male_infertility_clinvar_variants_best.tsv (best row per AlleleID; infertility phenotype filtered)
- male_infertility_clinvar_variants_high_confidence_pathogenic.tsv (subset; high-confidence P/LP)

Outputs:
- overlap summary with tiering and counts
- testis-only gene list
- clinvar-only gene list
- optional long evidence table (gene-variant rows with tier)

All outputs are TSV.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple


HIGH_CONF_REVIEW = {
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
}

PATHOGENICISH = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic, low penetrance",
}


def parse_args() -> argparse.Namespace:
    """
    Parse arguments.

    Returns:
        Parsed args.
    """
    parser = argparse.ArgumentParser(description="Overlap testis-specific genes with ClinVar infertility tiers.")
    parser.add_argument("--testis_gene_tsv", required=True, help="TSV with column 'gene_symbol'.")
    parser.add_argument(
        "--clinvar_best_tsv",
        required=True,
        help="TSV: male_infertility_clinvar_variants_best.tsv",
    )
    parser.add_argument(
        "--clinvar_hc_pathogenic_tsv",
        required=True,
        help="TSV: male_infertility_clinvar_variants_high_confidence_pathogenic.tsv",
    )
    parser.add_argument(
        "--gene_symbol_column",
        default="gene_symbol",
        help=(
            "Column name in testis_gene_tsv that contains gene symbols. "
            "Example: gene_symbol_norm or hgnc_symbol."
        ),
    )

    parser.add_argument("--out_dir", required=True, help="Output directory.")
    return parser.parse_args()


def _norm(value: str) -> str:
    """
    Normalise a value.

    Args:
        value: Raw value.

    Returns:
        Stripped value.
    """
    return (value or "").strip()


def _norm_upper(value: str) -> str:
    """
    Normalise and uppercase.

    Args:
        value: Raw value.

    Returns:
        Uppercased value.
    """
    return _norm(value).upper()


def _norm_lower(value: str) -> str:
    """
    Normalise and lowercase.

    Args:
        value: Raw value.

    Returns:
        Lowercased value.
    """
    return _norm(value).lower()



def load_gene_symbols(path: Path, gene_symbol_column: str) -> Set[str]:
    """
    Load gene symbols from a TSV.

    Args:
        path: TSV path.
        gene_symbol_column: Column containing gene symbols.

    Returns:
        Set of gene symbols (uppercased).
    """
    genes: Set[str] = set()
    with path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")

        if gene_symbol_column not in reader.fieldnames:
            raise ValueError(
                f"testis_gene_tsv must contain a column named '{gene_symbol_column}'. "
                f"Available columns: {reader.fieldnames}"
            )

        for row in reader:
            g = _norm_upper(row.get(gene_symbol_column, ""))
            if g:
                genes.add(g)

    return genes


def is_high_conf(review_status: str) -> bool:
    """
    High-confidence ClinVar review status.

    Args:
        review_status: ReviewStatus string.

    Returns:
        True if high confidence.
    """
    return _norm_lower(review_status) in HIGH_CONF_REVIEW


def is_pathogenicish(clinsig: str) -> bool:
    """
    Determine pathogenic-ish clinical significance.

    Args:
        clinsig: ClinicalSignificance string.

    Returns:
        True if pathogenic-ish.
    """
    return _norm_lower(clinsig) in PATHOGENICISH


def tier_for_row(review_status: str, clinsig: str) -> str:
    """
    Assign evidence tier for a ClinVar row.

    Args:
        review_status: ClinVar ReviewStatus.
        clinsig: ClinVar ClinicalSignificance.

    Returns:
        Tier label.
    """
    hc = is_high_conf(review_status)
    pth = is_pathogenicish(clinsig)

    if hc and pth:
        return "tier1_hc_pathogenic"
    if hc and not pth:
        return "tier2_hc_other"
    if (not hc) and pth:
        return "tier3_lc_pathogenic"
    return "tier4_lc_other"


TIER_ORDER = {
    "tier1_hc_pathogenic": 1,
    "tier2_hc_other": 2,
    "tier3_lc_pathogenic": 3,
    "tier4_lc_other": 4,
}


def choose_best_tier(tiers: Set[str]) -> str:
    """
    Choose the best tier present.

    Args:
        tiers: Set of tiers.

    Returns:
        Best tier label.
    """
    if not tiers:
        return ""
    return sorted(tiers, key=lambda t: TIER_ORDER.get(t, 99))[0]


def load_clinvar_rows(path: Path) -> List[Dict[str, str]]:
    """
    Load ClinVar TSV rows.

    Args:
        path: TSV path.

    Returns:
        Rows.
    """
    rows: List[Dict[str, str]] = []
    with path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in: {path}")
        for row in reader:
            rows.append(row)
    return rows


def main() -> None:
    """
    Main.
    """
    args = parse_args()

    testis_genes = load_gene_symbols(
            path=Path(args.testis_gene_tsv),
            gene_symbol_column=str(args.gene_symbol_column),
        )

    clinvar_best_rows = load_clinvar_rows(path=Path(args.clinvar_best_tsv))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Aggregate evidence per gene from best rows (infertility phenotype-filtered already)
    gene_to_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    gene_to_tiers: Dict[str, Set[str]] = defaultdict(set)
    gene_to_example: Dict[str, Dict[str, str]] = {}

    for row in clinvar_best_rows:
        gene = _norm_upper(row.get("GeneSymbol", ""))
        if not gene:
            continue

        review = _norm(row.get("ReviewStatus", ""))
        clinsig = _norm(row.get("ClinicalSignificance", ""))
        tier = tier_for_row(review_status=review, clinsig=clinsig)

        gene_to_counts[gene][tier] += 1
        gene_to_tiers[gene].add(tier)

        # Store a representative example: prefer best tier
        if gene not in gene_to_example:
            gene_to_example[gene] = row
        else:
            current_best = choose_best_tier({tier_for_row(gene_to_example[gene].get("ReviewStatus", ""),
                                                         gene_to_example[gene].get("ClinicalSignificance", ""))})
            if TIER_ORDER.get(tier, 99) < TIER_ORDER.get(current_best, 99):
                gene_to_example[gene] = row

    clinvar_genes = set(gene_to_counts.keys())

    overlap = sorted(testis_genes.intersection(clinvar_genes))
    testis_only = sorted(testis_genes.difference(clinvar_genes))
    clinvar_only = sorted(clinvar_genes.difference(testis_genes))

    # Write overlap summary
    overlap_path = out_dir / "overlap_testis_vs_clinvar_genes.tsv"
    with overlap_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_symbol",
                "tier_best",
                "n_hc_pathogenic_variants",
                "n_hc_other_variants",
                "n_lc_pathogenic_variants",
                "n_lc_other_variants",
                "example_name",
                "example_clinsig",
                "example_review_status",
                "example_phenotypes",
            ]
        )
        for gene in overlap:
            counts = gene_to_counts[gene]
            tiers = gene_to_tiers[gene]
            best_tier = choose_best_tier(tiers)

            example = gene_to_example.get(gene, {})
            writer.writerow(
                [
                    gene,
                    best_tier,
                    str(counts.get("tier1_hc_pathogenic", 0)),
                    str(counts.get("tier2_hc_other", 0)),
                    str(counts.get("tier3_lc_pathogenic", 0)),
                    str(counts.get("tier4_lc_other", 0)),
                    _norm(example.get("Name", "")),
                    _norm(example.get("ClinicalSignificance", "")),
                    _norm(example.get("ReviewStatus", "")),
                    _norm(example.get("PhenotypeList", "")),
                ]
            )

    # Write gene lists
    def write_gene_list(path: Path, genes: List[str]) -> None:
        with path.open(mode="w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene_symbol"])
            for g in genes:
                writer.writerow([g])

    write_gene_list(out_dir / "testis_only_genes.tsv", testis_only)
    write_gene_list(out_dir / "clinvar_only_genes.tsv", clinvar_only)

    # Write long evidence table (optional but useful for tracking)
    long_path = out_dir / "clinvar_gene_evidence_long.tsv"
    with long_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_symbol",
                "tier",
                "AlleleID",
                "VariationID",
                "Name",
                "ClinicalSignificance",
                "ReviewStatus",
                "LastEvaluated",
                "RCVaccession",
                "RS# (dbSNP)",
                "PhenotypeList",
                "Assembly",
                "Chromosome",
                "Start",
                "Stop",
                "ReferenceAlleleVCF",
                "AlternateAlleleVCF",
            ]
        )
        for row in clinvar_best_rows:
            gene = _norm_upper(row.get("GeneSymbol", ""))
            if not gene:
                continue
            tier = tier_for_row(
                review_status=_norm(row.get("ReviewStatus", "")),
                clinsig=_norm(row.get("ClinicalSignificance", "")),
            )
            writer.writerow(
                [
                    gene,
                    tier,
                    _norm(row.get("AlleleID", row.get("#AlleleID", ""))),
                    _norm(row.get("VariationID", "")),
                    _norm(row.get("Name", "")),
                    _norm(row.get("ClinicalSignificance", "")),
                    _norm(row.get("ReviewStatus", "")),
                    _norm(row.get("LastEvaluated", "")),
                    _norm(row.get("RCVaccession", "")),
                    _norm(row.get("RS# (dbSNP)", "")),
                    _norm(row.get("PhenotypeList", "")),
                    _norm(row.get("Assembly", "")),
                    _norm(row.get("Chromosome", "")),
                    _norm(row.get("Start", "")),
                    _norm(row.get("Stop", "")),
                    _norm(row.get("ReferenceAlleleVCF", "")),
                    _norm(row.get("AlternateAlleleVCF", "")),
                ]
            )


if __name__ == "__main__":
    main()
