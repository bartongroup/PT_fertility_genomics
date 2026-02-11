#!/usr/bin/env python3
"""
Create a concise, boss-ready report from high-confidence pathogenic ClinVar TSV.

Inputs:
- male_infertility_clinvar_variants_high_confidence_pathogenic.tsv

Outputs:
- report TSV with selected columns
- per-gene summary TSV with counts and phenotype strings

All outputs are TSV.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


REPORT_COLUMNS: List[str] = [
    "GeneSymbol",
    "Name",
    "ClinicalSignificance",
    "ReviewStatus",
    "PhenotypeList",
    "LastEvaluated",
    "RCVaccession",
    "RS# (dbSNP)",
    "Assembly",
    "Chromosome",
    "Start",
    "Stop",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "VariationID",
    "AlleleID",
]


def parse_args() -> argparse.Namespace:
    """
    Parse arguments.

    Returns:
        Parsed args.
    """
    parser = argparse.ArgumentParser(description="Make a concise report from ClinVar high-confidence pathogenic TSV.")
    parser.add_argument("--in_tsv", required=True, help="Input TSV (high confidence pathogenic).")
    parser.add_argument("--out_report_tsv", required=True, help="Output TSV (concise report).")
    parser.add_argument("--out_gene_summary_tsv", required=True, help="Output TSV (per-gene summary).")
    parser.add_argument(
        "--max_phenotypes_per_gene",
        type=int,
        default=10,
        help="Maximum number of phenotype strings to keep per gene in summary.",
    )
    return parser.parse_args()


def _norm(value: str) -> str:
    """
    Normalise values.

    Args:
        value: Raw string.

    Returns:
        Stripped string.
    """
    return (value or "").strip()


def _safe_get(row: Dict[str, str], key: str) -> str:
    """
    Safe field extraction.

    Args:
        row: TSV row.
        key: Column name.

    Returns:
        Field value, or empty string.
    """
    return _norm(row.get(key, ""))


def write_report(
    rows: List[Dict[str, str]],
    out_path: Path,
) -> None:
    """
    Write concise report TSV.

    Args:
        rows: Input rows.
        out_path: Output TSV path.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(REPORT_COLUMNS)
        for row in rows:
            writer.writerow([_safe_get(row, c) for c in REPORT_COLUMNS])


def summarise_by_gene(
    rows: List[Dict[str, str]],
    max_phenotypes_per_gene: int,
) -> List[Dict[str, str]]:
    """
    Summarise rows by gene.

    Args:
        rows: Input rows.
        max_phenotypes_per_gene: Max phenotype strings to retain.

    Returns:
        List of summary dicts.
    """
    gene_to_variants: Dict[str, Set[str]] = defaultdict(set)
    gene_to_phenotypes: Dict[str, Counter[str]] = defaultdict(Counter)
    gene_to_review: Dict[str, Counter[str]] = defaultdict(Counter)

    for row in rows:
        gene = _safe_get(row, "GeneSymbol")
        if not gene:
            continue

        allele_id = _safe_get(row, "AlleleID")
        variation_id = _safe_get(row, "VariationID")
        key = allele_id or variation_id or ""

        if key:
            gene_to_variants[gene].add(key)

        phenos = _safe_get(row, "PhenotypeList")
        if phenos:
            # Split on pipes, keep unique strings
            for p in [p.strip() for p in phenos.split("|") if p.strip()]:
                gene_to_phenotypes[gene][p] += 1

        review = _safe_get(row, "ReviewStatus")
        if review:
            gene_to_review[gene][review] += 1

    summaries: List[Dict[str, str]] = []
    for gene in sorted(gene_to_variants.keys()):
        top_phenos = []
        for p, _ in gene_to_phenotypes[gene].most_common():
            if p not in top_phenos:
                top_phenos.append(p)
            if len(top_phenos) >= max_phenotypes_per_gene:
                break

        
        top_reviews = [r for r, _ in gene_to_review[gene].most_common(5)]

        summaries.append(
            {
                "GeneSymbol": gene,
                "n_variants": str(len(gene_to_variants[gene])),
                "top_phenotypes": ";".join(top_phenos),
                "review_statuses": ";".join(top_reviews),
            }
        )

    return summaries


def write_gene_summary(rows: List[Dict[str, str]], out_path: Path, max_phenotypes_per_gene: int) -> None:
    """
    Write per-gene summary TSV.

    Args:
        rows: Input rows.
        out_path: Output TSV path.
        max_phenotypes_per_gene: Max phenotypes per gene.
    """
    summary_rows = summarise_by_gene(rows=rows, max_phenotypes_per_gene=max_phenotypes_per_gene)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["GeneSymbol", "n_variants", "top_phenotypes", "review_statuses"])
        for r in summary_rows:
            writer.writerow([r["GeneSymbol"], r["n_variants"], r["top_phenotypes"], r["review_statuses"]])


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()
    in_tsv = Path(args.in_tsv)
    out_report_tsv = Path(args.out_report_tsv)
    out_gene_summary_tsv = Path(args.out_gene_summary_tsv)

    if not in_tsv.exists():
        raise FileNotFoundError(f"Missing input: {in_tsv}")

    rows: List[Dict[str, str]] = []
    with in_tsv.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")
        for row in reader:
            rows.append(row)

    # Sort for readability: GeneSymbol then Name
    rows_sorted = sorted(rows, key=lambda r: (_safe_get(r, "GeneSymbol"), _safe_get(r, "Name")))

    write_report(rows=rows_sorted, out_path=out_report_tsv)
    write_gene_summary(rows=rows_sorted, out_path=out_gene_summary_tsv, max_phenotypes_per_gene=int(args.max_phenotypes_per_gene))


if __name__ == "__main__":
    main()

