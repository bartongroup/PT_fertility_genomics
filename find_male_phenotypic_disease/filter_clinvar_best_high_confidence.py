#!/usr/bin/env python3
"""
Filter a best-per-variant ClinVar TSV to high-confidence records.

Input:
- male_infertility_clinvar_variants_best.tsv (one row per variant, chosen by ranking)

Outputs:
- high_confidence.tsv: rows with strong ReviewStatus values
- high_confidence_pathogenic.tsv: high confidence AND pathogenic-ish
- counts.tsv: summary counts

All outputs are TSV.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple


HIGH_CONF_REVIEW: Set[str] = {
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
}


PATHOGENICISH: Set[str] = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "likely pathogenic, low penetrance",
}


def parse_args() -> argparse.Namespace:
    """
    Parse CLI arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Filter ClinVar best TSV to high-confidence records.")
    parser.add_argument("--in_tsv", required=True, help="Input TSV (best per variant).")
    parser.add_argument("--out_high_conf_tsv", required=True, help="Output TSV (high confidence).")
    parser.add_argument(
        "--out_high_conf_pathogenic_tsv",
        required=True,
        help="Output TSV (high confidence + pathogenic-ish).",
    )
    parser.add_argument("--out_counts_tsv", required=True, help="Output TSV (summary counts).")
    return parser.parse_args()


def _norm(value: str) -> str:
    """
    Normalise strings.

    Args:
        value: Input string.

    Returns:
        Normalised string.
    """
    return (value or "").strip()


def _norm_lower(value: str) -> str:
    """
    Normalise strings to lower-case.

    Args:
        value: Input string.

    Returns:
        Normalised lower-case string.
    """
    return _norm(value).lower()


def is_high_confidence(review_status: str) -> bool:
    """
    Determine if ReviewStatus is considered high confidence.

    Args:
        review_status: ReviewStatus field.

    Returns:
        True if in high-confidence set.
    """
    return _norm_lower(review_status) in HIGH_CONF_REVIEW


def is_pathogenicish(clinsig: str) -> bool:
    """
    Determine if ClinicalSignificance is pathogenic-ish.

    Args:
        clinsig: ClinicalSignificance field.

    Returns:
        True if pathogenic-ish.
    """
    return _norm_lower(clinsig) in PATHOGENICISH


def write_rows(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    """
    Write TSV output.

    Args:
        path: Output path.
        fieldnames: Column names.
        rows: Row dicts.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_counts(path: Path, label: str, rows: List[Dict[str, str]]) -> None:
    """
    Write summary counts.

    Args:
        path: Output TSV path.
        label: Label for this dataset.
        rows: Rows to summarise.
    """
    clinsig_counts: Counter[str] = Counter()
    review_counts: Counter[str] = Counter()
    gene_counts: Counter[str] = Counter()

    for row in rows:
        clinsig_counts[_norm(row.get("ClinicalSignificance", ""))] += 1
        review_counts[_norm(row.get("ReviewStatus", ""))] += 1
        gene_counts[_norm(row.get("GeneSymbol", ""))] += 1

    with path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["dataset", "metric", "category", "count"])
        writer.writerow([label, "rows", "kept", len(rows)])

        for k, v in clinsig_counts.most_common():
            writer.writerow([label, "ClinicalSignificance", k, v])

        for k, v in review_counts.most_common():
            writer.writerow([label, "ReviewStatus", k, v])

        writer.writerow([label, "genes", "unique_genes", len([g for g in gene_counts if g])])


def main() -> None:
    """
    Main entry point.
    """
    args = parse_args()
    in_tsv = Path(args.in_tsv)

    if not in_tsv.exists():
        raise FileNotFoundError(f"Missing input: {in_tsv}")

    with in_tsv.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")
        fieldnames = reader.fieldnames

        high_conf: List[Dict[str, str]] = []
        high_conf_path: List[Dict[str, str]] = []

        for row in reader:
            if not is_high_confidence(row.get("ReviewStatus", "")):
                continue
            high_conf.append(row)

            if is_pathogenicish(row.get("ClinicalSignificance", "")):
                high_conf_path.append(row)

    write_rows(path=Path(args.out_high_conf_tsv), fieldnames=fieldnames, rows=high_conf)
    write_rows(path=Path(args.out_high_conf_pathogenic_tsv), fieldnames=fieldnames, rows=high_conf_path)

    Path(args.out_counts_tsv).parent.mkdir(parents=True, exist_ok=True)
    write_counts(path=Path(args.out_counts_tsv), label="high_confidence", rows=high_conf)

    # Also append counts for pathogenic subset in same file (simple append)
    with Path(args.out_counts_tsv).open(mode="a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([])
    write_counts(path=Path(args.out_counts_tsv), label="high_confidence_pathogenic", rows=high_conf_path)


if __name__ == "__main__":
    main()
