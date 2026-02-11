#!/usr/bin/env python3
"""
Collapse a ClinVar infertility-filtered TSV to one best row per variant.

Input:
- TSV from the phenotype filter step, e.g. male_infertility_clinvar_variants_by_phenotype.tsv
  (can be the full or minimal output; this script expects key columns to exist)

Method:
- Choose a single "best" row per AlleleID (default grouping key).
- Prefer GRCh38 over GRCh37 if both exist.
- Rank ReviewStatus by strength.
- Rank ClinicalSignificance by interpretation priority.

Outputs:
- best_per_variant.tsv (one row per variant)
- best_per_variant_pathogenic.tsv (Pathogenic/Likely pathogenic only)
- summary_counts.tsv

All outputs are TSV.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


REVIEW_STATUS_RANK: Dict[str, int] = {
    "practice guideline": 1,
    "reviewed by expert panel": 2,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, conflicting classifications": 4,
    "criteria provided, single submitter": 5,
    "no assertion criteria provided": 6,
    "no assertion provided": 7,
    "no interpretation for the single variant": 8,
    "": 99,
}

# Lower is better
CLINSIG_RANK: Dict[str, int] = {
    "pathogenic": 1,
    "likely pathogenic": 2,
    "pathogenic/likely pathogenic": 3,
    "uncertain significance": 4,
    "conflicting classifications of pathogenicity": 5,
    "likely benign": 6,
    "benign": 7,
    "benign/likely benign": 8,
    "drug response": 20,
    "risk factor": 21,
    "association": 22,
    "protective": 23,
    "affects": 24,
    "not provided": 90,
    "": 99,
}


@dataclass(frozen=True, order=True)
class RowScore:
    """
    Score used for choosing the best row for a given variant.

    Lower values are preferred.

    Attributes:
        assembly_rank: Prefer GRCh38 (0) over GRCh37 (1) over other/blank (2).
        review_rank: Strength of review status.
        clinsig_rank: Priority of clinical significance.
    """
    assembly_rank: int
    review_rank: int
    clinsig_rank: int


def parse_args() -> argparse.Namespace:
    """
    Parse arguments.

    Returns:
        Parsed args.
    """
    parser = argparse.ArgumentParser(description="Collapse ClinVar TSV to best row per variant.")
    parser.add_argument("--in_tsv", required=True, help="Input TSV (variants_by_phenotype).")
    parser.add_argument("--out_best_tsv", required=True, help="Output TSV (best per variant).")
    parser.add_argument(
        "--out_best_pathogenic_tsv",
        required=True,
        help="Output TSV (best per variant, Pathogenic/Likely pathogenic only).",
    )
    parser.add_argument("--out_summary_tsv", required=True, help="Output TSV (counts).")
    parser.add_argument(
        "--group_key",
        default="AlleleID",
        help="Grouping key column (default: AlleleID). Alternative: VariationID.",
    )
    return parser.parse_args()


def _norm(value: str) -> str:
    """
    Normalise string values.

    Args:
        value: Input value.

    Returns:
        Stripped string.
    """
    return (value or "").strip()


def _assembly_rank(assembly: str) -> int:
    """
    Rank assembly preference.

    Args:
        assembly: Assembly string.

    Returns:
        Rank integer (lower is preferred).
    """
    a = _norm(assembly)
    if a == "GRCh38":
        return 0
    if a == "GRCh37":
        return 1
    return 2


def _review_rank(review_status: str) -> int:
    """
    Rank review status.

    Args:
        review_status: ReviewStatus field.

    Returns:
        Rank integer (lower is preferred).
    """
    rs = _norm(review_status).lower()
    return REVIEW_STATUS_RANK.get(rs, 50)


def _clinsig_rank(clinsig: str) -> int:
    """
    Rank clinical significance.

    Args:
        clinsig: ClinicalSignificance field.

    Returns:
        Rank integer (lower is preferred).
    """
    cs = _norm(clinsig).lower()
    return CLINSIG_RANK.get(cs, 40)


def score_row(row: Dict[str, str]) -> RowScore:
    """
    Compute a ranking score for a row.

    Args:
        row: TSV row.

    Returns:
        RowScore.
    """
    return RowScore(
        assembly_rank=_assembly_rank(row.get("Assembly", "")),
        review_rank=_review_rank(row.get("ReviewStatus", "")),
        clinsig_rank=_clinsig_rank(row.get("ClinicalSignificance", "")),
    )


def is_pathogenicish(clinsig: str) -> bool:
    """
    Determine if a clinical significance is pathogenic-ish.

    Args:
        clinsig: ClinicalSignificance value.

    Returns:
        True if pathogenic or likely pathogenic (including combined labels).
    """
    cs = _norm(clinsig).lower()
    return cs in {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}


def collapse_to_best(
    in_tsv: Path,
    group_key: str,
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Collapse input rows to the best row per group_key.

    Args:
        in_tsv: Input TSV path.
        group_key: Column to group by.

    Returns:
        A tuple of:
        - header fieldnames
        - dict mapping group value -> best row
    """
    best_rows: Dict[str, Dict[str, str]] = {}
    best_scores: Dict[str, RowScore] = {}

    with in_tsv.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")
        fieldnames = [f.lstrip("#") for f in reader.fieldnames]

        # DictReader uses original keys; we want to support '#AlleleID'
        raw_fieldnames = reader.fieldnames
        norm_map = {f.lstrip("#"): f for f in raw_fieldnames}

        if group_key not in norm_map:
            raise ValueError(f"Grouping key '{group_key}' not present in input header.")

        for row in reader:
            group_val = _norm(row.get(norm_map[group_key], ""))
            if not group_val:
                continue

            # normalise row keys for output consistency
            norm_row: Dict[str, str] = {}
            for norm, raw in norm_map.items():
                norm_row[norm] = row.get(raw, "")

            s = score_row(norm_row)

            if group_val not in best_rows:
                best_rows[group_val] = norm_row
                best_scores[group_val] = s
                continue

            if s < best_scores[group_val]:
                best_rows[group_val] = norm_row
                best_scores[group_val] = s

    return fieldnames, best_rows


def write_tsv(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    """
    Write rows to a TSV.

    Args:
        path: Output path.
        fieldnames: Output columns.
        rows: List of row dicts.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(path: Path, rows: List[Dict[str, str]]) -> None:
    """
    Write summary counts by ClinicalSignificance and ReviewStatus.

    Args:
        path: Output TSV path.
        rows: Rows to summarise.
    """
    clinsig_counts: Counter[str] = Counter()
    review_counts: Counter[str] = Counter()

    for row in rows:
        clinsig_counts[_norm(row.get("ClinicalSignificance", ""))] += 1
        review_counts[_norm(row.get("ReviewStatus", ""))] += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "category", "count"])
        writer.writerow(["rows", "kept", len(rows)])
        for k, v in clinsig_counts.most_common():
            writer.writerow(["ClinicalSignificance", k, v])
        for k, v in review_counts.most_common():
            writer.writerow(["ReviewStatus", k, v])


def main() -> None:
    """
    Entry point.
    """
    args = parse_args()

    in_tsv = Path(args.in_tsv)
    out_best_tsv = Path(args.out_best_tsv)
    out_best_pathogenic_tsv = Path(args.out_best_pathogenic_tsv)
    out_summary_tsv = Path(args.out_summary_tsv)

    if not in_tsv.exists():
        raise FileNotFoundError(f"Missing input: {in_tsv}")

    fieldnames, best_map = collapse_to_best(in_tsv=in_tsv, group_key=args.group_key)
    best_rows = list(best_map.values())

    # Sort for stable outputs (GeneSymbol then VariationID if present)
    def _sort_key(r: Dict[str, str]) -> Tuple[str, str, str]:
        return (_norm(r.get("GeneSymbol", "")), _norm(r.get("VariationID", "")), _norm(r.get("AlleleID", "")))

    best_rows_sorted = sorted(best_rows, key=_sort_key)
    write_tsv(path=out_best_tsv, fieldnames=fieldnames, rows=best_rows_sorted)
    write_summary(path=out_summary_tsv, rows=best_rows_sorted)

    pathogenic_rows = [r for r in best_rows_sorted if is_pathogenicish(r.get("ClinicalSignificance", ""))]
    write_tsv(path=out_best_pathogenic_tsv, fieldnames=fieldnames, rows=pathogenic_rows)


if __name__ == "__main__":
    main()
