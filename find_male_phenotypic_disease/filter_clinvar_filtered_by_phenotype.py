#!/usr/bin/env python3
"""
Filter a ClinVar variant_summary subset to variants whose phenotypes indicate male infertility.

This script expects a TSV produced from ClinVar variant_summary filtering by gene list.
It then filters rows based on matching keywords against the PhenotypeList column.

Outputs are TSV (tab-separated).

Typical workflow:
1) Filter ClinVar variant_summary.txt.gz by gene list -> large TSV
2) Filter that TSV by phenotype terms -> male infertility relevant variants
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_PHENOTYPE_KEYWORDS: Tuple[str, ...] = (
    "male infertility",
    "infertility, male",
    "infertility",
    "azoospermia",
    "non-obstructive azoospermia",
    "obstructive azoospermia",
    "oligozoospermia",
    "severe oligozoospermia",
    "cryptozoospermia",
    "asthenozoospermia",
    "asthenospermia",
    "teratozoospermia",
    "teratospermia",
    "oligoasthenoteratozoospermia",
    "globozoospermia",
    "macrozoospermia",
    "spermatogenic failure",
    "spermatogenesis maturation arrest",
    "maturation arrest",
    "hypospermatogenesis",
    "sertoli cell-only",
    "spermiation failure",
    "absent vas deferens",
    "congenital bilateral absence of the vas deferens",
    "cbavd",
    "retrograde ejaculation",
    "hypospermia",
    "low semen volume",
    "aspermia",
    "immotile sperm",
    "abnormal sperm motility",
    "reduced sperm motility",
    "reduced progressive sperm motility",
    "abnormal sperm morphology",
    "abnormal sperm head morphology",
    "abnormal sperm tail morphology",
    "absent sperm flagella",
    "short sperm flagella",
    "coiled sperm flagella",
    "sperm flagella",
    "axoneme",
    "acrosome",
    "acrosomal hypoplasia",
)

EXCLUDE_KEYWORDS: Tuple[str, ...] = (
    # To avoid pulling female-only fertility contexts when we match generic "infertility"
    "female infertility",
    "oocyte",
    "ovary",
    "uterus",
    "uterine",
    "endometr",
    "ovulation",
    "fallopian",
    "cervix",
    "vagina",
    "vaginal",
)


def parse_args() -> argparse.Namespace:
    """
    Parse arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(description="Filter ClinVar TSV by male infertility phenotype keywords.")
    parser.add_argument(
        "--in_tsv",
        required=True,
        help="Input TSV (ClinVar variants already filtered by gene list).",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV containing only rows matching phenotype keywords.",
    )
    parser.add_argument(
        "--out_minimal_tsv",
        required=True,
        help="Output TSV with a minimal set of columns for reporting.",
    )
    parser.add_argument(
        "--out_summary_tsv",
        required=True,
        help="Output TSV with counts (rows kept, ClinicalSignificance distribution, etc.).",
    )
    parser.add_argument(
        "--keyword_file",
        default="",
        help="Optional text file with one phenotype keyword per line (overrides defaults).",
    )
    parser.add_argument(
        "--exclude_keyword_file",
        default="",
        help="Optional text file with one exclude keyword per line (overrides defaults).",
    )
    parser.add_argument(
        "--case_sensitive",
        action="store_true",
        help="If set, use case-sensitive matching (default is case-insensitive).",
    )
    parser.add_argument(
        "--require_clinsig_substring",
        default="",
        help="Optional substring filter on ClinicalSignificance (e.g. 'Pathogenic').",
    )
    return parser.parse_args()


def load_keywords(path_str: str, defaults: Sequence[str]) -> List[str]:
    """
    Load keywords from file or return defaults.

    Args:
        path_str: Path to keyword file (one per line) or empty string.
        defaults: Default keywords.

    Returns:
        List of keywords.
    """
    if not path_str.strip():
        return list(defaults)

    path = Path(path_str)
    keywords: List[str] = []
    with path.open(mode="r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                keywords.append(line)
    return keywords


def compile_pattern(keywords: List[str], case_sensitive: bool) -> re.Pattern:
    """
    Compile a regex pattern that matches any keyword (escaped) as substring.

    Args:
        keywords: Keywords to match.
        case_sensitive: Whether matching is case-sensitive.

    Returns:
        Compiled regex pattern.
    """
    flags = 0 if case_sensitive else re.IGNORECASE
    escaped = [re.escape(k) for k in keywords if k.strip()]
    if not escaped:
        raise ValueError("No keywords provided for matching.")
    return re.compile("|".join(escaped), flags=flags)


def normalise_header(fieldnames: List[str]) -> List[str]:
    """
    Normalise headers so '#AlleleID' becomes 'AlleleID' for consistent lookup.

    Args:
        fieldnames: Raw fieldnames.

    Returns:
        Normalised fieldnames.
    """
    normed: List[str] = []
    for f in fieldnames:
        if f.startswith("#"):
            normed.append(f.lstrip("#"))
        else:
            normed.append(f)
    return normed


def filter_rows(
    in_tsv: Path,
    out_tsv: Path,
    out_minimal_tsv: Path,
    out_summary_tsv: Path,
    include_pat: re.Pattern,
    exclude_pat: re.Pattern,
    require_clinsig_substring: str,
    case_sensitive: bool,
) -> None:
    """
    Filter ClinVar TSV rows by phenotype keywords.

    Args:
        in_tsv: Input TSV path.
        out_tsv: Output TSV path (full columns).
        out_minimal_tsv: Output TSV path (minimal columns).
        out_summary_tsv: Output TSV path (counts).
        include_pat: Compiled include regex.
        exclude_pat: Compiled exclude regex.
        require_clinsig_substring: Optional substring filter for ClinicalSignificance.
        case_sensitive: Case sensitivity flag.
    """
    kept = 0
    total = 0

    clinsig_counts: Counter[str] = Counter()
    review_counts: Counter[str] = Counter()
    keyword_hit_counts: Counter[str] = Counter()

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_minimal_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    require_cs = require_clinsig_substring.strip()
    require_cs_cmp = require_cs if case_sensitive else require_cs.lower()

    with in_tsv.open(mode="r", encoding="utf-8", newline="") as handle_in:
        reader = csv.DictReader(handle_in, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")

        raw_fields = reader.fieldnames
        fields = normalise_header(raw_fields)

        # Build a mapping from normalised field name -> raw field name (for DictReader row keys)
        norm_to_raw: Dict[str, str] = {}
        for raw, norm in zip(raw_fields, fields):
            norm_to_raw[norm] = raw

        if "PhenotypeList" not in norm_to_raw:
            raise ValueError("Input TSV does not contain a PhenotypeList column.")
        if "ClinicalSignificance" not in norm_to_raw:
            raise ValueError("Input TSV does not contain a ClinicalSignificance column.")
        if "ReviewStatus" not in norm_to_raw:
            raise ValueError("Input TSV does not contain a ReviewStatus column.")

        pheno_key = norm_to_raw["PhenotypeList"]
        clinsig_key = norm_to_raw["ClinicalSignificance"]
        review_key = norm_to_raw["ReviewStatus"]

        # Minimal columns for the meeting-ready output
        minimal_fields = [
            "AlleleID",
            "VariationID",
            "GeneSymbol",
            "ClinicalSignificance",
            "ReviewStatus",
            "LastEvaluated",
            "RS# (dbSNP)",
            "RCVaccession",
            "PhenotypeList",
            "Assembly",
            "Chromosome",
            "Start",
            "Stop",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
            "Type",
            "Name",
        ]

        # Map minimal fieldnames to raw keys where possible; keep if present
        minimal_raw_keys: List[str] = []
        minimal_out_headers: List[str] = []
        for col in minimal_fields:
            raw_key = norm_to_raw.get(col, None)
            if raw_key is not None:
                minimal_raw_keys.append(raw_key)
                minimal_out_headers.append(col)

        with out_tsv.open(mode="w", encoding="utf-8", newline="") as handle_out_full, \
                out_minimal_tsv.open(mode="w", encoding="utf-8", newline="") as handle_out_min:

            writer_full = csv.DictWriter(handle_out_full, delimiter="\t", fieldnames=raw_fields, extrasaction="ignore")
            writer_full.writeheader()

            writer_min = csv.writer(handle_out_min, delimiter="\t")
            writer_min.writerow(minimal_out_headers)

            for row in reader:
                total += 1

                phenotype_list = (row.get(pheno_key) or "").strip()
                if not phenotype_list:
                    continue

                if exclude_pat.search(phenotype_list):
                    continue

                if not include_pat.search(phenotype_list):
                    continue

                clinsig = (row.get(clinsig_key) or "").strip()
                if require_cs:
                    clinsig_cmp = clinsig if case_sensitive else clinsig.lower()
                    if require_cs_cmp not in clinsig_cmp:
                        continue

                writer_full.writerow(row)

                min_row = [(row.get(k) or "").strip() for k in minimal_raw_keys]
                writer_min.writerow(min_row)

                kept += 1
                clinsig_counts[clinsig] += 1
                review_counts[(row.get(review_key) or "").strip()] += 1

                # For rough keyword hit counting, count which include keywords appear
                # (approximate; useful for auditing)
                for kw in include_pat.pattern.split("|"):
                    pass  # pattern is escaped; skip per-keyword counting in this compact version

    with out_summary_tsv.open(mode="w", encoding="utf-8", newline="") as handle_sum:
        writer = csv.writer(handle_sum, delimiter="\t")
        writer.writerow(["metric", "category", "count"])
        writer.writerow(["rows", "total", total])
        writer.writerow(["rows", "kept", kept])

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
    out_tsv = Path(args.out_tsv)
    out_minimal_tsv = Path(args.out_minimal_tsv)
    out_summary_tsv = Path(args.out_summary_tsv)

    if not in_tsv.exists():
        raise FileNotFoundError(f"Missing input TSV: {in_tsv}")

    include_keywords = load_keywords(args.keyword_file, DEFAULT_PHENOTYPE_KEYWORDS)
    exclude_keywords = load_keywords(args.exclude_keyword_file, EXCLUDE_KEYWORDS)

    include_pat = compile_pattern(include_keywords, case_sensitive=bool(args.case_sensitive))
    exclude_pat = compile_pattern(exclude_keywords, case_sensitive=bool(args.case_sensitive))

    filter_rows(
        in_tsv=in_tsv,
        out_tsv=out_tsv,
        out_minimal_tsv=out_minimal_tsv,
        out_summary_tsv=out_summary_tsv,
        include_pat=include_pat,
        exclude_pat=exclude_pat,
        require_clinsig_substring=args.require_clinsig_substring,
        case_sensitive=bool(args.case_sensitive),
    )


if __name__ == "__main__":
    main()
