#!/usr/bin/env python3
"""
make_gene_list_overlap_plots.py

Read gene lists from sheet 1 of an Excel workbook (one gene list per column)
and generate Venn diagrams (2- and 3-way) plus an UpSet-style intersection plot.

Usage
-----
python make_gene_list_overlap_plots.py \
    --input_excel gene_lists.xlsx \
    --output_dir overlap_plots

Optional:
    --sheet_name Sheet1
    --min_set_size 1
    --max_sets_for_upset 20
"""

from __future__ import annotations

import argparse
import itertools
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2, venn3


def setup_logging(verbose: bool) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate Venn diagrams and an UpSet-style plot "
        "from Excel gene lists stored in columns."
    )
    parser.add_argument(
        "--input_excel",
        required=True,
        help="Path to the Excel file.",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where plots will be written.",
    )
    parser.add_argument(
        "--sheet_name",
        default=0,
        help=(
            "Sheet name or index to read. "
            "Default is 0 (first sheet)."
        ),
    )
    parser.add_argument(
        "--min_set_size",
        type=int,
        default=1,
        help="Drop columns with fewer than this many genes. Default: 1",
    )
    parser.add_argument(
        "--max_sets_for_upset",
        type=int,
        default=20,
        help=(
            "Maximum number of sets to include in the UpSet-style plot. "
            "Default: 20"
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def normalise_gene(value: object) -> str:
    """Normalise a gene symbol / ID."""
    return str(value).strip()


def load_gene_lists(
    input_excel: str,
    sheet_name: object,
    min_set_size: int,
) -> Dict[str, Set[str]]:
    """Load gene lists from an Excel sheet, one list per column."""
    logging.info("Reading Excel file: %s", input_excel)
    df = pd.read_excel(input_excel, sheet_name=sheet_name)

    gene_lists: Dict[str, Set[str]] = {}

    for column in df.columns:
        series = df[column].dropna()
        genes = {
            normalise_gene(value)
            for value in series
            if normalise_gene(value) != ""
        }

        if len(genes) >= min_set_size:
            gene_lists[str(column)] = genes
            logging.info(
                "Loaded %d genes from column '%s'",
                len(genes),
                column,
            )
        else:
            logging.info(
                "Skipping column '%s' because it has only %d genes",
                column,
                len(genes),
            )

    if not gene_lists:
        raise ValueError("No usable gene lists were found.")

    return gene_lists


def safe_filename(text: str) -> str:
    """Convert text to a filesystem-safe filename fragment."""
    return "".join(
        char if char.isalnum() or char in {"-", "_"} else "_"
        for char in text
    )


def save_venn2(
    set_a_name: str,
    set_b_name: str,
    set_a: Set[str],
    set_b: Set[str],
    output_dir: Path,
) -> None:
    """Create and save a 2-set Venn diagram."""
    plt.figure(figsize=(7, 7))
    venn2(
        [set_a, set_b],
        set_labels=(set_a_name, set_b_name),
    )
    plt.title(f"Venn diagram: {set_a_name} vs {set_b_name}")
    plt.tight_layout()

    filename = (
        f"venn2_{safe_filename(set_a_name)}_"
        f"{safe_filename(set_b_name)}.pdf"
    )
    out_path = output_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info("Wrote %s", out_path)


def save_venn3(
    names: Sequence[str],
    sets: Sequence[Set[str]],
    output_dir: Path,
) -> None:
    """Create and save a 3-set Venn diagram."""
    plt.figure(figsize=(8, 8))
    venn3(
        [sets[0], sets[1], sets[2]],
        set_labels=(names[0], names[1], names[2]),
    )
    plt.title(f"Venn diagram: {' vs '.join(names)}")
    plt.tight_layout()

    filename = (
        f"venn3_{safe_filename(names[0])}_"
        f"{safe_filename(names[1])}_"
        f"{safe_filename(names[2])}.pdf"
    )
    out_path = output_dir / filename
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info("Wrote %s", out_path)


def all_nonempty_intersections(
    gene_lists: Dict[str, Set[str]],
) -> List[Tuple[Tuple[str, ...], int]]:
    """
    Compute all non-empty exact intersections across input sets.

    Returns a list of:
        ((set_name_1, set_name_2, ...), size)
    where the tuple indicates the exact combination of sets
    in which the genes are present.
    """
    set_names = list(gene_lists.keys())
    all_genes = sorted(set().union(*gene_lists.values()))

    pattern_counts: Dict[Tuple[str, ...], int] = {}

    for gene in all_genes:
        present_in = tuple(
            name for name in set_names if gene in gene_lists[name]
        )
        if present_in:
            pattern_counts[present_in] = (
                pattern_counts.get(present_in, 0) + 1
            )

    intersections = sorted(
        pattern_counts.items(),
        key=lambda item: item[1],
        reverse=True,
    )
    return intersections


def make_upset_style_plot(
    gene_lists: Dict[str, Set[str]],
    output_dir: Path,
    max_sets_for_upset: int,
) -> None:
    """
    Create a simple UpSet-style plot using matplotlib only.

    This avoids relying on the upsetplot package.
    """
    if len(gene_lists) > max_sets_for_upset:
        logging.warning(
            "Skipping UpSet-style plot because there are %d sets and "
            "--max_sets_for_upset is %d. "
            "Reduce the number of columns or increase the limit.",
            len(gene_lists),
            max_sets_for_upset,
        )
        return

    set_names = list(gene_lists.keys())
    intersections = all_nonempty_intersections(gene_lists)

    if not intersections:
        logging.warning("No intersections found for UpSet-style plot.")
        return

    # Keep top intersections to avoid an unreadable plot.
    max_intersections_to_plot = 25
    intersections = intersections[:max_intersections_to_plot]

    labels = [
        "\n".join(combo) for combo, _ in intersections
    ]
    counts = [count for _, count in intersections]

    fig = plt.figure(figsize=(max(10, len(intersections) * 0.5), 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 2], hspace=0.05)

    # Top panel: bar plot of intersection sizes
    ax_bar = fig.add_subplot(gs[0])
    x_positions = list(range(len(intersections)))
    ax_bar.bar(x_positions, counts)
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_xticks([])
    ax_bar.set_title("UpSet-style plot of exact set intersections")

    # Bottom panel: matrix of set membership
    ax_mat = fig.add_subplot(gs[1], sharex=ax_bar)

    for row_index, set_name in enumerate(reversed(set_names)):
        for col_index, (combo, _) in enumerate(intersections):
            if set_name in combo:
                ax_mat.plot(col_index, row_index, "o")
        # connect dots in each column
    for col_index, (combo, _) in enumerate(intersections):
        rows = [
            idx
            for idx, set_name in enumerate(reversed(set_names))
            if set_name in combo
        ]
        if len(rows) >= 2:
            ax_mat.plot(
                [col_index] * len(rows),
                rows,
                "-",
                linewidth=1,
            )

    ax_mat.set_yticks(range(len(set_names)))
    ax_mat.set_yticklabels(list(reversed(set_names)))
    ax_mat.set_xlabel("Exact intersections (top 25 by size)")
    ax_mat.set_xlim(-0.5, len(intersections) - 0.5)

    # Put shortened labels under the matrix if desired.
    ax_mat.set_xticks(x_positions)
    ax_mat.set_xticklabels(
        [str(i + 1) for i in x_positions],
        rotation=0,
    )

    # Optional legend text mapping index -> combination
    mapping_lines = [
        f"{i + 1}: {' & '.join(combo)} ({count})"
        for i, (combo, count) in enumerate(intersections)
    ]
    fig.text(
        0.01,
        0.01,
        "\n".join(mapping_lines),
        fontsize=8,
        va="bottom",
        ha="left",
    )

    out_path = output_dir / "upset_style_plot.pdf"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info("Wrote %s", out_path)


def write_set_sizes(
    gene_lists: Dict[str, Set[str]],
    output_dir: Path,
) -> None:
    """Write a TSV of set sizes."""
    rows = [
        {"set_name": name, "gene_count": len(genes)}
        for name, genes in gene_lists.items()
    ]
    df = pd.DataFrame(rows).sort_values(
        by="gene_count",
        ascending=False,
    )
    out_path = output_dir / "set_sizes.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info("Wrote %s", out_path)


def write_intersections(
    gene_lists: Dict[str, Set[str]],
    output_dir: Path,
) -> None:
    """Write exact intersections as a TSV."""
    intersections = all_nonempty_intersections(gene_lists)
    rows = []
    for combo, count in intersections:
        rows.append(
            {
                "sets": " & ".join(combo),
                "n_sets": len(combo),
                "intersection_size": count,
            }
        )

    df = pd.DataFrame(rows).sort_values(
        by="intersection_size",
        ascending=False,
    )
    out_path = output_dir / "exact_intersections.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logging.info("Wrote %s", out_path)


def main() -> None:
    """Run the script."""
    args = parse_args()
    setup_logging(args.verbose)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_lists = load_gene_lists(
        input_excel=args.input_excel,
        sheet_name=args.sheet_name,
        min_set_size=args.min_set_size,
    )

    write_set_sizes(gene_lists, output_dir)
    write_intersections(gene_lists, output_dir)

    # 2-way Venn diagrams
    for name_a, name_b in itertools.combinations(gene_lists.keys(), 2):
        save_venn2(
            name_a,
            name_b,
            gene_lists[name_a],
            gene_lists[name_b],
            output_dir,
        )

    # 3-way Venn diagrams
    for combo_names in itertools.combinations(gene_lists.keys(), 3):
        combo_sets = [gene_lists[name] for name in combo_names]
        save_venn3(combo_names, combo_sets, output_dir)

    # UpSet-style plot
    make_upset_style_plot(
        gene_lists=gene_lists,
        output_dir=output_dir,
        max_sets_for_upset=args.max_sets_for_upset,
    )

    logging.info("All done.")


if __name__ == "__main__":
    main()
