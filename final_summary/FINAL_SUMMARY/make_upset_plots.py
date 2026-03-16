#!/usr/bin/env python3
"""
make_clean_upset.py

Create a clean publication-style UpSet plot from an Excel sheet where each
column contains a gene list.

Features
--------
- Reads sheet 1 (or a user-specified sheet) from an Excel file
- Treats each column as one gene list
- Allows selection of a subset of columns
- Removes blanks and duplicate genes within each list
- Writes PDF and SVG output
- Uses upsetplot for cleaner publication-quality visuals

Example
-------
python make_clean_upset.py \
    --input_excel gene_lists.xlsx \
    --output_prefix figures/upset_main \
    --columns "HPO infertility genes" "Testis-enriched genes" "Sperm RNA" \
              "Combined proteome" "Prioritised candidates" "Druggable targets" \
    --min_degree 1 \
    --max_intersections 15
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Set

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_contents


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
        description="Create a clean UpSet plot from Excel column-based gene lists."
    )
    parser.add_argument(
        "--input_excel",
        required=True,
        help="Path to input Excel file.",
    )
    parser.add_argument(
        "--sheet_name",
        default=0,
        help="Sheet name or index. Default: 0 (first sheet).",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Output prefix (without extension). PDF and SVG will be written.",
    )
    parser.add_argument(
        "--columns",
        nargs="+",
        required=True,
        help="Exact column names to include in the UpSet plot.",
    )
    parser.add_argument(
        "--max_intersections",
        type=int,
        default=15,
        help="Maximum number of intersections to display. Default: 15",
    )
    parser.add_argument(
        "--min_degree",
        type=int,
        default=1,
        help="Minimum intersection degree to display. Default: 1",
    )
    parser.add_argument(
        "--fig_width",
        type=float,
        default=10.0,
        help="Figure width in inches. Default: 10",
    )
    parser.add_argument(
        "--fig_height",
        type=float,
        default=6.5,
        help="Figure height in inches. Default: 6.5",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def normalise_gene(value: object) -> str:
    """Normalise gene symbols/IDs."""
    return str(value).strip()


def load_selected_gene_lists(
    input_excel: str,
    sheet_name: object,
    selected_columns: List[str],
) -> Dict[str, Set[str]]:
    """
    Load selected gene-list columns from the Excel sheet.

    Each selected column becomes one set of genes.
    """
    logging.info("Reading Excel file: %s", input_excel)
    df = pd.read_excel(input_excel, sheet_name=sheet_name)

    missing = [col for col in selected_columns if col not in df.columns]
    if missing:
        raise ValueError(
            "The following requested columns were not found in the Excel file: "
            + ", ".join(missing)
        )

    gene_lists: Dict[str, Set[str]] = {}

    for column in selected_columns:
        series = df[column].dropna()
        genes = {
            normalise_gene(value)
            for value in series
            if normalise_gene(value) != ""
        }
        gene_lists[column] = genes
        logging.info("Loaded %d genes from '%s'", len(genes), column)

    return gene_lists


def write_summary_tables(
    gene_lists: Dict[str, Set[str]],
    output_prefix: str,
) -> None:
    """Write simple TSV summaries of set sizes."""
    output_path = Path(output_prefix)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = [
        {"set_name": name, "gene_count": len(genes)}
        for name, genes in gene_lists.items()
    ]
    df = pd.DataFrame(rows).sort_values(by="gene_count", ascending=False)
    out_tsv = output_dir / f"{output_path.name}_set_sizes.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    logging.info("Wrote %s", out_tsv)


def make_upset_plot(
    gene_lists: Dict[str, Set[str]],
    output_prefix: str,
    max_intersections: int,
    min_degree: int,
    fig_width: float,
    fig_height: float,
) -> None:
    """Create and save a clean UpSet plot."""
    output_path = Path(output_prefix)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    contents = {
        set_name: sorted(genes)
        for set_name, genes in gene_lists.items()
    }

    upset_data = from_contents(contents)

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
        }
    )

    fig = plt.figure(figsize=(fig_width, fig_height))

    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by=None,
        min_degree=min_degree,
        max_subset_rank=max_intersections,
        intersection_plot_elements=max_intersections,
    )

    upset.plot(fig=fig)

    # Clean title handling
    fig.suptitle(
        "Overlap of prioritised sperm-associated gene sets",
        y=0.98,
        fontsize=12,
    )

    pdf_path = output_dir / f"{output_path.name}.pdf"
    svg_path = output_dir / f"{output_path.name}.svg"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    logging.info("Wrote %s", pdf_path)
    logging.info("Wrote %s", svg_path)


def main() -> None:
    """Entry point."""
    args = parse_args()
    setup_logging(args.verbose)

    gene_lists = load_selected_gene_lists(
        input_excel=args.input_excel,
        sheet_name=args.sheet_name,
        selected_columns=args.columns,
    )

    write_summary_tables(
        gene_lists=gene_lists,
        output_prefix=args.output_prefix,
    )

    make_upset_plot(
        gene_lists=gene_lists,
        output_prefix=args.output_prefix,
        max_intersections=args.max_intersections,
        min_degree=args.min_degree,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )

    logging.info("Done.")


if __name__ == "__main__":
    main()