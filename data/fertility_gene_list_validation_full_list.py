import logging
from pathlib import Path

import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


INPUT_EXCEL = Path("SUMMARY_fertility_evidence_gene_lists_only.xlsx")
REFERENCE_GENES = Path("curated_male_infertility_reference_genes.tsv")
OUTPUT_DIR = Path("validation_outputs")

OUTPUT_DIR.mkdir(exist_ok=True)


def load_reference_genes(path: Path) -> set:
    logging.info("Loading infertility reference gene list")

    df = pd.read_csv(path, sep="\t")
    genes = set(df["gene"].dropna().str.upper())

    logging.info("Loaded %d reference genes", len(genes))
    return genes


def load_gene_sets(path: Path) -> dict:
    logging.info("Loading gene sets from Excel")

    df = pd.read_excel(path, sheet_name=0)

    gene_sets = {}

    for column in df.columns:
        genes = set(
            df[column]
            .dropna()
            .astype(str)
            .str.strip()
            .str.upper()
        )

        gene_sets[column] = genes
        logging.info(
            "Loaded %d genes for set %s",
            len(genes),
            column
        )

    return gene_sets


def compute_overlap(reference: set, gene_sets: dict) -> pd.DataFrame:
    logging.info("Computing overlap with infertility reference genes")

    records = []

    for name, genes in gene_sets.items():
        overlap = genes.intersection(reference)

        records.append(
            {
                "gene_set": name,
                "gene_count": len(genes),
                "reference_overlap": len(overlap),
                "percent_overlap": round(
                    (len(overlap) / len(genes)) * 100,
                    3
                ),
            }
        )

    return pd.DataFrame(records)


def write_overlap_genes(
    reference: set,
    gene_sets: dict
) -> None:

    logging.info("Writing overlapping gene lists")

    for name, genes in gene_sets.items():
        overlap = sorted(genes.intersection(reference))

        output_file = OUTPUT_DIR / f"{name}_validated_overlap.tsv"

        pd.DataFrame(
            {"gene": overlap}
        ).to_csv(
            output_file,
            sep="\t",
            index=False
        )


def main():

    reference_genes = load_reference_genes(REFERENCE_GENES)

    gene_sets = load_gene_sets(INPUT_EXCEL)

    summary_df = compute_overlap(reference_genes, gene_sets)

    summary_file = OUTPUT_DIR / "infertility_gene_recovery_summary.tsv"

    summary_df.to_csv(
        summary_file,
        sep="\t",
        index=False
    )

    logging.info("Summary written to %s", summary_file)

    write_overlap_genes(reference_genes, gene_sets)


if __name__ == "__main__":
    main()
