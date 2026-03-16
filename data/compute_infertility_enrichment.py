import logging
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


INPUT_SUMMARY_TSV = Path("validation_outputs/infertility_gene_recovery_summary.tsv")
OUTPUT_TSV = Path("infertility_gene_enrichment_statistics.tsv")

BACKGROUND_TOTAL = 13621
REFERENCE_TOTAL = 268

EXCLUDE_GENE_SETS = {
    "HPO_gene_set",
    "ClinVar_best_any",
    "ClinVar_best_pathogenic",
    "ClinVar_high_confidence_any",
    "ClinVar_high_confidence_pathogenic",
    "Literature_gene_set",
}


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """
    Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : list[float]
        Raw p-values.

    Returns
    -------
    list[float]
        Adjusted q-values in original order.
    """
    n_tests = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])

    adjusted = [0.0] * n_tests
    prev_q = 1.0

    for rank, (original_index, p_value) in enumerate(
        reversed(indexed),
        start=1
    ):
        bh_value = p_value * n_tests / (n_tests - rank + 1)
        q_value = min(prev_q, bh_value, 1.0)
        adjusted[original_index] = q_value
        prev_q = q_value

    return adjusted


def load_summary(path: Path) -> pd.DataFrame:
    """
    Load the infertility recovery summary TSV.

    Parameters
    ----------
    path : Path
        Path to the summary TSV.

    Returns
    -------
    pd.DataFrame
        Loaded summary dataframe.
    """
    logging.info("Loading summary table from %s", path)
    df = pd.read_csv(path, sep="\t")
    logging.info("Loaded %d rows", len(df))
    return df


def compute_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute fold-enrichment and Fisher's exact test statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Summary dataframe with gene counts and overlaps.

    Returns
    -------
    pd.DataFrame
        Enrichment results dataframe.
    """
    logging.info("Computing enrichment statistics")

    records = []

    for _, row in df.iterrows():
        gene_set = row["gene_set"]

        if gene_set in EXCLUDE_GENE_SETS:
            logging.info("Skipping reference/source set: %s", gene_set)
            continue

        gene_count = int(row["gene_count"])
        reference_overlap = int(row["reference_overlap"])

        non_overlap_in_set = gene_count - reference_overlap
        overlap_outside_set = REFERENCE_TOTAL - reference_overlap
        non_overlap_outside_set = (
            BACKGROUND_TOTAL
            - gene_count
            - overlap_outside_set
        )

        contingency_table = [
            [reference_overlap, non_overlap_in_set],
            [overlap_outside_set, non_overlap_outside_set],
        ]

        odds_ratio, p_value = fisher_exact(
            contingency_table,
            alternative="greater"
        )

        set_proportion = reference_overlap / gene_count
        background_proportion = REFERENCE_TOTAL / BACKGROUND_TOTAL
        fold_enrichment = set_proportion / background_proportion

        records.append(
            {
                "gene_set": gene_set,
                "gene_count": gene_count,
                "reference_overlap": reference_overlap,
                "set_proportion": round(set_proportion, 6),
                "background_proportion": round(background_proportion, 6),
                "fold_enrichment": round(fold_enrichment, 6),
                "odds_ratio": round(odds_ratio, 6),
                "fisher_p_value": p_value,
            }
        )

    result_df = pd.DataFrame(records)

    q_values = benjamini_hochberg(
        result_df["fisher_p_value"].tolist()
    )
    result_df["fdr_q_value"] = q_values

    return result_df.sort_values(
        by="fold_enrichment",
        ascending=False
    ).reset_index(drop=True)


def write_results(df: pd.DataFrame, path: Path) -> None:
    """
    Write enrichment results to TSV.

    Parameters
    ----------
    df : pd.DataFrame
        Results dataframe.
    path : Path
        Output TSV path.
    """
    logging.info("Writing enrichment results to %s", path)
    df.to_csv(path, sep="\t", index=False)
    logging.info("Done")


def main() -> None:
    summary_df = load_summary(INPUT_SUMMARY_TSV)
    result_df = compute_enrichment(summary_df)
    write_results(result_df, OUTPUT_TSV)


if __name__ == "__main__":
    main()
