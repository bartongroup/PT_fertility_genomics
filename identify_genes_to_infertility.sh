#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# identify_genes_to_infertility.sh
#
# End-to-end pipeline:
#   HPO expansion -> ClinVar filtering -> Evidence ranking -> Overlap
#
# All outputs written to a single master directory: results/
###############################################################################

###############################################################################
# Base directories
###############################################################################

BASE_DIR="${HOME}/data/2026_sperm_Gates"
DATA_DIR="${BASE_DIR}/data"
SCRIPTS_DIR="${BASE_DIR}/scripts/find_male_phenotypic_disease"
RESULTS_DIR="${BASE_DIR}/results"

###############################################################################
# Input data (from download_data.sh)
###############################################################################

HPO_DIR="${DATA_DIR}/hpo"
CLINVAR_DIR="${DATA_DIR}/clinvar"
CLINVAR_VARIANT_SUMMARY_GZ="${CLINVAR_DIR}/variant_summary.txt.gz"

# External gene list for overlap (edit as needed)
TESTIS_GENE_TSV="${DATA_DIR}/gtex/gtex_v11_testis_gene_symbols.tsv"

###############################################################################
# Output directories
###############################################################################

HPO_OUT="${RESULTS_DIR}/01_hpo"
CLINVAR_GENE_OUT="${RESULTS_DIR}/02_clinvar_gene_filter"
CLINVAR_PHENO_OUT="${RESULTS_DIR}/03_clinvar_phenotype_filter"
CLINVAR_COLLAPSE_OUT="${RESULTS_DIR}/04_clinvar_collapsed"
CLINVAR_HC_OUT="${RESULTS_DIR}/05_clinvar_high_confidence"
REPORT_OUT="${RESULTS_DIR}/06_reports"
OVERLAP_OUT="${RESULTS_DIR}/07_overlap"

mkdir -p "${HPO_OUT}"
mkdir -p "${CLINVAR_GENE_OUT}"
mkdir -p "${CLINVAR_PHENO_OUT}"
mkdir -p "${CLINVAR_COLLAPSE_OUT}"
mkdir -p "${CLINVAR_HC_OUT}"
mkdir -p "${REPORT_OUT}"
mkdir -p "${OVERLAP_OUT}"

###############################################################################
# 1) HPO expansion
###############################################################################

python "${SCRIPTS_DIR}/extract_male_infertility_genes_hpo_ontology.py" \
  --hpo_dir "${HPO_DIR}" \
  --out_dir "${HPO_OUT}" \
  --include_diseases \
  --require_male \
  --use_phrase_seeds

python "${SCRIPTS_DIR}/extract_gene_list_from_hpo_outputs.py" \
  --genes_summary_tsv "${HPO_OUT}/male_infertility_genes_summary.tsv" \
  --out_dir "${HPO_OUT}/gene_lists"

HPO_GENE_SYMBOLS="${HPO_OUT}/gene_lists/male_infertility_gene_symbols.tsv"

###############################################################################
# 2) ClinVar: filter by HPO gene list
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_variant_summary_by_gene_list.py" \
  --variant_summary_gz "${CLINVAR_VARIANT_SUMMARY_GZ}" \
  --gene_symbols_tsv "${HPO_GENE_SYMBOLS}" \
  --out_tsv "${CLINVAR_GENE_OUT}/clinvar_gene_filtered.tsv" \
  --out_summary_tsv "${CLINVAR_GENE_OUT}/clinvar_gene_filtered_counts.tsv"

###############################################################################
# 3) ClinVar: phenotype filter
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_filtered_by_phenotype.py" \
  --in_tsv "${CLINVAR_GENE_OUT}/clinvar_gene_filtered.tsv" \
  --out_tsv "${CLINVAR_PHENO_OUT}/clinvar_infertility_variants.tsv" \
  --out_minimal_tsv "${CLINVAR_PHENO_OUT}/clinvar_infertility_variants_minimal.tsv" \
  --out_summary_tsv "${CLINVAR_PHENO_OUT}/clinvar_infertility_counts.tsv"

###############################################################################
# 4) Collapse to best per AlleleID
###############################################################################

python "${SCRIPTS_DIR}/collapse_clinvar_variants_to_best.py" \
  --in_tsv "${CLINVAR_PHENO_OUT}/clinvar_infertility_variants.tsv" \
  --out_best_tsv "${CLINVAR_COLLAPSE_OUT}/clinvar_best.tsv" \
  --out_best_pathogenic_tsv "${CLINVAR_COLLAPSE_OUT}/clinvar_best_pathogenic.tsv" \
  --out_summary_tsv "${CLINVAR_COLLAPSE_OUT}/clinvar_best_counts.tsv" \
  --group_key AlleleID

###############################################################################
# 5) High-confidence review status filter
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_best_high_confidence.py" \
  --in_tsv "${CLINVAR_COLLAPSE_OUT}/clinvar_best.tsv" \
  --out_high_conf_tsv "${CLINVAR_HC_OUT}/clinvar_high_confidence.tsv" \
  --out_high_conf_pathogenic_tsv "${CLINVAR_HC_OUT}/clinvar_high_confidence_pathogenic.tsv" \
  --out_counts_tsv "${CLINVAR_HC_OUT}/clinvar_high_confidence_counts.tsv"

###############################################################################
# 6) High-confidence pathogenic report
###############################################################################

python "${SCRIPTS_DIR}/make_clinvar_high_confidence_report.py" \
  --in_tsv "${CLINVAR_HC_OUT}/clinvar_high_confidence_pathogenic.tsv" \
  --out_report_tsv "${REPORT_OUT}/clinvar_pathogenic_high_confidence_report.tsv" \
  --out_gene_summary_tsv "${REPORT_OUT}/clinvar_pathogenic_high_confidence_gene_summary.tsv"

###############################################################################
# 7) Overlap with external gene list
###############################################################################

python "${SCRIPTS_DIR}/overlap_testis_specific_with_clinvar_tiers.py" \
  --testis_gene_tsv "${TESTIS_GENE_TSV}" \
  --clinvar_best_tsv "${CLINVAR_COLLAPSE_OUT}/clinvar_best.tsv" \
  --clinvar_hc_pathogenic_tsv "${CLINVAR_HC_OUT}/clinvar_high_confidence_pathogenic.tsv" \
  --out_dir "${OVERLAP_OUT}"

echo "Pipeline completed successfully."
echo "Results written to: ${RESULTS_DIR}"
