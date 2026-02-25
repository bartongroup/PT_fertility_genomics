#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# identify_genes_to_infertility_results.sh
#
# Inputs are taken from the existing folder layout under:
#   ~/data/2026_sperm_Gates/
#
# Outputs are written to:
#   ~/data/2026_sperm_Gates/results/
#
# No functions. Linear commands only.
###############################################################################

BASE_DIR="${HOME}/data/2026_sperm_Gates"

SCRIPTS_DIR="${BASE_DIR}/scripts/find_male_phenotypic_disease"

HPO_DIR="${BASE_DIR}/hpo_data"
CLINVAR_DIR="${BASE_DIR}/clinvar"
CLINVAR_VARIANT_SUMMARY_GZ="${CLINVAR_DIR}/variant_summary.txt.gz"

RESULTS_DIR="${BASE_DIR}/results"

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
# Basic input checks (fail fast)
###############################################################################

test -s "${HPO_DIR}/phenotype.hpoa"
test -s "${HPO_DIR}/genes_to_phenotype.txt"
test -s "${HPO_DIR}/hp.obo" || test -s "${HPO_DIR}/hp.json"
test -s "${CLINVAR_VARIANT_SUMMARY_GZ}"

###############################################################################
# 1) HPO expansion (phrase seeds)
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

HPO_GENE_SYMBOLS_TSV="${HPO_OUT}/gene_lists/male_infertility_gene_symbols.tsv"
test -s "${HPO_GENE_SYMBOLS_TSV}"

###############################################################################
# 2) ClinVar: filter variant_summary to HPO genes
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_variant_summary_by_gene_list.py" \
  --variant_summary_gz "${CLINVAR_VARIANT_SUMMARY_GZ}" \
  --gene_symbols_tsv "${HPO_GENE_SYMBOLS_TSV}" \
  --out_tsv "${CLINVAR_GENE_OUT}/male_infertility_clinvar_variant_summary.tsv" \
  --out_summary_tsv "${CLINVAR_GENE_OUT}/male_infertility_clinvar_variant_summary_counts.tsv"

test -s "${CLINVAR_GENE_OUT}/male_infertility_clinvar_variant_summary.tsv"

###############################################################################
# 3) ClinVar: phenotype filter (infertility-related)
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_filtered_by_phenotype.py" \
  --in_tsv "${CLINVAR_GENE_OUT}/male_infertility_clinvar_variant_summary.tsv" \
  --out_tsv "${CLINVAR_PHENO_OUT}/male_infertility_clinvar_variants_by_phenotype.tsv" \
  --out_minimal_tsv "${CLINVAR_PHENO_OUT}/male_infertility_clinvar_variants_by_phenotype_minimal.tsv" \
  --out_summary_tsv "${CLINVAR_PHENO_OUT}/male_infertility_clinvar_variants_by_phenotype_counts.tsv"

test -s "${CLINVAR_PHENO_OUT}/male_infertility_clinvar_variants_by_phenotype.tsv"

###############################################################################
# 4) Collapse: best row per AlleleID
###############################################################################

python "${SCRIPTS_DIR}/collapse_clinvar_variants_to_best.py" \
  --in_tsv "${CLINVAR_PHENO_OUT}/male_infertility_clinvar_variants_by_phenotype.tsv" \
  --out_best_tsv "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best.tsv" \
  --out_best_pathogenic_tsv "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best_pathogenic.tsv" \
  --out_summary_tsv "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best_counts.tsv" \
  --group_key AlleleID

test -s "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best.tsv"

###############################################################################
# 5) High-confidence: review-status filter
###############################################################################

python "${SCRIPTS_DIR}/filter_clinvar_best_high_confidence.py" \
  --in_tsv "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best.tsv" \
  --out_high_conf_tsv "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence.tsv" \
  --out_high_conf_pathogenic_tsv "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_counts_tsv "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence_counts.tsv"

test -s "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence.tsv"

###############################################################################
# 6) Reports: high-confidence pathogenic report + gene summary
###############################################################################

python "${SCRIPTS_DIR}/make_clinvar_high_confidence_report.py" \
  --in_tsv "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_report_tsv "${REPORT_OUT}/male_infertility_clinvar_pathogenic_high_confidence_report.tsv" \
  --out_gene_summary_tsv "${REPORT_OUT}/male_infertility_clinvar_pathogenic_high_confidence_gene_summary.tsv"

test -s "${REPORT_OUT}/male_infertility_clinvar_pathogenic_high_confidence_gene_summary.tsv"

###############################################################################
# 7) Overlap: use your existing combined omics table
# Your file has gene_symbol_norm, not gene_symbol, so this requires the overlap
# script to support --gene_symbol_column (as we discussed earlier).
###############################################################################

TESTIS_GENE_TSV="${BASE_DIR}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv"
test -s "${TESTIS_GENE_TSV}"

python "${SCRIPTS_DIR}/overlap_testis_specific_with_clinvar_tiers.py" \
  --testis_gene_tsv "${TESTIS_GENE_TSV}" \
  --gene_symbol_column gene_symbol_norm \
  --clinvar_best_tsv "${CLINVAR_COLLAPSE_OUT}/male_infertility_clinvar_variants_best.tsv" \
  --clinvar_hc_pathogenic_tsv "${CLINVAR_HC_OUT}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_dir "${OVERLAP_OUT}"

echo "Pipeline completed successfully."
echo "Results written to: ${RESULTS_DIR}"



