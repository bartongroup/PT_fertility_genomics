#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RUN_ALL.sh
#
# Orchestrates the sperm genomics prioritisation pipeline with NO downloads.
# Fails fast if expected inputs are missing.
###############################################################################

###############################################################################
# CONFIG
###############################################################################

BASE_DIR="${HOME}/data/2026_sperm_Gates"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Script folders (repo)
RNASEQ_DIR="${REPO_DIR}/RNAseq"
PHENO_DIR="${REPO_DIR}/find_male_phenotypic_disease"
GENOMIC_DIR="${REPO_DIR}/Genomic_interogation"
FINAL_DIR="${REPO_DIR}/final_summary"

# Data folders (base)
GTEX_DIR="${BASE_DIR}/GTEx"
CLINVAR_DIR="${BASE_DIR}/clinvar"
HPO_DIR="${BASE_DIR}/hpo_data"
SPERM_DIR="${BASE_DIR}/sperm_RNAseq"
PROT_DIR="${BASE_DIR}/proteomics"
RESOURCES_DIR="${BASE_DIR}/genome_resources"
UNIPROT_DIR="${RESOURCES_DIR}/uniprot"

# Key input files (edit these if your filenames differ)
GTEX_GCT="${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"
GTEX_SAMPLE_ATTR="${GTEX_DIR}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
HGNC_TSV="${BASE_DIR}/hgnc_complete_set.txt"

CLINVAR_GENE_COND="${CLINVAR_DIR}/gene_condition_source_id"
CLINVAR_VARIANT_SUMMARY_GZ="${CLINVAR_DIR}/variant_summary.txt.gz"

HPO_GENES_TO_PHENO="${HPO_DIR}/genes_to_phenotype.txt"

SPERM_TPM_GZ="${SPERM_DIR}/GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
PROT_XLSX="${PROT_DIR}/SPZ-HomoSapiens.xlsx"

# UniProt resources for biochemical accessibility
UNIPROT_ANNOTATION_TSV="${UNIPROT_DIR}/uniprotkb_reviewed_true.tsv"
GENE_TO_UNIPROT_TSV="${BASE_DIR}/results/FULL_SUMMARY/gene_to_uniprot.tsv"

# Thresholds
MIN_TPM_PRESENT="5"
SPERM_TPM_THRESHOLD="0.1"
MIN_TAU="0.95"
MIN_TESTIS_TPM="5"
MIN_SPERM_TPM="0.1"

# Run naming (used to create output folders)
TESTIS_RUN_ID="testis_tau${MIN_TAU}_testisTPM${MIN_TESTIS_TPM}_presentTPM${MIN_TPM_PRESENT}_spermTPM${MIN_SPERM_TPM}"

# Output locations
RESULTS_DIR="${BASE_DIR}/results"
TESTIS_RESULTS_DIR="${RESULTS_DIR}/${TESTIS_RUN_ID}"
FULL_SUMMARY_DIR="${RESULTS_DIR}/FULL_SUMMARY"

LOG_DIR="${TESTIS_RESULTS_DIR}/logs"
mkdir -p "${LOG_DIR}" "${FULL_SUMMARY_DIR}"

###############################################################################
# Helpers
###############################################################################

timestamp() { date -u +"%Y-%m-%dT%H:%M:%SZ"; }

log() {
  echo "[$(timestamp)] $*"
}

require_file() {
  local path="${1}"
  if [ ! -s "${path}" ]; then
    echo "ERROR: missing/empty required file: ${path}" 1>&2
    exit 1
  fi
}

require_dir() {
  local path="${1}"
  if [ ! -d "${path}" ]; then
    echo "ERROR: missing required directory: ${path}" 1>&2
    exit 1
  fi
}

###############################################################################
# 0) Sanity checks (no downloads, so validate everything up front)
###############################################################################

log "Repo dir: ${REPO_DIR}"
log "Base dir: ${BASE_DIR}"
log "Testis run id: ${TESTIS_RUN_ID}"

require_dir "${RNASEQ_DIR}"
require_dir "${PHENO_DIR}"
require_dir "${GENOMIC_DIR}"
require_dir "${FINAL_DIR}"

require_file "${GTEX_GCT}"
require_file "${GTEX_SAMPLE_ATTR}"
require_file "${HGNC_TSV}"
require_file "${CLINVAR_GENE_COND}"
require_file "${CLINVAR_VARIANT_SUMMARY_GZ}"
require_file "${HPO_GENES_TO_PHENO}"
require_file "${SPERM_TPM_GZ}"
require_file "${PROT_XLSX}"

###############################################################################
# 1) Male infertility gene list + ClinVar filtering (HPO -> ClinVar tiers)
###############################################################################

log "Step 1: HPO expansion + ClinVar filtering"

PHENO_OUT_01="${RESULTS_DIR}/01_hpo"
PHENO_OUT_02="${RESULTS_DIR}/02_clinvar_gene_filter"
PHENO_OUT_03="${RESULTS_DIR}/03_clinvar_phenotype_filter"
PHENO_OUT_04="${RESULTS_DIR}/04_clinvar_collapsed"
PHENO_OUT_05="${RESULTS_DIR}/05_clinvar_high_confidence"
PHENO_OUT_06="${RESULTS_DIR}/06_reports"

mkdir -p "${PHENO_OUT_01}" "${PHENO_OUT_02}" "${PHENO_OUT_03}" \
  "${PHENO_OUT_04}" "${PHENO_OUT_05}" "${PHENO_OUT_06}"

python "${PHENO_DIR}/extract_male_infertility_genes_hpo_ontology.py" \
  --hpo_dir "${HPO_DIR}" \
  --out_dir "${PHENO_OUT_01}" \
  --include_diseases \
  --require_male \
  --use_phrase_seeds

python "${PHENO_DIR}/extract_gene_list_from_hpo_outputs.py" \
  --genes_summary_tsv "${PHENO_OUT_01}/male_infertility_genes_summary.tsv" \
  --out_dir "${PHENO_OUT_01}/gene_lists"

HPO_GENE_SYMBOLS_TSV="${PHENO_OUT_01}/gene_lists/male_infertility_gene_symbols.tsv"
require_file "${HPO_GENE_SYMBOLS_TSV}"

python "${PHENO_DIR}/filter_clinvar_variant_summary_by_gene_list.py" \
  --variant_summary_gz "${CLINVAR_VARIANT_SUMMARY_GZ}" \
  --gene_symbols_tsv "${HPO_GENE_SYMBOLS_TSV}" \
  --out_tsv "${PHENO_OUT_02}/male_infertility_clinvar_variant_summary.tsv" \
  --out_summary_tsv "${PHENO_OUT_02}/male_infertility_clinvar_variant_summary_counts.tsv"

python "${PHENO_DIR}/filter_clinvar_filtered_by_phenotype.py" \
  --in_tsv "${PHENO_OUT_02}/male_infertility_clinvar_variant_summary.tsv" \
  --out_tsv "${PHENO_OUT_03}/male_infertility_clinvar_variants_by_phenotype.tsv" \
  --out_minimal_tsv "${PHENO_OUT_03}/male_infertility_clinvar_variants_by_phenotype_minimal.tsv" \
  --out_summary_tsv "${PHENO_OUT_03}/male_infertility_clinvar_variants_by_phenotype_counts.tsv"

python "${PHENO_DIR}/collapse_clinvar_variants_to_best.py" \
  --in_tsv "${PHENO_OUT_03}/male_infertility_clinvar_variants_by_phenotype.tsv" \
  --out_best_tsv "${PHENO_OUT_04}/male_infertility_clinvar_variants_best.tsv" \
  --out_best_pathogenic_tsv "${PHENO_OUT_04}/male_infertility_clinvar_variants_best_pathogenic.tsv" \
  --out_summary_tsv "${PHENO_OUT_04}/male_infertility_clinvar_variants_best_counts.tsv" \
  --group_key "AlleleID"

python "${PHENO_DIR}/filter_clinvar_best_high_confidence.py" \
  --in_tsv "${PHENO_OUT_04}/male_infertility_clinvar_variants_best.tsv" \
  --out_high_conf_tsv "${PHENO_OUT_05}/male_infertility_clinvar_variants_high_confidence.tsv" \
  --out_high_conf_pathogenic_tsv "${PHENO_OUT_05}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_counts_tsv "${PHENO_OUT_05}/male_infertility_clinvar_variants_high_confidence_counts.tsv"

python "${PHENO_DIR}/make_clinvar_high_confidence_report.py" \
  --in_tsv "${PHENO_OUT_05}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_report_tsv "${PHENO_OUT_06}/male_infertility_clinvar_pathogenic_high_confidence_report.tsv" \
  --out_gene_summary_tsv "${PHENO_OUT_06}/male_infertility_clinvar_pathogenic_high_confidence_gene_summary.tsv"

###############################################################################
# 2) GTEx testis specificity + sperm overlap + annotation table
###############################################################################

log "Step 2: GTEx testis specificity + sperm overlap + omics annotation"

OUT_01_GTEX="${TESTIS_RESULTS_DIR}/01_gtex_specificity"
OUT_02_HGNC="${TESTIS_RESULTS_DIR}/02_hgnc_mapping"
OUT_03_CLINVAR_ANN="${TESTIS_RESULTS_DIR}/03_clinvar_annotation"
OUT_04_SPERM="${TESTIS_RESULTS_DIR}/04_sperm_intersect"
OUT_05_HC="${TESTIS_RESULTS_DIR}/05_high_confidence"
OUT_06_FUNCTION="${TESTIS_RESULTS_DIR}/06_functional_annotation"
OUT_07_HPO="${TESTIS_RESULTS_DIR}/07_hpo_annotation"
OUT_08_PROT="${TESTIS_RESULTS_DIR}/08_proteomics_annotation"
OUT_09_FINAL_HC="${TESTIS_RESULTS_DIR}/09_high_confidence_final"
OUT_10_OVERLAP="${TESTIS_RESULTS_DIR}/10_overlap_clinvar_tiers"

mkdir -p "${OUT_01_GTEX}" "${OUT_02_HGNC}" "${OUT_03_CLINVAR_ANN}" \
  "${OUT_04_SPERM}" "${OUT_05_HC}" "${OUT_06_FUNCTION}" "${OUT_07_HPO}" \
  "${OUT_08_PROT}" "${OUT_09_FINAL_HC}" "${OUT_10_OVERLAP}"

python "${RNASEQ_DIR}/gtex_testis_specificity_logged.py" \
  --gct_path "${GTEX_GCT}" \
  --sample_attributes_path "${GTEX_SAMPLE_ATTR}" \
  --target_tissue "Testis" \
  --min_tpm_present "${MIN_TPM_PRESENT}" \
  --out_ranked_tsv "${OUT_01_GTEX}/gtex_v11_testis_specificity_ranked.tsv" \
  --out_tissue_medians_tsv "${OUT_01_GTEX}/gtex_v11_tissue_medians.tsv" \
  --log_level "INFO" \
  --log_path "${LOG_DIR}/gtex_testis_specificity_presentTPM${MIN_TPM_PRESENT}.log"

python "${RNASEQ_DIR}/add_hgnc_symbols.py" \
  --ranked_genes_tsv "${OUT_01_GTEX}/gtex_v11_testis_specificity_ranked.tsv" \
  --hgnc_tsv "${HGNC_TSV}" \
  --output_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --log_path "${LOG_DIR}/hgnc_mapping.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/annotate_testis_genes_with_clinvar.py" \
  --ranked_genes_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --gene_symbol_col "hgnc_symbol" \
  --clinvar_gene_condition_path "${CLINVAR_GENE_COND}" \
  --clinvar_variant_summary_gz "${CLINVAR_VARIANT_SUMMARY_GZ}" \
  --include_variant_counts \
  --output_tsv "${OUT_03_CLINVAR_ANN}/gtex_testis_genes_clinvar_annotated_hgnc.tsv" \
  --log_path "${LOG_DIR}/clinvar_annotated_hgnc.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/map_sperm_entrez_to_hgnc_and_intersect.py" \
  --sperm_tpm_path "${SPERM_TPM_GZ}" \
  --hgnc_tsv "${HGNC_TSV}" \
  --gtex_ranked_with_hgnc_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --output_sperm_mapped_tsv "${OUT_04_SPERM}/GSE40181_sperm_tpm_mapped_hgnc.tsv" \
  --output_gtex_with_sperm_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --sperm_tpm_threshold "${SPERM_TPM_THRESHOLD}" \
  --log_path "${LOG_DIR}/sperm_entrez_to_hgnc.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/extract_high_confidence_testis_sperm_genes.py" \
  --input_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "${OUT_05_HC}/high_confidence_testis_sperm_genes.tsv" \
  --min_tau "${MIN_TAU}" \
  --min_testis_tpm "${MIN_TESTIS_TPM}" \
  --min_sperm_tpm "${MIN_SPERM_TPM}" \
  --log_path "${LOG_DIR}/high_confidence_testis_sperm.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/add_functional_annotation_mygene.py" \
  --input_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "${OUT_06_FUNCTION}/gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --symbol_col "hgnc_symbol" \
  --cache_jsonl "${TESTIS_RESULTS_DIR}/mygene_cache.jsonl" \
  --batch_size "100" \
  --sleep_seconds "0.3" \
  --log_path "${LOG_DIR}/mygene_annotation.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/add_hpo_annotation.py" \
  --input_tsv "${OUT_06_FUNCTION}/gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --hpo_genes_to_phenotype_tsv "${HPO_GENES_TO_PHENO}" \
  --output_tsv "${OUT_07_HPO}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --log_path "${LOG_DIR}/hpo_annotation.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/add_proteomics_annotation.py" \
  --input_gene_tsv "${OUT_07_HPO}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --proteomics_xlsx "${PROT_XLSX}" \
  --output_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --require_fdr_confidence "High" \
  --restrict_species_name "Homo sapiens" \
  --log_path "${LOG_DIR}/proteomics_annotation.log" \
  --log_level "INFO"

python "${RNASEQ_DIR}/extract_high_confidence_testis_sperm_genes.py" \
  --input_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --output_tsv "${OUT_09_FINAL_HC}/high_confidence_testis_sperm_genes_function_hpo_prot.tsv" \
  --min_tau "${MIN_TAU}" \
  --min_testis_tpm "${MIN_TESTIS_TPM}" \
  --min_sperm_tpm "${MIN_SPERM_TPM}" \
  --log_path "${LOG_DIR}/high_confidence_testis_sperm_function_hpo_prot.log" \
  --log_level "INFO"

python "${PHENO_DIR}/overlap_testis_specific_with_clinvar_tiers.py" \
  --testis_gene_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --gene_symbol_column "gene_symbol_norm" \
  --clinvar_best_tsv "${PHENO_OUT_04}/male_infertility_clinvar_variants_best.tsv" \
  --clinvar_hc_pathogenic_tsv "${PHENO_OUT_05}/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv" \
  --out_dir "${OUT_10_OVERLAP}"

###############################################################################
# 3) Optional: OpenTargets tractability annotation (only if universe table exists)
###############################################################################

log "Step 3: OpenTargets tractability annotation (optional)"

GENE_CONTEXT_UNIVERSE="${RESOURCES_DIR}/gene_context_features_universe.tsv"
GENE_CONTEXT_UNIVERSE_TRACT="${RESOURCES_DIR}/gene_context_features_universe_plus_tractability.tsv"

if [ -s "${GENE_CONTEXT_UNIVERSE}" ]; then
  python "${GENOMIC_DIR}/annotate_opentargets_tractability.py" \
    --in_tsv "${GENE_CONTEXT_UNIVERSE}" \
    --ensembl_column "gene_id" \
    --out_tsv "${GENE_CONTEXT_UNIVERSE_TRACT}" \
    --verbose
else
  log "Skipping tractability: missing ${GENE_CONTEXT_UNIVERSE}"
fi

###############################################################################
# 4) Build master Excel summary + plots
###############################################################################

log "Step 4: Build master summary workbook + plots"

MASTER_XLSX="${FULL_SUMMARY_DIR}/SUMMARY_fertility_evidence.xlsx"

python "${FINAL_DIR}/make_master_fertility_results_summary.py" \
  --base_dir "${BASE_DIR}" \
  --testis_run_id "${TESTIS_RUN_ID}" \
  --out_xlsx "${MASTER_XLSX}"

PLOTS_DIR="${FULL_SUMMARY_DIR}/plots"
mkdir -p "${PLOTS_DIR}"

python "${FINAL_DIR}/make_fertility_gene_set_plots.py" \
  --in_xlsx "${MASTER_XLSX}" \
  --sheet_name "Tier_Summary_With_Omics" \
  --out_dir "${PLOTS_DIR}"

###############################################################################
# 5) Add biochemical accessibility (only if UniProt inputs exist)
###############################################################################

log "Step 5: Biochemical accessibility annotation (optional)"

BIOCHEM_XLSX="${FULL_SUMMARY_DIR}/SUMMARY_fertility_evidence.biochem.xlsx"

if [ -s "${UNIPROT_ANNOTATION_TSV}" ] && [ -s "${GENE_TO_UNIPROT_TSV}" ]; then
  python "${FINAL_DIR}/annotate_sperm_biochemical_accessibility.py" \
    --excel_in "${MASTER_XLSX}" \
    --excel_out "${BIOCHEM_XLSX}" \
    --uniprot_annotation_tsv "${UNIPROT_ANNOTATION_TSV}" \
    --gene_to_uniprot_tsv "${GENE_TO_UNIPROT_TSV}" \
    --verbose
else
  log "Skipping biochemical accessibility: missing UniProt inputs"
  log "  expected: ${UNIPROT_ANNOTATION_TSV}"
  log "  expected: ${GENE_TO_UNIPROT_TSV}"
  BIOCHEM_XLSX="${MASTER_XLSX}"
fi

###############################################################################
# 6) Final target prioritisation (only if tractability table exists)
###############################################################################

log "Step 6: Final target prioritisation (optional)"

if [ -s "${GENE_CONTEXT_UNIVERSE_TRACT}" ]; then
  python "${FINAL_DIR}/prioritise_druggable_sperm_targets.py" \
    --excel_in "${BIOCHEM_XLSX}" \
    --tractability_tsv "${GENE_CONTEXT_UNIVERSE_TRACT}" \
    --tractability_gene_id_column "gene_id" \
    --strip_ensembl_version \
    --min_memberships "2" \
    --require_testis_dominance \
    --top_n "200" \
    --out_prefix "${FULL_SUMMARY_DIR}/sperm_target_priorities_from_master" \
    --verbose
else
  log "Skipping prioritisation: missing ${GENE_CONTEXT_UNIVERSE_TRACT}"
fi

log "All done."
log "Testis results: ${TESTIS_RESULTS_DIR}"
log "Master workbook: ${MASTER_XLSX}"
log "Final workbook used for prioritisation: ${BIOCHEM_XLSX}"