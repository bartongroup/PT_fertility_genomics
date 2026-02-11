#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# run_testis_specificity_and_overlap.sh
#
# Linear pipeline:
#  1) GTEx testis specificity (tau + TPM present)
#  2) Add HGNC symbols
#  3) ClinVar annotation (gene_condition + variant_summary)
#  4) Map sperm RNA-seq (GSE40181) to HGNC and intersect with GTEx
#  5) High-confidence gene extraction (tau/testis TPM/sperm TPM thresholds)
#  6) Add functional annotation (mygene)
#  7) Add HPO annotation
#  8) Add proteomics annotation
#  9) High-confidence extraction on the fully annotated table
# 10) Overlap with ClinVar male-infertility tiers
#
# Outputs:
#   All outputs go to: ~/data/2026_sperm_Gates/results/<run_id>/
#
# No functions. No positional arguments.
###############################################################################

###############################################################################
# Base paths (current on-disk layout)
###############################################################################

BASE_DIR="${HOME}/data/2026_sperm_Gates"

SCRIPTS_DIR="${BASE_DIR}/scripts"
PHENO_SCRIPTS_DIR="${SCRIPTS_DIR}/find_male_phenotypic_disease"

GTEX_DIR="${BASE_DIR}/GTEx"
CLINVAR_DIR="${BASE_DIR}/clinvar"
HPO_DIR="${BASE_DIR}/hpo_data"
SPERM_DIR="${BASE_DIR}/sperm_RNAseq"
PROT_DIR="${BASE_DIR}/proteomics"

###############################################################################
# Key input files (from download script + your existing files)
###############################################################################

GTEX_GCT="${GTEX_DIR}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"
GTEX_SAMPLE_ATTR="${GTEX_DIR}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"

HGNC_TSV="${BASE_DIR}/hgnc_complete_set.txt"

CLINVAR_GENE_COND="${CLINVAR_DIR}/gene_condition_source_id"
CLINVAR_VARIANT_SUMMARY_GZ="${CLINVAR_DIR}/variant_summary.txt.gz"

SPERM_TPM_GZ="${SPERM_DIR}/GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
PROT_XLSX="${PROT_DIR}/SPZ-HomoSapiens.xlsx"

# From your male infertility ClinVar pipeline (already generated)
CLINVAR_MALE_BEST="${CLINVAR_DIR}/clinvar_filtered/male_infertility_clinvar_variants_best.tsv"
CLINVAR_MALE_HC_PATH="${CLINVAR_DIR}/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv"

###############################################################################
# Thresholds (edit here)
###############################################################################

MIN_TPM_PRESENT="5"
SPERM_TPM_THRESHOLD="0.1"

MIN_TAU="0.95"
MIN_TESTIS_TPM="5"
MIN_SPERM_TPM="0.1"

###############################################################################
# Results folder naming (include thresholds)
###############################################################################

RUN_ID="testis_tau${MIN_TAU}_testisTPM${MIN_TESTIS_TPM}_presentTPM${MIN_TPM_PRESENT}_spermTPM${MIN_SPERM_TPM}"
RESULTS_DIR="${BASE_DIR}/results/${RUN_ID}"

OUT_01_GTEX="${RESULTS_DIR}/01_gtex_specificity"
OUT_02_HGNC="${RESULTS_DIR}/02_hgnc_mapping"
OUT_03_CLINVAR_ANN="${RESULTS_DIR}/03_clinvar_annotation"
OUT_04_SPERM="${RESULTS_DIR}/04_sperm_intersect"
OUT_05_HC="${RESULTS_DIR}/05_high_confidence"
OUT_06_FUNCTION="${RESULTS_DIR}/06_functional_annotation"
OUT_07_HPO="${RESULTS_DIR}/07_hpo_annotation"
OUT_08_PROT="${RESULTS_DIR}/08_proteomics_annotation"
OUT_09_FINAL_HC="${RESULTS_DIR}/09_high_confidence_final"
OUT_10_OVERLAP="${RESULTS_DIR}/10_overlap_clinvar_tiers"
LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p "${OUT_01_GTEX}" "${OUT_02_HGNC}" "${OUT_03_CLINVAR_ANN}" \
  "${OUT_04_SPERM}" "${OUT_05_HC}" "${OUT_06_FUNCTION}" "${OUT_07_HPO}" \
  "${OUT_08_PROT}" "${OUT_09_FINAL_HC}" "${OUT_10_OVERLAP}" "${LOG_DIR}"

###############################################################################
# Input checks (fail fast)
###############################################################################

test -s "${GTEX_GCT}"
test -s "${GTEX_SAMPLE_ATTR}"
test -s "${HGNC_TSV}"
test -s "${CLINVAR_VARIANT_SUMMARY_GZ}"
test -d "${CLINVAR_GENE_COND}"
test -s "${SPERM_TPM_GZ}"
test -s "${PROT_XLSX}"
test -s "${CLINVAR_MALE_BEST}"
test -s "${CLINVAR_MALE_HC_PATH}"
test -s "${HPO_DIR}/genes_to_phenotype.txt"

###############################################################################
# 0) Activate environment (optional)
###############################################################################
# If conda activation is required in your cluster job shell, keep this.
# If you run inside an already-activated env, you can comment it out.

#source "${HOME}/miniconda3/etc/profile.d/conda.sh"
#conda activate python3.9

###############################################################################
# 1) GTEx testis specificity ranking
###############################################################################

python "${SCRIPTS_DIR}/gtex_testis_specificity_logged.py" \
  --gct_path "${GTEX_GCT}" \
  --sample_attributes_path "${GTEX_SAMPLE_ATTR}" \
  --target_tissue "Testis" \
  --min_tpm_present "${MIN_TPM_PRESENT}" \
  --out_ranked_tsv "${OUT_01_GTEX}/gtex_v11_testis_specificity_ranked.tsv" \
  --out_tissue_medians_tsv "${OUT_01_GTEX}/gtex_v11_tissue_medians.tsv" \
  --log_level "INFO" \
  --log_path "${LOG_DIR}/gtex_testis_specificity_presentTPM${MIN_TPM_PRESENT}.log"

###############################################################################
# 2) Add HGNC mapping
###############################################################################

python "${SCRIPTS_DIR}/add_hgnc_symbols.py" \
  --ranked_genes_tsv "${OUT_01_GTEX}/gtex_v11_testis_specificity_ranked.tsv" \
  --hgnc_tsv "${HGNC_TSV}" \
  --output_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --log_path "${LOG_DIR}/hgnc_mapping.log" \
  --log_level "INFO"

###############################################################################
# 3) Annotate with ClinVar
###############################################################################

python "${SCRIPTS_DIR}/annotate_testis_genes_with_clinvar.py" \
  --ranked_genes_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --gene_symbol_col "hgnc_symbol" \
  --clinvar_gene_condition_path "${CLINVAR_GENE_COND}" \
  --clinvar_variant_summary_gz "${CLINVAR_VARIANT_SUMMARY_GZ}" \
  --include_variant_counts \
  --output_tsv "${OUT_03_CLINVAR_ANN}/gtex_testis_genes_clinvar_annotated_hgnc.tsv" \
  --log_path "${LOG_DIR}/clinvar_annotated_hgnc.log" \
  --log_level "INFO"

###############################################################################
# 4) Map sperm RNA-seq to HGNC and intersect
###############################################################################

python "${SCRIPTS_DIR}/map_sperm_entrez_to_hgnc_and_intersect.py" \
  --sperm_tpm_path "${SPERM_TPM_GZ}" \
  --hgnc_tsv "${HGNC_TSV}" \
  --gtex_ranked_with_hgnc_tsv "${OUT_02_HGNC}/gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --output_sperm_mapped_tsv "${OUT_04_SPERM}/GSE40181_sperm_tpm_mapped_hgnc.tsv" \
  --output_gtex_with_sperm_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --sperm_tpm_threshold "${SPERM_TPM_THRESHOLD}" \
  --log_path "${LOG_DIR}/sperm_entrez_to_hgnc.log" \
  --log_level "INFO"

###############################################################################
# 5) High-confidence extraction (GTEx + sperm)
###############################################################################

python "${SCRIPTS_DIR}/extract_high_confidence_testis_sperm_genes.py" \
  --input_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "${OUT_05_HC}/high_confidence_testis_sperm_genes.tsv" \
  --min_tau "${MIN_TAU}" \
  --min_testis_tpm "${MIN_TESTIS_TPM}" \
  --min_sperm_tpm "${MIN_SPERM_TPM}" \
  --log_path "${LOG_DIR}/high_confidence_testis_sperm.log" \
  --log_level "INFO"

###############################################################################
# 6) Add functional annotation (mygene)
###############################################################################

python "${SCRIPTS_DIR}/add_functional_annotation_mygene.py" \
  --input_tsv "${OUT_04_SPERM}/gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "${OUT_06_FUNCTION}/gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --symbol_col "hgnc_symbol" \
  --cache_jsonl "${RESULTS_DIR}/mygene_cache.jsonl" \
  --batch_size "100" \
  --sleep_seconds "0.3" \
  --log_path "${LOG_DIR}/mygene_annotation.log" \
  --log_level "INFO"

python "${SCRIPTS_DIR}/extract_high_confidence_testis_sperm_genes.py" \
  --input_tsv "${OUT_06_FUNCTION}/gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --output_tsv "${OUT_06_FUNCTION}/high_confidence_testis_sperm_genes_function.tsv" \
  --min_tau "${MIN_TAU}" \
  --min_testis_tpm "${MIN_TESTIS_TPM}" \
  --min_sperm_tpm "${MIN_SPERM_TPM}" \
  --log_path "${LOG_DIR}/high_confidence_testis_sperm_function.log" \
  --log_level "INFO"

###############################################################################
# 7) Add HPO annotation
###############################################################################

python "${SCRIPTS_DIR}/add_hpo_annotation.py" \
  --input_tsv "${OUT_06_FUNCTION}/gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --hpo_genes_to_phenotype_tsv "${HPO_DIR}/genes_to_phenotype.txt" \
  --output_tsv "${OUT_07_HPO}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --log_path "${LOG_DIR}/hpo_annotation.log" \
  --log_level "INFO"

###############################################################################
# 8) Add proteomics annotation
###############################################################################

python "${SCRIPTS_DIR}/add_proteomics_annotation.py" \
  --input_gene_tsv "${OUT_07_HPO}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --proteomics_xlsx "${PROT_XLSX}" \
  --output_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --require_fdr_confidence "High" \
  --restrict_species_name "Homo sapiens" \
  --log_path "${LOG_DIR}/proteomics_annotation.log" \
  --log_level "INFO"

###############################################################################
# 9) Final high-confidence extraction (after HPO + proteomics)
###############################################################################

python "${SCRIPTS_DIR}/extract_high_confidence_testis_sperm_genes.py" \
  --input_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --output_tsv "${OUT_09_FINAL_HC}/high_confidence_testis_sperm_genes_function_hpo_prot.tsv" \
  --min_tau "${MIN_TAU}" \
  --min_testis_tpm "${MIN_TESTIS_TPM}" \
  --min_sperm_tpm "${MIN_SPERM_TPM}" \
  --log_path "${LOG_DIR}/high_confidence_testis_sperm_function_hpo_prot.log" \
  --log_level "INFO"

###############################################################################
# 10) Overlap with ClinVar male infertility tiers
###############################################################################

python "${PHENO_SCRIPTS_DIR}/overlap_testis_specific_with_clinvar_tiers.py" \
  --testis_gene_tsv "${OUT_08_PROT}/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --gene_symbol_column "gene_symbol_norm" \
  --clinvar_best_tsv "${CLINVAR_MALE_BEST}" \
  --clinvar_hc_pathogenic_tsv "${CLINVAR_MALE_HC_PATH}" \
  --out_dir "${OUT_10_OVERLAP}"

echo "Completed successfully."
echo "Results directory: ${RESULTS_DIR}"
