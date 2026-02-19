#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# download_data_linear.sh
#
# Downloads:
#   - GTEx v11 gene TPM GCT (RNASeQC)
#   - GSE40181 sperm RNA-seq normalised counts (TPM)
#   - ClinVar tab-delimited files
#
# Behaviour:
#   - Creates directories if missing
#   - Skips download if file exists and is non-empty
#   - Uses robust wget options
###############################################################################

# Check wget exists
if ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: wget not found on PATH"
    exit 1
fi

echo "Starting downloads..."

echo "some of these fail and need to be downloaded manually -  sorry. See the script for details."

###############################################################################
# GTEx
###############################################################################

mkdir -p GTEx

GTEX_FILE="GTEx/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
GTEX_URL="https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"

if [[ -s "${GTEX_FILE}" ]]; then
    echo "GTEx file already present. Skipping."
else
    echo "Downloading GTEx data..."
    wget \
        --show-progress \
        --progress=dot:giga \
        --tries=5 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        --continue \
        --output-document="${GTEX_FILE}" \
        "${GTEX_URL}"
fi




###############################################################################
# more proteomics
###############################################################################

mkdir -p PXD014618_proteomics

cd PXD014618_proteomics

wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014618/3506_AS_7738_RedAlkC18_Acclaim50_14_ms2.OUTPUT_TABLE
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD014618/3506_AS_7738_RedAlkC18_Acclaim50_07_ms2.OUTPUT_TABLE

cd ../

###############################################################################
# Sperm RNA-seq (GSE40181)
###############################################################################

mkdir -p sperm_RNAseq

SPERM_FILE="sperm_RNAseq/GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
SPERM_URL="https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE40181&format=file&file=GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"

if [[ -s "${SPERM_FILE}" ]]; then
    echo "GSE40181 file already present. Skipping."
else
    echo "Downloading GSE40181 sperm RNA-seq data..."
    wget \
        --show-progress \
        --progress=dot:giga \
        --tries=5 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        --continue \
        --output-document="${SPERM_FILE}" \
        "${SPERM_URL}"
fi

###############################################################################
# HPO (manual step)
###############################################################################

mkdir -p hpo_data

echo ""
echo "HPO files must be downloaded manually into hpo_data/:"
echo "  - genes_to_phenotype.txt"
echo "  - hp.json"
echo "  - hp.obo"
echo "  - phenotype.hpoa"
echo "Sources:"
echo "  https://hpo.jax.org/data/ontology"
echo "  https://hpo.jax.org/data/annotations"
echo ""

###############################################################################
# ClinVar
###############################################################################

mkdir -p clinvar

CLINVAR_SUMMARY_FILE="clinvar/variant_summary.txt.gz"
CLINVAR_SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

if [[ -s "${CLINVAR_SUMMARY_FILE}" ]]; then
    echo "ClinVar variant_summary already present. Skipping."
else
    echo "Downloading ClinVar variant_summary..."
    wget \
        --show-progress \
        --progress=dot:giga \
        --tries=5 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        --continue \
        --output-document="${CLINVAR_SUMMARY_FILE}" \
        "${CLINVAR_SUMMARY_URL}"
fi

CLINVAR_ALLELE_FILE="clinvar/variation_allele.txt.gz"
CLINVAR_ALLELE_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variation_allele.txt.gz"

if [[ -s "${CLINVAR_ALLELE_FILE}" ]]; then
    echo "ClinVar variation_allele already present. Skipping."
else
    echo "Downloading ClinVar variation_allele..."
    wget \
        --show-progress \
        --progress=dot:giga \
        --tries=5 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        --continue \
        --output-document="${CLINVAR_ALLELE_FILE}" \
        "${CLINVAR_ALLELE_URL}"
fi

echo ""
echo "All downloads complete."
