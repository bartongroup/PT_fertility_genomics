#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# download_data.sh
#
# Downloads:
#   - GTEx v11 gene TPM GCT (RNASeQC)
#   - GSE40181 sperm RNA-seq normalised counts (TPM)
#   - ClinVar tab-delimited files
#
# Behaviour:
#   - Creates target directories if missing
#   - Skips downloads if the target file already exists and is non-empty
#   - Uses wget with sensible retry options
#
# Notes:
#   - HPO files still need manual download (see section below).
###############################################################################

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${1}"
}

die() {
    printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${1}" >&2
    exit 1
}

need_cmd() {
    command -v "${1}" >/dev/null 2>&1 || die "Required command not found on PATH: ${1}"
}

download_if_missing() {
    """
    Download a URL to a destination path if the file is missing or empty.

    Parameters
    ----------
    url : str
        Source URL to download.
    dest : str
        Destination file path.
    """
    local url="${1}"
    local dest="${2}"

    if [[ -s "${dest}" ]]; then
        log "Present, skipping: ${dest}"
        return 0
    fi

    log "Downloading: ${url}"
    log "To: ${dest}"

    wget \
        --quiet \
        --show-progress \
        --progress=dot:giga \
        --tries=5 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        --continue \
        --output-document="${dest}" \
        "${url}"
}

main() {
    need_cmd wget

    # -------------------------------------------------------------------------
    # GTEx
    # -------------------------------------------------------------------------
    mkdir -p "GTEx"
    download_if_missing \
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz" \
        "GTEx/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"

    # -------------------------------------------------------------------------
    # Sperm RNA-seq (GSE40181)
    # -------------------------------------------------------------------------
    mkdir -p "sperm_RNAseq"
    download_if_missing \
        "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE40181&format=file&file=GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz" \
        "sperm_RNAseq/GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"

    # -------------------------------------------------------------------------
    # HPO data (manual)
    # -------------------------------------------------------------------------
    mkdir -p "hpo_data"
    log "HPO files are not downloaded automatically (manual step)."
    log "Download into: hpo_data/"
    log "Needed files:"
    log "  - genes_to_phenotype.txt"
    log "  - hp.json"
    log "  - hp.obo"
    log "  - phenotype.hpoa"
    log "Sources:"
    log "  - https://hpo.jax.org/data/ontology"
    log "  - https://hpo.jax.org/data/annotations"

    # -------------------------------------------------------------------------
    # ClinVar
    # -------------------------------------------------------------------------
    mkdir -p "clinvar"
    download_if_missing \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" \
        "clinvar/variant_summary.txt.gz"

    download_if_missing \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variation_allele.txt.gz" \
        "clinvar/variation_allele.txt.gz"

    log "All done."
}

main "$@"
