#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# download_genome_resources.sh
#
# Robust downloader for GRCh38 resources (UCSC + GENCODE).
#
# Why this is more robust on flaky/proxied HPC networks:
# - Uses curl with resume (-C -) so reruns continue partial downloads
# - Uses a "stall detector" so transfers don't hang forever
# - Retries essentially forever with backoff
# - Downloads into /tmp first (fast local scratch) then copies to target
#   (prevents GPFS metadata/write stalls from killing long downloads)
#
# Requirements:
# - curl
# - gzip (for integrity checks of .gz)
#
# Notes:
# - GENCODE GTF is not hard-coded; set GENCODE_GTF_URL below.
# - If your proxy dislikes HTTPS, switch UCSC URLs to http://
###############################################################################

TARGET_DIR="${1:-genome_resources}"
TMP_DIR="${TMPDIR:-/tmp}/${USER}_genome_resources_download"

# Stall detector settings:
# If speed < 1 KB/s for 60 seconds, abort and retry.
SPEED_LIMIT_BYTES=1024
SPEED_TIME_SECONDS=60

# Retry settings
RETRY_MAX=9999
RETRY_DELAY=5

mkdir -p "${TARGET_DIR}"
mkdir -p "${TMP_DIR}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${1}"
}

die() {
  printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${1}" >&2
  exit 1
}

need_cmd() {
  command -v "${1}" >/dev/null 2>&1 || die "Missing required command: ${1}"
}

download() {
  # Named args only (bash doesn't enforce, but keep it explicit)
  local url="${1}"
  local out_path="${2}"

  local out_dir
  out_dir="$(dirname "${out_path}")"
  mkdir -p "${out_dir}"

  if [[ -s "${out_path}" ]]; then
    log "Already present (non-empty), skipping: ${out_path}"
    return 0
  fi

  log "Downloading: ${url}"
  log "To: ${out_path}"

  # -L follow redirects
  # -C - resume
  # --speed-* abort if connection stalls (prevents infinite hang)
  # --retry... keep retrying on transient issues
  curl -L \
    --retry "${RETRY_MAX}" \
    --retry-delay "${RETRY_DELAY}" \
    --retry-connrefused \
    --speed-limit "${SPEED_LIMIT_BYTES}" \
    --speed-time "${SPEED_TIME_SECONDS}" \
    -C - \
    -o "${out_path}" \
    "${url}"

  if [[ ! -s "${out_path}" ]]; then
    die "Download failed or produced empty file: ${out_path}"
  fi
}

check_gzip() {
  local path="${1}"
  if [[ "${path}" == *.gz ]]; then
    log "Checking gzip integrity: ${path}"
    gzip -t "${path}" || die "gzip integrity check failed: ${path}"
  fi
}

copy_to_target() {
  local tmp_path="${1}"
  local final_path="${2}"

  if [[ -s "${final_path}" ]]; then
    log "Target already present (non-empty), skipping copy: ${final_path}"
    return 0
  fi

  log "Copying to target: ${final_path}"
  mkdir -p "$(dirname "${final_path}")"
  cp -v "${tmp_path}" "${final_path}"
}

###############################################################################
# 0) Basic checks
###############################################################################
need_cmd "curl"
need_cmd "gzip"

log "Target directory: ${TARGET_DIR}"
log "Temporary directory: ${TMP_DIR}"

###############################################################################
# 1) GENCODE GTF (set this explicitly)
#
# Go to: https://www.gencodegenes.org/human/
# Choose a release and copy the URL for the "comprehensive gene annotation (CHR)" GTF.
#
# Example (placeholder): export GENCODE_GTF_URL="https://.../gencode.vXX.annotation.chr.gtf.gz"
###############################################################################
GENCODE_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"

GENCODE_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"


if [[ -n "${GENCODE_GTF_URL}" ]]; then
  gtf_name="$(basename "${GENCODE_GTF_URL}")"
  tmp_gtf="${TMP_DIR}/${gtf_name}"
  final_gtf="${TARGET_DIR}/${gtf_name}"

  download "${GENCODE_GTF_URL}" "${tmp_gtf}"
  check_gzip "${tmp_gtf}"
  copy_to_target "${tmp_gtf}" "${final_gtf}"
else
  log "GENCODE_GTF_URL not set. Skipping GENCODE GTF download."
  log "Set it like:"
  log "  export GENCODE_GTF_URL=\"<paste the comprehensive.chr.gtf.gz URL>\""
fi

###############################################################################
# 2) hg38 FASTA (UCSC)
###############################################################################
UCSC_BASE="https://hgdownload.soe.ucsc.edu"

tmp_fa="${TMP_DIR}/hg38.fa.gz"
final_fa="${TARGET_DIR}/hg38.fa.gz"
download "${UCSC_BASE}/goldenpath/hg38/bigZips/latest/hg38.fa.gz" "${tmp_fa}"
check_gzip "${tmp_fa}"
copy_to_target "${tmp_fa}" "${final_fa}"

# Index (optional but useful). If UCSC index fails, pyfaidx can create one later.
tmp_fai="${TMP_DIR}/hg38.fa.gz.fai"
final_fai="${TARGET_DIR}/hg38.fa.gz.fai"
download "${UCSC_BASE}/goldenpath/hg38/bigZips/latest/hg38.fa.gz.fai" "${tmp_fai}"
copy_to_target "${tmp_fai}" "${final_fai}"

###############################################################################
# 3) RepeatMasker table dump
###############################################################################
tmp_rmsk="${TMP_DIR}/rmsk.txt.gz"
final_rmsk="${TARGET_DIR}/rmsk.txt.gz"
download "${UCSC_BASE}/goldenPath/hg38/database/rmsk.txt.gz" "${tmp_rmsk}"
check_gzip "${tmp_rmsk}"
copy_to_target "${tmp_rmsk}" "${final_rmsk}"

###############################################################################
# 4) Segmental duplications table dump
###############################################################################
tmp_segdup="${TMP_DIR}/genomicSuperDups.txt.gz"
final_segdup="${TARGET_DIR}/genomicSuperDups.txt.gz"
download "${UCSC_BASE}/goldenPath/hg38/database/genomicSuperDups.txt.gz" "${tmp_segdup}"
check_gzip "${tmp_segdup}"
copy_to_target "${tmp_segdup}" "${final_segdup}"

###############################################################################
# 5) Recombination map (deCODE average bigWig)
###############################################################################
tmp_recomb="${TMP_DIR}/recombAvg.bw"
final_recomb="${TARGET_DIR}/recombAvg.bw"
download "${UCSC_BASE}/gbdb/hg38/recombRate/recombAvg.bw" "${tmp_recomb}"
copy_to_target "${tmp_recomb}" "${final_recomb}"

log "All done."
log "Files in: ${TARGET_DIR}"
ls -lh "${TARGET_DIR}" || true


