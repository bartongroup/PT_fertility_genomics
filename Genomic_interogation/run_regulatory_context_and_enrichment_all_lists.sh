#!/usr/bin/env bash
set -euo pipefail

script_dir="/home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/Genomic_interogation"
resources_dir="${HOME}/data/2026_sperm_Gates/genome_resources"
gene_lists_dir="${resources_dir}/gene_lists"

encode_dir="${resources_dir}/encode_tracks"
results_dir="${resources_dir}/enrichment_results_regulatory"
mkdir -p "${encode_dir}" "${results_dir}"

base_features="${resources_dir}/gene_context_features_universe.tsv"
if [[ ! -s "${base_features}" ]]; then
  printf '[%s] ERROR: Missing %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${base_features}" >&2
  exit 1
fi

# 1) Download a small set of tracks (edit biosamples/targets as you like)
python "${script_dir}/download_encode_bigwig_tracks.py" \
  --out_dir "${encode_dir}" \
  --assembly "GRCh38" \
  --assay_title "Repli-seq" \
  --biosample "H1-hESC" \
  --max_files_per_query 2 \
  --verbose

# Optional: Lamin B1 ChIP-seq signal as a LAD proxy (not true LAD calls, but informative)
python "${script_dir}/download_encode_bigwig_tracks.py" \
  --out_dir "${encode_dir}" \
  --assembly "GRCh38" \
  --assay_title "ChIP-seq" \
  --target "Lamin B1" \
  --biosample "H1-hESC" \
  --max_files_per_query 2 \
  --verbose

manifest="${encode_dir}/encode_tracks_manifest.tsv"
aug_features="${resources_dir}/gene_context_features_universe_plus_tracks.tsv"

# 2) Add bigWig-derived features
python "${script_dir}/add_bigwig_signal_features.py" \
  --features_tsv "${base_features}" \
  --bigwig_manifest_tsv "${manifest}" \
  --tss_window_bp 50000 \
  --include_max \
  --out_tsv "${aug_features}" \
  --verbose

# 3) Iterate enrichment across all gene lists
shopt -s nullglob
for gene_list in "${gene_lists_dir}"/*.tsv; do
  base="$(basename "${gene_list}")"
  stem="${base%.tsv}"
  out_dir="${results_dir}/${stem}"
  mkdir -p "${out_dir}"
  out_prefix="${out_dir}/${stem}_vs_null_regulatory"

  python "${script_dir}/matched_set_enrichment.py" \
    --features_tsv "${aug_features}" \
    --target_gene_list_tsv "${gene_list}" \
    --gene_key_column "gene_key" \
    --n_permutations 10000 \
    --n_length_bins 10 \
    --match_on_gc \
    --n_gc_bins 10 \
    --out_prefix "${out_prefix}"
done
