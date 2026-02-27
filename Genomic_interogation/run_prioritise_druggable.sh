

python prioritise_druggable_sperm_targets.py \
  --gene_lists_dir "${HOME}/data/2026_sperm_Gates/genome_resources/gene_lists" \
  --features_tsv "${HOME}/data/2026_sperm_Gates/genome_resources/gene_context_features_universe_plus_tracks.tsv" \
  --proteomics_tsv "${HOME}/data/2026_sperm_Gates/genome_resources/gene_lists/proteomics_internal__public__any__A_all.tsv" \
  --tractability_tsv "${HOME}/data/2026_sperm_Gates/genome_resources/gene_context_features_universe_plus_tractability.tsv" \
  --novelty_exclude_list_regex "literature" \
  --min_lists 1 \
  --out_tsv "${HOME}/data/2026_sperm_Gates/genome_resources/sperm_target_priorities.tsv" \
  --verbose
