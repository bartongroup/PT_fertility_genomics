


cd ~/data/2026_sperm_Gates/hpo_data/male_infertility_out 

python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/extract_male_infertility_genes_hpo_ontology.py   --hpo_dir ~/data/2026_sperm_Gates/hpo_data   --out_dir ~/data/2026_sperm_Gates/hpo_data/male_infertility_out   --include_diseases   --require_male


python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/extract_male_infertility_genes_hpo_ontology.py \
  --hpo_dir ~/data/2026_sperm_Gates/hpo_data \
  --out_dir ~/data/2026_sperm_Gates/hpo_data/male_infertility_out_phrase_seeds \
  --include_diseases \
  --require_male \
  --use_phrase_seeds

python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/extract_gene_list_from_hpo_outputs.py \
  --genes_summary_tsv ~/data/2026_sperm_Gates/hpo_data/male_infertility_out_phrase_seeds/male_infertility_genes_summary.tsv \
  --out_dir ~/data/2026_sperm_Gates/hpo_data/male_infertility_out_phrase_seeds/gene_lists


 python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/filter_clinvar_variant_summary_by_gene_list.py --variant_summary_gz ~/data/2026_sperm_Gates/clinvar/variant_summary.txt.gz \
  --gene_symbols_tsv ~/data/2026_sperm_Gates/hpo_data/male_infertility_out_phrase_seeds/male_infertility_gene_symbols.tsv \
  --out_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variant_summary.tsv \
  --out_summary_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variant_summary_counts.tsv

 python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/filter_clinvar_filtered_by_phenotype.py \
  --in_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variant_summary.tsv \
  --out_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_by_phenotype.tsv \
  --out_minimal_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_by_phenotype_minimal.tsv \
  --out_summary_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_by_phenotype_counts.tsv


python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/collapse_clinvar_variants_to_best.py \
  --in_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_by_phenotype.tsv \
  --out_best_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best.tsv \
  --out_best_pathogenic_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best_pathogenic.tsv \
  --out_summary_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best_counts.tsv \
  --group_key AlleleID

python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/filter_clinvar_best_high_confidence.py \
  --in_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best.tsv \
  --out_high_conf_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence.tsv \
  --out_high_conf_pathogenic_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv \
  --out_counts_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_counts.tsv

python ~/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/make_clinvar_high_confidence_report.py \
  --in_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv \
  --out_report_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_pathogenic_high_confidence_report.tsv \
  --out_gene_summary_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_pathogenic_high_confidence_gene_summary.tsv


python overlap_testis_specific_with_clinvar_tiers.py \
  --testis_gene_tsv /path/to/testis_specific_gene_symbols.tsv \
  --clinvar_best_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best.tsv \
  --clinvar_hc_pathogenic_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv \
  --out_dir ~/data/2026_sperm_Gates/overlap_outputs
