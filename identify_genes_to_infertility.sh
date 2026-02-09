


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
