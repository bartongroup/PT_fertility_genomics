#!/bin/bash
set -euo pipefail

cd ~/data/2026_sperm_Gates

conda activate python3.9

#######################################################
# Find testis-specific genes with sperm expression
#######################################################
# data downloaded from: (RNAseq)
# https://gtexportal.org/home/downloads/adult-gtex/metadata
# https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
# clinvar : https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
# Tau tissue-specificity score (0–1, higher = more tissue-specific)
# sperm RNAseq/:  GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz      
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40181



python  scripts/gtex_testis_specificity_logged.py \
  --gct_path "GTEx/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct" \
  --sample_attributes_path "GTEx/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt" \
  --target_tissue "Testis" \
  --min_tpm_present 5 \
  --out_ranked_tsv "gtex_v11_testis_specificity_ranked_tpm5.tsv" \
  --out_tissue_medians_tsv "gtex_v11_tissue_medians.tsv" \
  --log_level "INFO" \
  --log_path "logs/gtex_testis_specificity_tpm5.log"



python  scripts/add_hgnc_symbols.py \
  --ranked_genes_tsv "gtex_v11_testis_specificity_ranked.tsv" \
  --hgnc_tsv "hgnc_complete_set.txt" \
  --output_tsv "gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --log_path "logs/hgnc_mapping.log" \
  --log_level "INFO"



python  scripts/annotate_testis_genes_with_clinvar.py \
  --ranked_genes_tsv "gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --gene_symbol_col "hgnc_symbol" \
  --clinvar_gene_condition_path "clinvar/gene_condition_source_id" \
  --clinvar_variant_summary_gz "clinvar/variant_summary.txt.gz" \
  --include_variant_counts \
  --output_tsv "gtex_testis_genes_clinvar_annotated_hgnc.tsv" \
  --log_path "logs/clinvar_annotated_hgnc.log" \
  --log_level "INFO"


python  scripts/map_sperm_entrez_to_hgnc_and_intersect.py \
  --sperm_tpm_path "sperm_RNAseq/GSE40181_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz" \
  --hgnc_tsv "hgnc_complete_set.txt" \
  --gtex_ranked_with_hgnc_tsv "gtex_v11_testis_specificity_ranked_with_hgnc.tsv" \
  --output_sperm_mapped_tsv "GSE40181_sperm_tpm_mapped_hgnc.tsv" \
  --output_gtex_with_sperm_tsv "gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --sperm_tpm_threshold 0.1 \
  --log_path "logs/sperm_entrez_to_hgnc.log" \
  --log_level "INFO"


python  scripts/extract_high_confidence_testis_sperm_genes.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "high_confidence_testis_sperm_genes.tsv" \
  --min_tau 0.95 \
  --min_testis_tpm 5 \
  --min_sperm_tpm 0.1 \
  --log_path "logs/high_confidence_testis_sperm.log" \
  --log_level "INFO"



python  scripts/add_functional_annotation_mygene.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence.tsv" \
  --output_tsv "gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --symbol_col "hgnc_symbol" \
  --cache_jsonl "mygene_cache.jsonl" \
  --batch_size 100 \
  --sleep_seconds 0.3 \
  --log_path "logs/mygene_annotation.log" \
  --log_level "INFO"


python  scripts/extract_high_confidence_testis_sperm_genes.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --output_tsv "high_confidence_testis_sperm_genes_function.tsv" \
  --min_tau 0.95 \
  --min_testis_tpm 5 \
  --min_sperm_tpm 0.1 \
  --log_path "logs/high_confidence_testis_sperm_function.log" \
  --log_level "INFO"


python  scripts/add_hpo_annotation.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence_function.tsv" \
  --hpo_genes_to_phenotype_tsv "./hpo_data/genes_to_phenotype.txt" \
  --output_tsv "gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --log_path "logs/hpo_annotation.log" \
  --log_level "INFO"



python  scripts/extract_high_confidence_testis_sperm_genes.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --output_tsv "high_confidence_testis_sperm_genes_function.tsv" \
  --min_tau 0.95 \
  --min_testis_tpm 5 \
  --min_sperm_tpm 0.1 \
  --log_path "logs/high_confidence_testis_sperm_function.log" \
  --log_level "INFO"



python  scripts/add_proteomics_annotation.py \
  --input_gene_tsv "gtex_v11_testis_ranked_with_sperm_presence_function_hpo.tsv" \
  --proteomics_xlsx "proteomics/SPZ-HomoSapiens.xlsx" \
  --output_tsv "gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --gene_symbol_col "hgnc_symbol_norm" \
  --require_fdr_confidence "High" \
  --restrict_species_name "Homo sapiens" \
  --log_path "logs/proteomics_annotation.log" \
  --log_level "INFO"


python  scripts/extract_high_confidence_testis_sperm_genes.py \
  --input_tsv "gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv" \
  --output_tsv "high_confidence_testis_sperm_genes_function.tsv" \
  --min_tau 0.95 \
  --min_testis_tpm 5 \
  --min_sperm_tpm 0.1 \
  --log_path "logs/high_confidence_testis_sperm_function.log" \
  --log_level "INFO"



python /home/pthorpe001/data/2026_sperm_Gates/scripts/find_male_phenotypic_disease/overlap_testis_specific_with_clinvar_tiers.py \
  --testis_gene_tsv  ~/data/2026_sperm_Gates/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv \
  --clinvar_best_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_best.tsv \
  --clinvar_hc_pathogenic_tsv ~/data/2026_sperm_Gates/clinvar/clinvar_filtered/male_infertility_clinvar_variants_high_confidence_pathogenic.tsv \
  --out_dir ~/data/2026_sperm_Gates/overlap_outputs \
   --gene_symbol_column gene_symbol_norm 