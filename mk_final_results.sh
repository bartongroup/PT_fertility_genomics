#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${HOME}/data/2026_sperm_Gates"


SCRIPTS_DIR="${BASE_DIR}/PT_fertility_genomics/"


testis_annotated_tsv="${base_dir}/results/${testis_run_id}/08_proteomics_annotation/gtex_v11_testis_ranked_with_sperm_presence_function_hpo_prot.tsv"
testis_annotated_backup="${testis_annotated_tsv%.tsv}.pre_accessibility_backup.tsv"
testis_annotated_access="${testis_annotated_tsv%.tsv}.with_accessibility.tsv"

# Inputs for accessibility annotation (you will set these)
gene_to_uniprot_tsv="${resources_dir}/gene_to_uniprot.tsv"
uniprot_annotation_tsv="${resources_dir}/uniprot_annotations.tsv"

# 1) Backup once (safe)
if [ ! -s "${testis_annotated_backup}" ]; then
  cp "${testis_annotated_tsv}" "${testis_annotated_backup}"
fi

# 2) Produce annotated version
python "${script_dir}/annotate_sperm_biochemical_accessibility.py" \
  --in_tsv "${testis_annotated_backup}" \
  --out_tsv "${testis_annotated_access}" \
  --gene_column "gene_key" \
  --gene_to_uniprot_tsv "${gene_to_uniprot_tsv}" \
  --uniprot_annotation_tsv "${uniprot_annotation_tsv}" \
  --verbose

# 3) Replace the original path with the annotated file (so workbook builder sees it)
cp "${testis_annotated_access}" "${testis_annotated_tsv}"

# 4) Run your workbook builder as normal (unchanged)
python "${script_dir}/make_master_fertility_workbook_from_results.py" \
  --base_dir "${base_dir}" \
  --testis_run_id "${testis_run_id}" \
  --out_xlsx "${out_xlsx}" \
  ${include_variant_sheets_flag} \
  ${public_proteomics_flag} \
  ${no_internal_proteomics_flag} \
  ${literature_flags}



# add public proteomics data too .. 
python ~/data/2026_sperm_Gates/PT_fertility_genomics/final_summary/make_master_fertility_results_summary.py  \
  --base_dir ~/data/2026_sperm_Gates \
  --testis_run_id testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --public_proteomics_tsv /home/pthorpe001/data/2026_sperm_Gates/PXD037531/PXD037531_gene_level_proteomics.tsv \
  --out_xlsx ~/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx



# make plots

python ~/data/2026_sperm_Gates/PT_fertility_genomics/final_summary/make_fertility_gene_set_plots.py \
  --in_xlsx /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx \
  --sheet_name Tier_Summary_With_Omics \
  --out_dir /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/plots


# add biochem properties

python ~/data/2026_sperm_Gates/PT_fertility_genomics/annotate_sperm_biochemical_accessibility.py \
  --excel_in /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx \
  --excel_out /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.biochem.xlsx \
  --uniprot_annotation_tsv /home/pthorpe001/data/2026_sperm_Gates/genome_resources/uniprot/uniprotkb_reviewed_true_2026_02_25.tsv \
  --online_uniprot_mapping \
  --verbose
