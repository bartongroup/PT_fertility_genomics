#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${HOME}/data/2026_sperm_Gates"



python ~/data/2026_sperm_Gates/PT_fertility_genomics/final_summary/make_master_fertility_results_summary.py \
  --base_dir ~/data/2026_sperm_Gates \
  --testis_run_id testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --out_xlsx ~/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx



# add public proteomics data too .. 
python ~/data/2026_sperm_Gates/PT_fertility_genomics/final_summary/make_master_fertility_results_summary.py  \
  --base_dir ~/data/2026_sperm_Gates \
  --testis_run_id testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --public_proteomics_tsv /home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/PXD014618_proteomics/PXD014618_gene_level_proteomics.tsv \
  --out_xlsx ~/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx



# make plots

python ~/data/2026_sperm_Gates/PT_fertility_genomics/final_summary/make_fertility_gene_set_plots.py \
  --in_xlsx /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/SUMMARY_fertility_evidence.xlsx \
  --sheet_name Tier_Summary_With_Omics \
  --out_dir /home/pthorpe001/data/2026_sperm_Gates/results/FULL_SUMMARY/plots
