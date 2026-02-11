#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${HOME}/data/2026_sperm_Gates"


python ~/data/2026_sperm_Gates/scripts/make_master_fertility_workbook.py \
  --base_dir ~/data/2026_sperm_Gates \
  --testis_run_id testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --male_infertility_run_relpath hpo_data/male_infertility_out_phrase_seeds \
  --clinvar_results_relpath results/05_clinvar_high_confidence \
  --out_xlsx ~/data/2026_sperm_Gates/results/master_workbook/master_fertility_genes.xlsx \
  --include_variant_sheets

