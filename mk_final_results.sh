#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${HOME}/data/2026_sperm_Gates"



python ~/data/2026_sperm_Gates/scripts/make_master_fertility_workbook_from_results.py \
  --base_dir ~/data/2026_sperm_Gates \
  --testis_run_id testis_tau0.95_testisTPM5_presentTPM5_spermTPM0.1 \
  --out_xlsx ~/data/2026_sperm_Gates/results/master_workbook/master_fertility_evidence.xlsx
