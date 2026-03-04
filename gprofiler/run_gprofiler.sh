


  # or a folder of tsv gene list

  conda activate R_stats

python run_gprofiler.py \
  --gene_list_tsv sperm_only_genes.tsv \
  --out_dir gprofiler_out/sperm_only_genes \
  --organism hsapiens \
  --gene_column gene_key \
  --single_set_name sperm_only_genes



mkdir -p gprofiler_results

for f in *.tsv; do
  base=$(basename "$f" .tsv)

  python /home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/gprofiler/run_gprofiler.py \
    --gene_list_tsv "$f" \
    --out_dir "gprofiler_results/$base" \
    --organism hsapiens \
    --gene_column gene_key \
    --single_set_name "$base"
done