

python run_gprofiler.py \
  --excel /path/to/sperm_gene_sets.xlsx \
  --out_dir /path/to/gprofiler_out \
  --organism hsapiens

  # or a folder of tsv gene list

mkdir -p gprofiler_results

for f in *.tsv; do
    base=$(basename "$f" .tsv)

    python - <<PY
import pandas as pd
df = pd.read_csv("$f", sep="\t")
df.to_excel("${base}.xlsx", index=False)
PY

    python /home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/gprofiler/run_gprofiler.py \
        --excel "${base}.xlsx" \
        --out_dir gprofiler_results/"$base" \
        --organism hsapiens

    rm "${base}.xlsx"
done