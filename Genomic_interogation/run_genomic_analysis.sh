

script_dir="/home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/Genomic_interogation/"

# build the features


python build_gene_universe_from_gtf.py \
  --gencode_gtf gencode.v49.primary_assembly.annotation.gtf.gz \
  --out_tsv gencode_v49_gene_keys_universe.tsv \
  --verbose

# gzip of the  human fa will not work directly with pyfaidx 
python ${script_dir}/build_gene_context_table.py \
  --gene_list_tsv gencode_v49_gene_keys_universe.tsv  \
  --gencode_gtf gencode.v49.primary_assembly.annotation.gtf.gz \
  --hg38_fasta hg38.fa \
  --rmsk_tsv_gz rmsk.txt.gz \
  --segdup_tsv_gz genomicSuperDups.txt.gz \
  --recomb_bw recombAvg.bw \
  --out_tsv gene_context_features.tsv \
  --mapping_report_tsv gene_mapping_report.tsv



# Permutation tests (10,000 matched sets), matching on chromosome + gene length bins, and (optionally) GC bins:

python ${script_dir}/matched_set_enrichment.py \
  --features_tsv gene_context_features.tsv \
  --target_gene_list_tsv sperm_genes.tsv \
  --gene_key_column gene_key \
  --n_permutations 10000 \
  --n_length_bins 10 \
  --match_on_gc \
  --n_gc_bins 10 \
  --out_prefix sperm_vs_null

