



# build the features

python build_gene_context_table.py \
  --gene_list_tsv sperm_genes.tsv \
  --gencode_gtf gencode.vXX.annotation.chr.gtf.gz \
  --hg38_fasta hg38.fa.gz \
  --rmsk_tsv_gz rmsk.txt.gz \
  --segdup_tsv_gz genomicSuperDups.txt.gz \
  --recomb_bw recombAvg.bw \
  --out_tsv gene_context_features.tsv \
  --mapping_report_tsv gene_mapping_report.tsv



# Permutation tests (10,000 matched sets), matching on chromosome + gene length bins, and (optionally) GC bins:

python matched_set_enrichment.py \
  --features_tsv gene_context_features.tsv \
  --target_gene_list_tsv sperm_genes.tsv \
  --gene_key_column gene_key \
  --n_permutations 10000 \
  --n_length_bins 10 \
  --match_on_gc \
  --n_gc_bins 10 \
  --out_prefix sperm_vs_null

