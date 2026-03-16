

script_dir="/home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/Genomic_interogation/"

# build the features


python ${script_dir}/build_gene_universe_from_gtf.py \
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


# add Open Targets tractability annotations to the gene universe file
  python /home/pthorpe001/data/2026_sperm_Gates/PT_fertility_genomics/Genomic_interogation/annotate_opentargets_tractability.py \
  --in_tsv ~/data/2026_sperm_Gates/genome_resources/gene_context_features_universe.tsv \
  --ensembl_column gene_id \
  --out_tsv ~/data/2026_sperm_Gates/genome_resources/gene_context_features_universe_plus_tractability.tsv \
  --verbose

# This will run the enrichment analysis on the chr the genes are on and nearest nieghbours too. 

 python split_gene_lists_wide_to_tsv.py  --wide_tsv final_lists.txt --out_dir gene_lists_split

for f in gene_lists_split/*_genes.tsv; do
  base=$(basename "$f" _genes.tsv)
  python genome_organisation_tests.py \
    --features_tsv gene_context_features.tsv \
    --target_gene_list_tsv "$f" \
    --gene_key_column gene_key \
    --out_prefix "genome_org_${base}" \
    --n_permutations 10000 \
    --n_length_bins 10 \
    --match_on_gc \
    --n_gc_bins 10
done

python genome_organisation_tests.py \
  --features_tsv gene_context_features.tsv \
  --target_gene_list_tsv sperm_genes.tsv \
  --gene_key_column gene_key \
  --out_prefix sperm_genes \
  --n_permutations 10000 \
  --n_length_bins 10 \
  --match_on_gc \
  --n_gc_bins 10




# create summary of significant results for nearest neighbour test
  echo -e "gene_list\tn_targets\tn_permutations\tmatch_on_gc\tobserved_median_nn_bp\tobserved_mean_nn_bp\tp_median\tp_mean\tpattern_median" \
  > significant_nearest_neighbour.tsv

awk -F'\t' 'NR>1 && $9 < 0.05 {
    name = FILENAME
    gsub("genome_org_|\\.nearest_neighbour_summary.tsv","",name)

    pattern = ($7 < $8) ? "clustered" : "dispersed"

    print name "\t" $1 "\t" $2 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" pattern
}' genome_org_*.nearest_neighbour_summary.tsv >> significant_nearest_neighbour.tsv



# Create summary of significant results for chromosome enrichment test
echo -e "gene_list\tobserved_distance\tnull_mean_distance\tp_value\tpattern" > nearest_neighbour_significant.tsv

awk -F'\t' 'NR>1 && $4 < 0.05 {
    pattern = ($2 < $3) ? "clustered" : "dispersed"
    name = FILENAME
    gsub("genome_org_|\\.nearest_neighbour_summary.tsv","",name)
    print name "\t" $2 "\t" $3 "\t" $4 "\t" pattern
}' genome_org_*.nearest_neighbour_summary.tsv >> nearest_neighbour_significant.tsv