python replace_expression_ids_with_gene_names.py \
    --expression_tsv gtex_v11_tissue_medians_trinity.tsv \
    --mapping_tsv gtex_v11_tissue_medians_gene_names.tsv.mapping_table.tsv \
    --output_tsv gtex_v11_tissue_medians_gene_symbols.tsv


python replace_expression_ids.py \
  --expression_tsv gtex_v11_tissue_medians_trinity.tsv \
  --mapping_tsv gtex_v11_tissue_medians_gene_names.tsv.mapping_table.tsv \
  --output_tsv gtex_v11_tissue_medians_gene_symbols.tsv



head -1 gtex_v11_tissue_medians_gene_symbols_renamed.tsv  | tr '\t' '\n' | tail -n +2 | awk '{print $1"\t"$1}' > samples_described.txt


conda activate R

 perl  ~/apps//trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/PtR -s samples_described.txt -m prioritised_expression.txt --save --heatmap --log2 --center_rows





######### REDUCED COLOUMNS !!!!!!!!

cut -f 1,2,3,4,5,50- gtex_v11_tissue_medians_gene_symbols_renamed.tsv > tmp.txt
cat tmp.txt | grep -Ff  prioritised.txt >> tmp.txt


head -1 tmp.txt  | tr '\t' '\n' | tail -n +2 | awk '{print $1"\t"$1}' > samples_described.txt
 
head -n 1 tmp.txt   | tr '\t' '\n' | tail -n +2 | awk '{print $1"\t"$1}' > samples_described.txt
 
 perl  ~/apps//trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/PtR -s samples_described.txt -m prioritised_expression.txt --save --heatmap \ 
      --log2 --center_rows --heatmap_scale_limits "-5,5"   --top_genes 100


 perl  ~/apps//trinityrnaseq-v2.15.1/Analysis/DifferentialExpression/PtR -s samples_described.txt -m prioritised_expression.txt --save --heatmap \
      --log2 --center_rows --heatmap_scale_limits "-5,5" \
--Zscale_rows
