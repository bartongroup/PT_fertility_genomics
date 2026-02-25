#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages(library(biomaRt))

args <- commandArgs(trailingOnly = TRUE)
out_path <- "hsapiens_protein_coding_background.tsv"
if (length(args) >= 2 && args[1] == "--out") {
  out_path <- args[2]
}

mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

df <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = mart
)

genes <- unique(df$hgnc_symbol)
genes <- genes[!is.na(genes) & genes != ""]

write.table(
  data.frame(gene = genes),
  file = out_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message(sprintf("Wrote %d protein-coding HGNC symbols to %s", length(genes), out_path))
