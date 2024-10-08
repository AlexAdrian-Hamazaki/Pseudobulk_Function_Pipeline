#!/usr/bin/env Rscript


if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", version = 3.17)
}

if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

library(biomaRt)
library(tidyverse)

mart <- useMart("ENSEMBL_MART_ENSEMBL")


mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

# Get all gene ENGS and gene symbols of relevance
gene_info <- getBM(attributes = c("ensembl_gene_id",  'hgnc_symbol', 'uniprotswissprot', 'entrezgene_id',  'description'),
                   mart = mart)


# Filter for only PC genes
pc_genes <- gene_info %>%
  filter('hgnc_symbol' != "")


write_delim(gene_info, file = "humanAllGenes.csv", delim = ",")
write_delim(pc_genes, file = "humanPCGenes.csv", delim = ",")
