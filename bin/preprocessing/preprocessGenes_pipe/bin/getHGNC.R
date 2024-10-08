#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", version = 3.17)
}

if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

library("biomaRt")
library("tidyverse")

refseqPath <- commandArgs(trailingOnly = TRUE)[1]
refseqPath <-  "refseq_pcoding_human.tsv"
refseq <-  read.csv(refseqPath, header = TRUE, sep = "\t")

removeDecimal <- function(vec) {
  modified <- sub("\\..*", "", vec)
  return(modified)
}

ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

values<- removeDecimal(refseq$Refseq_ID)

master <-  getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = values, mart= ensembl)

write_delim(master, file = "humanPCMaster.tsv", delim = "\t")


