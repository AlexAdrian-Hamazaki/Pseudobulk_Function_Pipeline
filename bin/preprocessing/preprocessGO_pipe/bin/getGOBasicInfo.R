#!/usr/bin/env Rscript

library(ontologyIndex)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

go_obo_ch <- args[[1]]


ontology <- ontologyIndex::get_OBO(go_obo_ch,
                                  extract_tags = "everything")

columns_to_keep <- c("id", "name", 'namespace', 'def')

df = as.data.frame(ontology)

df_subset <- df %>%
    select(all_of(columns_to_keep))

write.table(df_subset, file = "go_metadata.tsv", row.names = FALSE, quote = FALSE, sep="\t") 
