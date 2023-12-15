#!/usr/bin/env Rscript


################ 1.1 : Load Inputs


# set the expression matrix path
args <- commandArgs(trailingOnly=TRUE)
expression_matrix_path <- args[[1]]
expression_matrix_name <- args[[2]]
GO_annot_path <- args[[3]]
GO_annot_path_name <- args[[4]]

# The column you need to select in the gene annotations based on
# what the gene names are in the expression matrix g
GO_annot_gene_column <- args[[5]]
bootstrap <- args[[6]]

# expression_matrix_path <- "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkEGADPipeline/data/splitOPs1/splits/subsamples/Blood_split.csv_1.csv.gz"
#expression_matrix_path <- "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkSimulations/data/simulated_bulk_dataset/brain_sc_with_metadata_cpm_pc.h5ad/simulated_ss100_var0.01_nsim10.csv"
#GO_annot_path <- "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/preprocessing/preprocessGO_pipe/data/GO_annotationsWithENSGandPC/bp_annotations_withGeneData.csv"
#GO_annot_gene_column <- "ensembl_gene_id"
#GO_annot_gene_column <- "hgnc_symbol"
#GO_annot_gene_column <- "DB_Object_Symbol"

################ 1.1 : Load Packages
library(EGAD)
library(tidyverse)
library(stringr)


################ 2.1 : MAKING DATA SETS

########### MAKE COEXPRESSION NETWORK

print("Loading Expression Dataset")
expression_data <- read.table(expression_matrix_path, sep = ',', header= TRUE)

# Transpose the data frame such that columns are Genes and rows are samples
#expression_data <- t(expression_data)

expression_data <- expression_data %>%
  select(-c('SAMPID'))

# Convert values to numeric while preserving NAs
#expression_data <-  apply(expression_data, 2, function(col) as.numeric(as.character(col)))

######### Build Coexpression Network

print("Building Coexpression Network")

coexpression_network <- cor(expression_data)
coexpression_network[is.na(coexpression_network)] <- 0

head(coexpression_network)
dim(coexpression_network)
############ BUILDING ANNOTATION SET
print("Building Annotation Set")

### Load GO annotations
GO_annotations <- read.delim(file = GO_annot_path, sep = ",", stringsAsFactors = TRUE, row.names=NULL)

dim(GO_annotations)

usingENSG<- function(expression_data) {
  ###
  # Returns true if all the gene symbols are ENSG
  ###
  allStartsWithENSG <- all(grepl("^ENSG", colnames(df)))
  return(allStartsWithENSG)
}


filterGOAnnotations <- function(GO_annotations, expression_data, GO_annot_gene_column) {
  # We only want GO terms with 20>= MEASURED GENES
  # Filter the GO_annotations, for only Genes that are in our expression data
  
  # isENSG: -> bool. True if using ENSG gene names

  # Returns a dataframe of GO annotations with only genes that are in our exapression dataset
  expression_genes <- colnames(expression_data)

  #print(GO_annotations[, GO_annot_gene_column] %in% expression_genes)
  print(GO_annot_gene_column)
  print(head(expression_genes))
  GO_annotations <- filter(GO_annotations, GO_annotations[, GO_annot_gene_column] %in% expression_genes)

  return (GO_annotations)
}

keep20PlusGOAnnotations <- function(GO_annotations, expression_data, GO_annot_gene_column){
  ###
  # We only want GO terms with 20>= MEASURED GENES
  # Filter out GO terms that contain less than 20 genes
  # Returns a dataframe of GO annotations that have 20 or more genes
  ###
  
  # Check if the expression data use ENSG or Gene symbols
  #isENSG <- usingENSG(expression_data)
  

  # Get a GO annotatino dataset containing only genes that are expressed
  GO_annotations <- filterGOAnnotations(GO_annotations = GO_annotations,
                                        expression_data = expression_data,
                                        GO_annot_gene_column = GO_annot_gene_column
  )
  # Remove rows where the GO_annot_gene_column is not unique
  GO_annotations_grouped = GO_annotations %>%
    group_by(GO.ID) %>%
    distinct(DB_Object_Symbol, .keep_all = TRUE) ## TODO Use GO annot gene column?

  # Get the size of each GO
  group_sizes <- GO_annotations_grouped  %>%
    summarise(n = n())
  print(group_sizes)

  # Select only GO Terms with more or equal to 20 genes 
  large_GO_Terms <- group_sizes %>%
    filter(n >= 20)
  print(large_GO_Terms)

  # Filter the measured GO annotations for only those with 20 or more genes
  GO_annotations_large = GO_annotations %>%
  filter(GO.ID %in% large_GO_Terms$GO.ID)


  # Remove rows with duplicate values in col1 and col2
  df_no_duplicates <- GO_annotations_large[!duplicated(GO_annotations_large[c("DB_Object_Symbol", "GO.ID")]), ]

  # stopifnot(length(unique(GO_annotations_large$GO.ID)) == nrow(large_GO_Terms))

  return (df_no_duplicates)
}

# Remove GO annotations that don't have at least 20 Genes in our expression_data matrix.
# Note the expression_data matrix should only have PC genes at this point
GO_annotations20 <- keep20PlusGOAnnotations(GO_annotations = GO_annotations,
                                          expression_data = expression_data,
                                          GO_annot_gene_column = GO_annot_gene_column)
### Write summary statistics to file ###
nrow(GO_annotations)
nrow(GO_annotations) - nrow(GO_annotations20) # The number of GO annotations we removed

# The average amount of Genes in a GO term

getMeanCount <- function(GO_annotations) {
  ###
  # Returns the average amount of genes in a GO id as an int
  ###
  grouped_df <- GO_annotations %>%
    group_by(GO.ID) %>%
    summarise(count = n())

  return (mean(grouped_df$count))
}


# TODO More summary stats
# A graph of histogram before and after
# The GO terms with the most genes in them that wre removed because the genes were not measured

###Make one hot encoding matrix
colnames(coexpression_network)

colnames(GO_annotations20)
annotations <- make_annotations(GO_annotations20[,c(GO_annot_gene_column, 'GO.ID')], unique(GO_annotations20[,GO_annot_gene_column]), unique(GO_annotations20$GO.ID))

annotations[1:10, 1:10]
head(coexpression_network)
sum(annotations)

colnames(annotations)

#coexpression_network <- coexpression_network[1:1000, 1:1000]
#annotations <- annotations[1:100, 1:300]

dim(coexpression_network)
dim(annotations)

################ Neighbor Voting
print("Performing Neighbor Voting. This can take a while")
auroc <- neighbor_voting(genes.labels = annotations,
                         network = coexpression_network,
                         nFold = 3,
                         output = "AUROC")

rm(coexpression_network)
rm(annotations)


########### Write EGAD results
# args[3] is organism part name
# args[4] is MF or BP name
write.table(x = auroc, paste0(expression_matrix_name,"_",GO_annot_path_name,"_",bootstrap,"_EGAD.csv"), sep = ",")

