print("RUNNING EGAD")



packages <-installed.packages()[,"Package"]

if (!"EGAD" %in% packages) {
    BiocManager::install("EGAD")
    }

library(EGAD)
library(tidyverse)
library(stringr)
#data_file <-snakemake@input[[1]]

#Change
save <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/EAGD/EGAD_sum_pc_OPfiltered.csv"
# data_file <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/pseudobulk/sum_pseudobulk.csv"
#data_file <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/bulk/bulk_pc.csv"

is_bulk <- TRUE

#Never change
pc_genes <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/pc_genes/processed_uniprot.csv"
sample_names <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/sample_names/bulk_pseudo.csv"
shared_genes <- "~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/pc_genes/scAndBulkOverlapGenes.txt"

#../../../../pipeline42/datasets/Gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct



################ 2.1 : MAKING DATA SETS

########### MAKE COEXPRESSION NETWORK

#~~~~~~~~~ If using file
#expression_data <- t(read.delim('data/LiverAverageCounts.csv', sep = ',', header= TRUE, row.names = 1))


#~~~~~~~~ Diagnostic WD

#getwd()

#~~~~~~~~~If using command line

#args = commandArgs(trailingOnly=TRUE)
#data_file <- args[1]


print("Loading Expression Dataset")

expression_data <- (read.delim(data_file, sep = ',', header= TRUE, row.names = 1))
sample_names_df <- read.delim(sample_names, sep = ",", row.names = 1, header = TRUE)




<<<<<<< HEAD
### Remove Organism Parts that are NOT shared between the two Tabula And Gtex Datasets
filter_by_OP <- function(OP, expression_data) {
  expression2 <- expression_data %>%
    filter(str_detect(rownames(expression_data), paste0(OP,".*")))
  return (expression2)
=======
coexpression_network[is.na(coexpression_network)] <- 0
>>>>>>> 3c1fb065b25c7ca5703de1a68e8cd379c6c7289b

}
if (is_bulk) {
  print('Subsampling for Bulk name Organism Parts')
  OP_names <- sample_names_df$bulk_names
  expression_data <- expression_data %>%
    filter(rownames(expression_data) %in% OP_names)
} else {
  print('Subsampling for Pseudobulk name Organism Parts')

  OP_names <- sample_names_df$pseudo_names


  expression_data_2 <- lapply(OP_names, filter_by_OP, expression_data)
  
  expression_data <- do.call(rbind, expression_data_2)
  }


print("Filtering Expression Data for PC genes")
pc <- read.delim(pc_genes, sep = ",", header = TRUE, row.names = 1)
pc_names <- pc$FirstUniprot
pc_expression <- expression_data[,colnames(expression_data) %in% pc_names]
print(paste("Removed",ncol(expression_data)- ncol(pc_expression) , "non-protein coding genes"))




######### Build Coexpression Network

coexpression_network <- cor(pc_expression)
coexpression_network[is.na(coexpression_network)] <- 0


############ BUILDING ANNOTATION SET
print("Building Annotation Set")

### With builtin GO
#annotations <- make_annotations(GO.human[,c('GO', 'evidence')], unique(GO.human$GO), unique(GO.human$evidence))

### With Custom GO with BP
GO <- read.delim(file = '~/Masters/Pseudobulk_Function_Pipeline_HighRes/data/GO/pro_GO.csv', sep = ",", stringsAsFactors = TRUE)

# Filter the GO to Gene pairings for only Genes measured in our expression data because we only want GO terms with 20>= genes
expression_genes <- colnames(expression_data)
GO <- filter(GO, GO$DB_Object_Symbol %in% expression_genes) 
#in our sc data this removes ~3,000 genes
#in bulk this removes ~16000 genes sheesh

# Filter for only the genes that are measured in both data types (note this is redudant with the prev step now)
sharedGenes <- read_csv(file =shared_genes, col_names = FALSE)
GO <- filter(GO, GO$DB_Object_Symbol %in%sharedGenes$X1)

GO_unique <- data.frame(table(GO$GO.ID))
colnames(GO_unique) <- c('GO', 'count')


# Create a histogram looking at how many genes are affiliated with each GO term
ggplot(data = GO_unique)+ 
  geom_histogram(mapping = aes(count)) +
  scale_y_continuous(trans = 'log10') +
  labs(title = 'Distribution of Genes in GO Terms') + 
  xlab('Number of Genes')+
  ylab('GO Terms with number of genes')
  

################ Remove GO Terms with less than 20 Genes in the expression data.

GO_unique_filtered <-  filter(GO_unique, count >=20)

1-nrow(GO_unique_filtered)/nrow(GO_unique) #92.8 % of GO terms were removed in sc. 92% in bulk

ggplot(data = GO_unique_filtered)+ 
  geom_histogram(mapping = aes(count), breaks = c(0, 19, 30, 40, 50, 60, 70, 80, 90, 100, 150)) + 
  scale_y_continuous(trans = 'log10') +
  labs(title = 'Distribution of Genes Assosiated with GO Terms') + 
  xlab('Number of Genes')+
  ylab('Count of GO Terms')

# With GO_unique_filtered, we now have all of the GO Terms we want to use in our analysis

GO_20_or_more <-dplyr::filter(GO, (GO$GO.ID %in% GO_unique_filtered$GO))

# 50694 GO To Gene assosiations were filtered out
1-nrow(GO_20_or_more)/nrow(GO) #44.9% of Gene to Go Term assosiations were for GO terms with less than 20 genes.  45% in bulk

#Note: We removed 86.2% of GO terms, this onyl removed 31.5% of GO to gene assosiations


#Make one hot encoding matrix
# Contains only GO terms with 20 genes or more that were measured in both datasets
annotations <- make_annotations(GO_20_or_more[,c('DB_Object_Symbol', 'GO.ID')],  unique(GO_20_or_more$DB_Object_Symbol), unique(GO_20_or_more$GO.ID))

################ Neighbor Voting
print("Performing Neighbor Voting. This can take a while")
auroc <- neighbor_voting(genes.labels = annotations, 
                          network = coexpression_network,
                          nFold = 3,
                          output = "AUROC")


#auroc <- run_GBA(network = coexpression_network,
#                labels = annotations)

rm(coexpression_network )
rm(annotations)


print(paste0("Wrote AURUC to ", save ))

write.table(x = auroc, paste0(save), sep = ",")



