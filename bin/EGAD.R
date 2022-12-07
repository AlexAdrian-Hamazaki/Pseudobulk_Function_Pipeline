
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


packages <-installed.packages()[,"Package"]

if (!"EGAD" %in% packages) {
    BiocManager::install("EGAD")
    }

library(EGAD)
data_file <-snakemake@input[[1]]
print(data_file)
save <- "data/EAGD/EGAD.csv"


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
expression_data <- t(read.delim(data_file, sep = ',', header= TRUE, row.names = 1))




coexpression_network <- build_coexp_network(expression_data, 
                                            gene.list = rownames(expression_data),
                                            method = 'spearman')

#coexpression_network[is.na(coexpression_network)] <- 0





############ BUILDING EDGE DATA FRAME

# build_edge_df<- function(colname, coexpression_network){
#   
#   interactions <- rownames(coexpression_network)[coexpression_network[,colname]>0]
#   
#   first_gene <- rep(colname, length(interactions))
#   
#   edges_df<- data.frame(gene1 = first_gene,
#                        gene2 = interactions)
#   
#   
#   return (edges_df)
# }

# 
#   
# edges <- bind_rows((lapply(colnames(coexpression_network), build_edge_df, coexpression_network)))


############ BUILDING ANNOTATION SET
print("Building Annotation Set")
annotations <- make_annotations(GO.human[,c('GO', 'evidence')], unique(GO.human$GO), unique(GO.human$evidence))
                       
#TODO GET UP TO DATE GO - > ENTREZ ID ONTOLOGY/ GENE SYMBOL ONTOLOGY

# any(GO.human$entrezID == "653635")
# rownames(expression_data) %in% GO.human

################ Neighbor Voting
print("Performing Neighbor Voting")
auroc <- neighbor_voting(genes.labels = annotations, 
                         network = coexpression_network,
                         nFold = 3,
                         output = "AUROC")

rm(coexpression_network )
rm(annotations)


print(paste0("Wrote AURUC to ", save ))

write.table(x = auroc, paste0(save), sep = ",")


