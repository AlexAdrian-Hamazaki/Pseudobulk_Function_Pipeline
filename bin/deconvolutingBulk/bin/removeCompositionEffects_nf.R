#!/usr/bin/env Rscript

if (!require(markerGeneProfile)) devtools::install_github('oganm/markerGeneProfile')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(gplots)) install.packages('gplots')
if (!require(viridis)) install.packages('viridis')
if (!require(dplyr)) install.packages('dplyr')
if (!require(knitr)) install.packages('knitr')

library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
tissue_path = args[1]
tissue_name = args[2]
ct_marker_path = args[3]
ct_marker_name = args[4]

# tissue_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/deconvolutingBulk/data/dev/july/tissues/Blood_gtex.csv"
# ct_marker_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/deconvolutingBulk/data/dev/july/markers/pbmc_sc_with_metadata_pc_cpm_markers.csv"

################## Loading Datasets #######################
# Load bulk Dataset. Coerce it so Rows are genes, cols are samples.
# but the column "Gene.Symbol" is where the gene information is
tissue_df = read.csv(tissue_path, row.names = 1, header = TRUE)
tissue_df  = as.data.frame(t(tissue_df ))
tissue_df  = rownames_to_column(tissue_df , var = 'Gene.Symbol') # this should be a dataframe where the first column is Gene Symbol and the other columns are samples


# load marker genes 
ct_markers = read.csv(ct_marker_path, header = TRUE)
# Coerce the brain_markers into a list of character vectors, where each character vector is the gene symbols of the marker genes for a cell type
# Split the data frame by 'group' column
lo_lo_CT_markers <- split(ct_markers$names, ct_markers$group) # Splits data into multiple dfs
# lo_lo_CT_markers <- lapply(names(grouped_data), function(ct) {
# # vector of genes assosiated with a cell type (group)
#   genes_vector <- grouped_data[[ct]]
#   return (genes_vector)
# })
# names(lo_lo_CT_markers) <- names(grouped_data)



################## Getting Marker Gene Profiles for each cell type we have #######################

# Estimate the first PC for list of cell type markers.
estimations =  mgpEstimate(exprData=tissue_df, #  is the expression matrix which should include gene identifiers as a column.
                           genes=lo_lo_CT_markers,
                           geneColName='Gene.Symbol',
                           geneTransform = NULL) # removes genes if they are the minority in terms of rotation sign from estimation process
estimations$main <- estimations$estimates %>% sapply(function(vals) { vals }) 

# Each vector in estimations$estimate is the first PC for the ct markers of a cell type. We still have the same number of instances because of PCA we are reducing the rank (genes) not the samples
# So each vector in estimations$estimate is of length equal to the number of samples we have. And basically each sample is given a numerical value about where it placed in its variance in the first PC
# The first PC is taken from only the expresion of the ct marker genes. It represents the indepdent axis that captures the most variance across those genes. It is previosly know that this has to do with cell type composition
# So now, we have principle components that represent variance in expression of cell type marker genes for all of our cell types
# Now we can regress out these effects, to remove bulk sample's variance in expression affiliated with these PCs (composition variance)

######################## Regress out the MGP Effects ######################

# First define a function that fits a linear model for each gene based on the MGPs for our cell types
fit_model_for_one_gene = function(currGene, tissue_df, estimations) {
    # First, get the bulk expression for the gene I'm currently fitting a model for
    expr_gene_in_bulk  = as.data.frame(t(tissue_df[currGene,]))
    colnames(expr_gene_in_bulk) = c(currGene)
    #Make a dataframe of all of the covariates (1rst PC of each cell type marker gene)
    covariates <- estimations$main[row.names(expr_gene_in_bulk), ] %>% as.data.frame() # ensure samples match exprmat
    # Combine the dataframes to have Y and ~ X together
    merged_df = merge(covariates, expr_gene_in_bulk, by = 'row.names') %>% 
        column_to_rownames(var = "Row.names")

    # Make model that tries to fit bulk expression's expression for this gene as a function of the composition of cell types in that bulk sample.
    # Here's the expectation. For example, if we have a cell type marker gene for neurons, we expect that that gene's bulk expression is
    # influenced by the proportion of neurons in each sample. Samples with lots of neurons will have higher values for that gene
    # So because our PCs roughly tell us how much of each cell type there is a sample (highest source of variance across cell type marker genes is percent composition),
    # We can fit a model that learns coefficient about how each gene's expression is a conseuqence of merely cell type compositions.
    # The residuals of that model would be the gene's expression that is not mrerly a consequence of cell type composition
    formula_str = paste(currGene, "~ .")
    model = lm(formula_str, merged_df)
    return(model)
}

lo_genes = tissue_df[["Gene.Symbol"]] # We will be fitting a model for each gene, so lets get their names
rownames(tissue_df) = tissue_df$Gene.Symbol # Set the gene symbol as the index for easy filtering in the lm
tissue_df = tissue_df %>% select(-Gene.Symbol) # Remove the gene symbol colm from the tissue _df to avoid complications in lm

##### Fit linear models for each gene that fot a gene's bulk expression as a consequence of the relative abundances of our cll types
lo_linear_models = lapply(lo_genes, fit_model_for_one_gene, tissue_df, estimations)
##### Extract residuals from our models. These models represent how much expression from bulk is NOT a consequence of cell type composition differences
exprmatRes <- lo_linear_models  %>% sapply(function(currLm) { currLm$residuals } ) %>% t() # extract residuals
# Set rows to gene names
row.names(exprmatRes) = lo_genes
# write
write.csv(exprmatRes, file = paste0(tissue_name,"_",ct_marker_name,"_regressed.csv"))
# write.csv(exprmatRes, file = "july_test_markers.csv")


#### Sanity Checks

# First, the expression of cell type marker genes should be much smaller in the residuals, but they'd have to be put on the same scale. For ct marker genes, if we plot the resids across samples and the real bulk across samples we should see less variance
#

# expr_filtered <- as.data.frame(exprmatRes) %>% 
#     rownames_to_column(var = 'genes') %>%
#     filter(genes %in% ct_markers$names) %>%
#     column_to_rownames('genes')
# dim(expr_filtered)
# tissue_filtered <- tissue_df %>%
#     rownames_to_column(var = 'genes') %>%
#     filter(genes %in% ct_markers$names)%>%
#     column_to_rownames('genes')
# dim(tissue_filtered)

# plot(x = 1:2642 ,y = expr_filtered[1,])
# plot(x = 1:2642 ,y = tissue_filtered[1,])

# # Calculate standard deviation of each row
# row_sds_expr <- apply(expr_filtered, MARGIN = 1, FUN = sd)
# row_sds_tissue <- apply(tissue_filtered, MARGIN = 1, FUN = sd)

# plot(row_sds_tissue, row_sds_expr)

