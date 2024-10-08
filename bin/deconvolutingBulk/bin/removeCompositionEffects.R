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

# load brain bulk dataset
brain_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/bulkEGADPipeline/data/splitOPs1/splits/Brain_split.csv.gz"
brain_df = read.csv(brain_path, row.names = 1, header = TRUE)
brain_df = as.data.frame(t(brain_df))
# make the index be integrers and the Gene.Symbol name be where all the gene information is
brain_df = rownames_to_column(brain_df, var = 'Gene.Symbol') # this should be a dataframe where the first column is Gene Symbol and the other columns are samples
head(brain_df)

# load marker genes /space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/deconvolutingBulk/data/brain_marker_genes.csv
#
marker_brain_path = "/space/grp/aadrian/Pseudobulk_Function_Pipeline_HighRes/bin/deconvolutingBulk/bin/notebooks/pbmc_markers.csv"
brain_markers = read.csv(marker_brain_path, header = TRUE, row.names = 1)
brain_markers_short = head(brain_markers, 10)
ncol(brain_markers_short)
brain_markers_short
grab_cols =  seq(2, nrow(brain_markers_short), by = 2)
grab_cols
brain_markers_col = brain_markers_short[,grab_cols]
brain_markers_col
# Coerce the brain_markers into a list of character vectors, where each character vector is the gene symbols of the marker genes for a cell type
lo_lo_CT_markers = lapply(brain_markers_col, as_vector)
lo_lo_CT_markers
# Check if lo_lo_hgnc_markers is a matrix or data frame
lo_lo_CT_markers = as.list(brain_markers_col)
lo_lo_CT_markers
dim(brain_df)
brain_df
brain_df$Gene.Symbol
typeof(lo_lo_CT_markers)
# Estimate the Marker Gene Proviles.
estimations =  mgpEstimate(exprData=brain_df,
                           genes=lo_lo_CT_markers,
                           geneColName='Gene.Symbol',
                           geneTransform = NULL) # removes genes if they are the minority in terms of rotation sign from estimation process
estimations$estimates
head(brain_df$Gene.Symbol)
estimations$main <- estimations$estimates %>% sapply(function(vals) { vals }) 

# Graph
for (ct in colnames(estimations$main)) {
    # Get the data for this cell type
    mgp = estimations$main[,ct]
    mgp_df = data.frame(mgp)
    # brain_markers
    ct_markers = brain_markers_short[,ct]

    ct_markers_expression = brain_df %>%
        filter(Gene.Symbol %in% ct_markers ) %>%
        column_to_rownames("Gene.Symbol")

    ct_markers_expression = as.data.frame(t(ct_markers_expression))  %>%
        rownames_to_column('gtex_id') %>%
        melt(id.vars = 'gtex_id')


    # print(head(ct_markers_expression))
    # print(head(data.frame(mgp)))
    merged_melted = merge(ct_markers_expression, mgp_df,  by.x = "gtex_id", by.y = "row.names")
    # print(head())

    faceted_plot = ggplot(merged_melted, aes(x = mgp , y = value)) + 
        geom_point() +
        facet_wrap(~variable, scales = "free_y") +
        labs(title = paste(ct, "mgp and marker gene correlations"))
    ggsave(paste0(ct, "_mgp_correlation.png"), plot = faceted_plot, units = "in")
}

########################
lo_genes = brain_df[,1]
rownames(brain_df) = brain_df$Gene.Symbol
brain_df = brain_df %>% select(-Gene.Symbol)

# Fit lm for each gene onto MGP
fit_model_for_one_gene = function(currGene, brain_df, estimations) {

    # First, get the bulk expression for the gene I'm currently fitting a model for
    expr_gene_in_bulk  = as.data.frame(t(brain_df[currGene,]))
    colnames(expr_gene_in_bulk) = c("currGene")
    #Make a dataframe of all of the covariates (1rst PC of each cell type marker gene)
    covariates <- estimations$main[row.names(expr_gene_in_bulk), ] %>% as.data.frame() # ensure samples match exprmat

    merged_df = merge(covariates, expr_gene_in_bulk, by = 'row.names') %>% 
        column_to_rownames(var = "Row.names")

    #make a dataframe with the coviariates and the bulk expression of a gene
    # Make model that tries to fit bulk from the coviariates
    model = lm(currGene ~ ., merged_df)

    return(model)

}


lo_linear_models = lapply(lo_genes, fit_model_for_one_gene, brain_df, estimations)

exprmatRes <- lo_linear_models  %>% sapply(function(currLm) { currLm$residuals } ) %>% t() # extract residuals

row.names(exprmatRes) = lo_genes

write.csv(exprmatRes, file = 'brain_residuals.csv')
