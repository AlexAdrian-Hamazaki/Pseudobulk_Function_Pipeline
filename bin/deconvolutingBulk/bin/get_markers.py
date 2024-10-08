#!/usr/bin/env python3

import sys
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def main():
    path_to_h5ad = sys.argv[1]
    tissue = sys.argv[2]
    # Load adata object
    adata =sc.read_h5ad(path_to_h5ad)
    
    # Log the data to get logcpm
    sc.pp.log1p(adata)
    
    # get the highly variable genes for later anlysis
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Subset adata for only highly variable genes
    # adata = adata[:,adata.var.highly_variable]
    
    # Perform 1 vs  unpaired all wilcoxon tests for each CT and gene
    sc.tl.rank_genes_groups(adata, "Cell type", method='wilcoxon', use_raw=False)
    
    # Fix the index names to gene symbols
    adata.var.reset_index(inplace=True)
    adata.var_names = adata.var.loc[:,'hgnc_symbol'].values.tolist()
    
    # Perform PCA and plot first 2 PCs
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='Cell type', save = f"{tissue}_pca.png")
    # Plot Scree
    sc.pl.pca_variance_ratio(adata, log=True, save = f"{tissue}_scree.png" )
    # Plot loadings of first 3 PCs
    sc.pl.pca_loadings(adata, save = f"{tissue}_loadings.png" )
    
    # Calculate UMAP manifold and plot
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='Cell type', save =f"{tissue}_umap.png" )
    
    # Calculate Cell Type marker via wilcoxon tets
    sc.tl.rank_genes_groups(adata, 'Cell type', method='wilcoxon')
    # Plot top differentially expressed genes
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = f"{tissue}_exp_marker_genes.png")
    
    # Extract markers
    df_ranked = sc.get.rank_genes_groups_df(adata = adata, group = None, log2fc_min=3)
    # Add ensemble gene ids back on
    df_ranked = df_ranked.merge(adata.var.loc[:,"ensembl_gene_id"], left_on = 'names', right_index=True)

    # Graph log fcs of marker genes
    plt.figure()
    g = sns.FacetGrid(df_ranked, col="group")
    g.map(sns.histplot, "logfoldchanges")
    g.savefig(f"{tissue}_marker_dist.png")

    plt.figure()
    # Graph marker counts per ct
    g2 = sns.barplot(df_ranked.groupby("group").size().reset_index(), x = 'group', y = 0)
    # Rotate x-axis labels by 30 degrees
    g2.set_xticklabels(g2.get_xticklabels(), rotation=30)
    g2.set_title("Num CT Marker Genes per CT")
    plt.savefig(f"{tissue}_marker_counts_per_ct.png")

    # save marker genes
    df_ranked.to_csv(f"{tissue}_markers.csv", index = None)
    
    # Get Random Marker Genes
    randoms = get_random_markers(path_to_h5ad, df_ranked=df_ranked)
    randoms.to_csv(f"{tissue}_random_markers.csv", index = None)

def get_random_markers(path_to_h5ad:str, df_ranked:pd.DataFrame) -> pd.DataFrame:
    """subsample a df_ranked number of genes randomly. And then attribute them to the cell types in df_ranked


    Args:
        path_to_h5ad (str): path to adata single cell data
        df_ranked (pd.DataFrame): df of cell type markers previously gotten

    Returns:
        pd.DataFrame: A dataframe of randomly selected cell type marker genes where each cell type has the same number of random cell type marker genes as it does real ct marker genes
    """
    adata = sc.read_h5ad(path_to_h5ad) # re-load the adata object because the previous was subset by highly variable genes
    sc.pp.filter_genes(adata, min_counts=2)
    sc.pp.filter_genes(adata, min_cells=5)
    # Get number of genes to sample
    n_genes = df_ranked.shape[0]
    # subsample genes without replacement
    df_genes = adata.var.reset_index().sample(n=n_genes).loc[:,["ensembl_gene_id","hgnc_symbol"]]
    # Set the genes to be affiliated with the cell types in df_ranked
    df_genes['group'] = df_ranked.loc[:,'group'].values
    df_genes.columns =  ['ensembl_gene_id', 'names', 'group']
    df_genes = df_genes.loc[:,['group', 'names', 'ensembl_gene_id']] # re-organize
    return df_genes.reset_index(drop=True)
if __name__=="__main__":
    main()