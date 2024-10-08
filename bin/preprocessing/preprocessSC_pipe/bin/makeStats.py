#!/usr/bin/env python3


import sys
import scanpy as sc
import pandas as pd

def main():
    adata_path = sys.argv[1]
    adata_name = sys.argv[2]
    
    adata = sc.read_h5ad(adata_path)
    
    # Log the data
    sc.pp.log1p(adata)

    # Get highly variable gene stats
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # PLot highly variable genes
    sc.pl.highly_variable_genes(adata, save = f"_{adata_name}.png")

    # Set the raw to the log cpm 
    adata.raw = adata

    # filter for highly variable genes
    adata = adata[:, adata.var.highly_variable]

    # Control for the effects of library size and mitochondrial percent
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

    # Scale to unit variance for each gene (devide by STdev)
    sc.pp.scale(adata, max_value=10)
    
    #Perform PCA on thehighlr variable genes
    sc.tl.pca(adata, svd_solver='arpack')

    # Make a screeplot 
    sc.pl.pca_variance_ratio(adata, log=True, save = f"_{adata_name}.png")

    # Calc neighbors
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # Calc umap embeddings
    sc.tl.umap(adata) 

    # Perform leiden clustering
    sc.tl.leiden(adata)
    
    # MAke the umap
    sc.pl.umap(adata, color=['leiden'], save = f"_{adata_name}.png")
    
    # Save the processed h5ad
    print(f"~~~~~~~~~~~~~~~~~~~{adata_name}_stats.h5ad")
    adata.write_h5ad(f"{adata_name}_stats.h5ad")
    




if __name__ == "__main__":
    main()