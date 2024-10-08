#!/usr/bin/env python3


import sys
import scanpy as sc
import pandas as pd

def main():
    adata_path = sys.argv[1]
    adata_name = sys.argv[2]
    
    adata = sc.read_h5ad(adata_path)
    


    # PLot highly variable genes
    sc.pl.highly_variable_genes(adata, save = f"_{adata_name}.png")
 
    # Make a screeplot 
    sc.pl.pca_variance_ratio(adata, log=True, save = f"_{adata_name}.png")

    
    # MAke the umap
    sc.pl.umap(adata, color=['Cell type', 'Tissue'], save = f"_{adata_name}.png")
    
    




if __name__ == "__main__":
    main()