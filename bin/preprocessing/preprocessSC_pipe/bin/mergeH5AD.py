#!/usr/bin/env python3


import sys
import scanpy as sc
import anndata as ad

def main():
    loadata_paths = sys.argv[1:]
    
    
    loadata = [sc.read_h5ad(adata_path) for adata_path in loadata_paths]
    
     # Concatenate the anndata objects along the observations axis
    adata_merge = ad.concat(loadata, axis=0, join = 'outer')
    
    adata_merge.write_h5ad("merged_adata.h5ad")

    

if __name__ == "__main__":
    main()