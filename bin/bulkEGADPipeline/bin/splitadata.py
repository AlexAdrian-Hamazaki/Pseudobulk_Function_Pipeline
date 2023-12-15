#!/usr/bin/env python3
import os
import sys
import pandas as pd
import anndata as ad


def main():
    merged_adata_path = sys.argv[1]
    organism_part_str = sys.argv[2]
    
    # Load adata
    adata = ad.read_h5ad(merged_adata_path)
    
    # Split adata into organism_part
    adata_organism_part = adata[adata.obs.SMTS == organism_part_str]
    
    # save split
    df = adata_organism_part.to_df()
    df.to_csv(f"{organism_part_str}_split.csv.gz")
    


if __name__ == '__main__':
    main()
    
    
