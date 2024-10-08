#!/usr/bin/env python3

import sys
import anndata as ad
import pandas as pd


def main():
    tissue_path = sys.argv[1]
    
    adata = ad.read_h5ad(tissue_path)
    
    lo_tissues = ['Brain']
    
    for tissue in lo_tissues:
        adata_tissue = adata[adata.obs.SMTS==tissue]
        adata_tissue.to_df().to_csv(f"{tissue}_gtex.csv")
    
if __name__ == "__main__":
    main()