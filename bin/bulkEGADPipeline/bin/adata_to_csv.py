#!/usr/bin/env python3

import anndata as ad
import sys

def main():
    path_to_h5ad = sys.argv[1]
    
    # load adata
    adata = ad.read_h5ad(path_to_h5ad)
    # Save
    adata.to_df().to_csv('gtex.csv.gz')

if __name__ == "__main__":
	main()