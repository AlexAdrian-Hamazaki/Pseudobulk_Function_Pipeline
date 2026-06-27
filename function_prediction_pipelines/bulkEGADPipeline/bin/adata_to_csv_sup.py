#!/usr/bin/env python3

import anndata as ad
import sys

def main():
    path_to_h5ad = sys.argv[1]
    
    # load adata
    adata = ad.read_h5ad(path_to_h5ad)

    # Group by 'pbs' column and sample 100 from each group
    sampled_indices = (
        adata.obs.groupby("SMTSD")
        .apply(lambda x: x.sample(n=100, replace=True, random_state=42))  # set replace=True if some groups <100
        .index.get_level_values(1)
    )
    # Subset the AnnData object
    adata = adata[sampled_indices].copy()

    # Save
    adata.to_df().to_csv('gtex.csv.gz')

if __name__ == "__main__":
	main()