#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import sys



def main():
    h5ad_path = sys.argv[1]
    h5ad_name = sys.argv[2]
    curated_cutoffs_path = sys.argv[3]
    
    adata = sc.read_h5ad(h5ad_path)
    curated_cutoffs = pd.read_csv(curated_cutoffs_path)
    
    print(adata)
    
    # Perform obvious filtering. Remove cells with less than 200 genes
    # Remove genes that are expressed in less than 4 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Add a column in adata.var indicating if the gene is mitochondrial or not
    adata.var['mt'] = adata.var.hgnc_symbol.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    
    # Fill the missing values with False
    adata.var.mt = adata.var.mt.fillna(False)

    # Calculate Qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    
    # Identify what organism part we are at
    organism_part = adata.obs.Tissue[0]
    
    # get the cutoffs for this organism part
    organism_part_cutoffs = curated_cutoffs[curated_cutoffs.organism_part == organism_part]
    
    # Get the various cutoffs 
    pct_mt_cutoff = organism_part_cutoffs.pct_mt_cutoff.iloc[0]
    n_genes_cutoff = organism_part_cutoffs.n_genes_cutoff.iloc[0]
    total_counts_cutoff = organism_part_cutoffs.total_counts_cutoff.iloc[0]
    
    # Slice the df based on these curated cutoffs
    adata = adata[adata.obs.pct_counts_mt <= pct_mt_cutoff,:] # Remove high mitochondrial cells
    adata = adata[adata.obs.n_genes_by_counts <= n_genes_cutoff,:] # Remove cells with too many genes
    adata = adata[adata.obs.total_counts <= total_counts_cutoff,:] # Remove cells with too many counts

    # write with NO CPM
    adata.write_h5ad(f"{h5ad_name}.h5ad")

    # Calculate CPM
    sc.pp.normalize_total(adata, target_sum=1e6, key_added="norm_factor",
                          inplace=True)
        
    adata.write_h5ad(f"{h5ad_name}_cpm.h5ad")

    






















if __name__ == "__main__": 
    main()