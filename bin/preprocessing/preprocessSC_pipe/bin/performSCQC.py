#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import sys



def main():
    h5ad_path = sys.argv[1]
    h5ad_name = sys.argv[2]
    
    adata = sc.read_h5ad(h5ad_path)
    
    # Get Graphs of the highest expressed genes in the dataset
    sc.pl.highest_expr_genes(adata, n_top=20, save = f'_{h5ad_name}.png')
    
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

    #total_counts: the number of genes expressed in the cell (library size)
    #n_genes_by_counts: Number of genes expressed in each cell
    #pct_counts_mt: the percentage of counts in mitochondrial genes
    
    # Plot the QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = f'_{h5ad_name}.png')
    
    # Plot Scatterplots of QC metrics
    
    # Figure out where the threshold is of sufficiently high mitochondrial percentage expression
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save = f'_{h5ad_name}_mito.png')
    
    # Figure out threshold of where to remove duplicates. These would have a very high number of genes by counts
    # and a very high number of total counts may indicate a duplet
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = f'_{h5ad_name}_duplet.png')
    

    











if __name__ == "__main__": 
    main()