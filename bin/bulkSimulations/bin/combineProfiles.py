#!/usr/bin/env python3


import sys
import scanpy as sc
import sys
import pandas as pd


def main():

    lo_adata_paths = sys.argv[1:]
    print(lo_adata_paths)
    # Get the name of the organism part

    
    # Get the adata names. We need this because it contains information about the simulation number
    lo_adata_names = [path.split("/")[-1] for path in lo_adata_paths]
    
    # load all of the anndata paths
    loadata = [sc.read_h5ad(adata_path) for adata_path in lo_adata_paths]
    
    # Get the sum of each genes expression
    loGeneSums = [adata.X.sum(axis = 0) for adata in loadata]
    
    # Create a series where the values are the gene expression
    # the index is the gene variable names
    # and the "name" includes information about the simulation number
    loGeneSumSeries = []
    for i,geneSum in enumerate(loGeneSums):
        ser = pd.Series(geneSum, 
                        index = loadata[i].var_names,
                        name = lo_adata_names[i])
        loGeneSumSeries.append(ser)
                       
                        
    # save to file
    for i, sum_series in enumerate(loGeneSumSeries):
        sum_series.to_csv(f'{lo_adata_names[i]}.csv', index = True)
    
    

if __name__ == "__main__":
    main()