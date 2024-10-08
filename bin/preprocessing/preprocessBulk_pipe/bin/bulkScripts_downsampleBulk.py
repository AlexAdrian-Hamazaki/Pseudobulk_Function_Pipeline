#!/usr/bin/env python3

# Instead of creating an expression profile for bulk. Just downsample

import anndata as ad
import numpy as np
import sys
import scanpy as sc
import pandas as pd
import os


def group_by_op(adata:ad.AnnData):
    """Group an anndata object by "Organism Part"
    #TODO When more data is expanded for this pipeline. We will need to ensure organism part column IDs are the same. And then subset by that colunm name
    

    Args:
        adata (ad.Anndata): Anndata object of bulk data

    Returns:
        adata :  Anndata object grouped by organism part. Note. "SMTS" may not always be the .obs column to group by
    """

    return adata.obs.groupby('SMTS')

def get_group_index(adata:ad.AnnData):
    """Get dictionary that looks like
    {OrganismPart: [Samples Ids that are from that OrganismPart]}

    Args:
        adata (ad.Anndata): bulk andnata

    Returns:
        dict: Dictionary that looks like {OrganismPart: [Samples Ids that are from that OrganismPart]}
    """
    
    grouped = group_by_op(adata)
    


    OPs_indexes = {}

    # You can loop over a grouped Anndata Object.
    # It will return 2 items. 
    # group_name: The name that the data is grouped by
    # group_df: The anndata object
    for group_name, group_df in grouped:

        OPs_indexes.update({group_name: group_df.index})
            
    return OPs_indexes

def subsampleByIndex(adata:ad.AnnData, indexes:list):
    """Subsample an anndata object by alist of observervations

    Args:
        adata (ad.AnnData): _description_
        indexes (list): _description_

    Returns:
        ad.Anndata: a subsampled anndata for one organism part
    """
    return adata[indexes,:]


# %%
def sampleKOPs(adata:ad.AnnData, k = 20):
    """Subsample an andata object. 
    Sample K samples for each organism part

    Args:
        adata (ad.Anndata): Anndata object
        k (int, optional): K to subsampel to. Defaults to 20.

    Returns:
        adata: an anndata object that has only K samples of each organism part
    """
    OPs_indexes = get_group_index(adata)
    
    # list of anndata objects
    loadata = []
    
    for OP_key in OPs_indexes.keys(): # For each Organism Part
        indexes = OPs_indexes[OP_key].tolist() #These are the sample Ids that belong to this organism part
        
        OP_adata = subsampleByIndex(adata, indexes)

        if OP_adata.shape[0] < k: #If there is not enough data to downsample to K, just take the whole organism part data
            print("WARNING: {OP_key} did not have sufficient data to sample K")
            loadata.append(OP_adata)
            continue
            
        
        #index_samples = np.random.choice(np.array(indexes), k, replace = True)
        subsampled_adata = sc.pp.subsample(OP_adata, n_obs=k, copy=True)
        
        loadata.append(subsampled_adata) # Add These sampled adata objects into the list of adata objects
        
    return ad.concat(loadata) # Concats the adata objects in a list into 1 adata object
        
   
   
def main():
    
    args = sys.argv

    adataPath = str(args[1])
    k = int(args[2])

    bootstrap = int(args[3])

    # Reads in anndata object
    adata = ad.read_h5ad(adataPath)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    
    # make directory for csvs to be stored
    os.mkdir('csvs')
    # bootstrap several subsampled
    n = 1
    while n < bootstrap+1:
        print(f'~~~~~~{k}~~~~{n}')
        adata_K = sampleKOPs(adata, k)
        adata_K.write_h5ad(f'K{k}_boot{n}_bulk.h5ad')
        adata_K.to_df().to_csv(f'csvs/K{k}_boot{n}_bulk.csv')
        n+=1
   
     
if __name__ == "__main__":
    main()
    





