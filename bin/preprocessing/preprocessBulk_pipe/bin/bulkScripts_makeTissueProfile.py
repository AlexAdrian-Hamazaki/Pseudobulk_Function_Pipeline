#!/usr/bin/env python3

# Instead of creating an expression profile for bulk. Just downsample

import anndata as ad
import numpy as np
import sys
import scanpy as sc
import pandas as pd


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
def collapseOPs(adata:ad.AnnData):
    """Collapse all values of a given organism part into 1 average value per gene.
    

    Args:
        adata (ad.Anndata): Anndata object

    Returns:
        adata: An anndata object that has a 1 row per organism part.
    """
    OPs_indexes = get_group_index(adata)
    
    # list of anndata objects. Each object will be for 1 organism part
    loadata = []
    
    for OP_key in OPs_indexes.keys(): # For each Organism Part
        indexes = OPs_indexes[OP_key].tolist() #These are the sample Ids that belong to this organism part
        
        # OP_adata is all rows for 1 organism part
        OP_adata = subsampleByIndex(adata, indexes)
        
        OP_adata.var_names_make_unique()
        
        
        loadata.append(OP_adata) # Add These sampled adata objects into the list of adata objects

    print(loadata[0])
    print(loadata[0].obs_names)
    
    concat_adata = ad.concat(loadata, axis = 0)

    return concat_adata # Concats the adata objects in a list into 1 adata object
        
   
   
def main():
    
    args = sys.argv

    adataPath = str(args[1])

    # Reads in anndata object
    adata = ad.read_h5ad(adataPath)
    adata.obs_names_make_unique()
    

    adata = collapseOPs(adata)
    adata.write_h5ad(f'organismPartProfile.h5ad')
   
     
if __name__ == "__main__":
    main()
    





