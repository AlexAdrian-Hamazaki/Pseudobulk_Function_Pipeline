#!/usr/bin/env python3

import sys
import anndata as ad
import pandas as pd
import numpy as np

def main():
    tissue_path = sys.argv[1] 
    # First, load the tissue
    adata = ad.read_h5ad(tissue_path)
    
    # Make a dict where keys are cell type, and values are average gene expression for that cell type
    dictCellTypeProfile = makeDictCellTypeProfile(adata = adata)
    
    # Make a pandas dataframe from the dictionary
    df = pd.DataFrame(dictCellTypeProfile).T
    # update columns
    df.columns = adata.var_names
    
    # Save the dataframe
    df.to_csv(f"all_human_CT_profiles.csv")
    
    # Make control ct profiles where the genes are permuted across the cell types
    # Shuffle around the CT Profile information (random control)
    # control_CT_profile_df = shuffle_CT_Profile_df(df)
    # #save 
    # control_CT_profile_df.to_csv(f"cntrl_{tissue}_cell_type_profiles.csv")

def shuffle_CT_Profile_df(CT_profile_df:pd.DataFrame) -> pd.DataFrame:
    """Randomly shuffle each column in the CT profile DF such that the CT profiles will no longer be comprised of actual true CT profile expression.
    
    THe expression of a gene in each CT profile will ranodmly become the expression of one of our CT profiles
    

    Args:
        CT_profile_df (pd.DataFrame): columns are genes, rows are CT profiles, values are expression

    Returns:
        pd.DataFrame: control_CT_profile_df: Expression of each gene has been shuffled between CT profiles
    """
    control_CT_profile_df = CT_profile_df.apply(np.random.permutation)
    return control_CT_profile_df
    
    
    


    
def makeDictCellTypeProfile(adata: ad.AnnData) -> dict:
    """Make a dictionary where keys are the cell type, and values is a list of gene expression values
    
    The gene expression values is a cell type profile that is generated by taking the average expression of all of the cells in a cell type

    Args:
        adata_grouped (ad.AnnData): anndata object that is grouped by cell type
        adata (ad.AnnData): anndata object that is not goruped

    Returns:
        dict: Keys are cell types, values are a list of gene expression values
    """
    #Group the anndata by cell type
    adata_grouped = adata.obs.groupby(["Cell type"])
       
    dictCellTypeProfile = {}

    for group, lo_cell_ids in adata_grouped.indices.items():
        
        # Retrive all the cells in a cell type
        grouped_cells = adata[lo_cell_ids,:]
        
        
        # fill NAs
        grouped_cells_df = grouped_cells.to_df().fillna(0)

        
        # Summarize the counts for that cell type
        cell_type_profile = grouped_cells_df.mean(axis = 0)
        
        
        # Add to the dict. key will be the cell type(group) and the value would be a list of summarized genes
        dictCellTypeProfile[group] = cell_type_profile
        
    return dictCellTypeProfile

if __name__ == "__main__":
    main()