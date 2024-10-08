#!/usr/bin/env python3

####
# Purpose:
# Turn all SC tsvs from the human protein atlas into .h5ad objects 
####


import sys
import pandas as pd
import anndata as ad
import os
import numpy as np

def unpackDir(curDir:str) -> tuple:
    files = sorted(os.listdir(curDir))
    
    cell_data,  organism, read_count = files
    
    # prepend directory
    cell_data = f'{curDir}/{cell_data}'
    organism = f'{curDir}/{organism}'
    read_count = f'{curDir}/{read_count}'
    return cell_data,  organism, read_count
 

def main():
    
    cell_data,  organism_part, read_count = unpackDir(sys.argv[1])
    
    organism_file_name = sys.argv[2]
    
    cluster_meta = sys.argv[3]
    
    genes_path = sys.argv[4]
    
    # load Dataframes
    sc_df = pd.read_csv(read_count, sep = "\t", index_col=0)
    
    sc_meta_df = pd.read_csv(cell_data, sep = "\t", index_col=0)
    sc_meta_df.index = sc_meta_df.index.astype(str)
    
    organism_part_df = pd.read_csv(organism_part, header = None)
    organism_part_name = organism_part_df.iloc[0,0] # This is the name of the organism part in the cluster_meta_Df
    cluster_meta_df = pd.read_csv(cluster_meta, sep= "\t")
    
    genes_df = pd.read_csv(genes_path, sep = ",")
    
    # initiate anndata object
    adata = ad.AnnData(sc_df.T, dtype=np.int64)    
    adata.obs.index = adata.obs.index.astype(str)
    
    # Add Gene information
    # add gene data. First drop ensembl_gene_id duplicate
    genes_df = genes_df.drop_duplicates('ensembl_gene_id')
    # merge the gene data to the adata.vars
    adata.var = adata.var.merge(genes_df, left_index = True, right_on = "ensembl_gene_id").set_index('ensembl_gene_id')
        
    # Add Cells information
    adata.obs = adata.obs.merge(sc_meta_df, left_index=True, right_index=True)    
    
    # filter the cluster metadata for only data of specific organism part
    filtered_cluster_metadata = filterClusterMeta(cluster_meta=cluster_meta_df, species = organism_part_name)
    
    # Make sure filtering occured
    if filtered_cluster_metadata.empty:
        raise ValueError(f"Did not filter for organism part properly for {organism_part_name}")
    
    # Create a new column with the cluster id
    filtered_cluster_metadata = extractClusterName(filtered_cluster_metadata=filtered_cluster_metadata)
    
    # Make cluster column str instead of int for the adata
    adata.obs.cluster = adata.obs.cluster.astype('str')
    
    # Merge the dfs
    adata.obs = adata.obs.merge(filtered_cluster_metadata, on = 'cluster', how = 'left')

    # Add gene names to the adata.var
    
    # Write anndata object
    adata.write_h5ad(f'{organism_file_name }_sc_with_metadata.h5ad')


    
    
def filterClusterMeta(cluster_meta:pd.DataFrame, species:str) -> pd.DataFrame:
    """Filter the cluster metadata for a specific species

    Args:
        cluster_meta (pd.DataFrame): _description_
        species (str): _description_

    Returns:
        pd.DataFrame: _description_
    """
    filtered = cluster_meta[cluster_meta.Tissue == species]
    
    return filtered


def extractClusterName(filtered_cluster_metadata: pd.DataFrame) -> pd.DataFrame:
    """Extract the cluster id as a new series column

    Args:
        filtered_cluster_metadata (pd.DataFrame): _description_

    Returns:
        pd.DataFrame: _description_
    """
    cluster = filtered_cluster_metadata.Cluster.apply(lambda x: x.split("-")[1] if "-" in x else "")
    
    filtered_cluster_metadata['cluster'] = cluster
    return filtered_cluster_metadata  



if __name__ == "__main__": 
    main()