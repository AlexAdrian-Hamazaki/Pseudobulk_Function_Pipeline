#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import numpy as np
import sys



def main():
    bulk_path = sys.argv[1]
    bulk_meta = sys.argv[2]
    df_name = sys.argv[3]
    
    # Load Datasets
    bulk_df = pd.read_csv(bulk_path,
                   skiprows=2,
                   sep = "\t")
    bulk_meta = pd.read_csv(bulk_meta, sep = "\t")

    # Make adata object
    adata = ad.AnnData(bulk_df.set_index("Description").iloc[:, 1:].T, dtype=float)
    
    adata.var_names_make_unique()


    # Set variable information
    adata.var = bulk_df.loc[:,["Name", "Description"]]
    adata.var_names = bulk_df.loc[:,"Description"]

    # Set obs informtaion
    adata.obs_names = bulk_df.iloc[:,2:].columns
    # join metadata onto adata.obs
    adata.obs = adata.obs.merge(bulk_meta, left_index = True, right_on = "SAMPID").set_index("SAMPID")  
      
    #This column has NA values thats messing up saving
    adata.obs = adata.obs.drop(columns = "SMGTC")


    # save
    adata.write(f'{df_name}.h5ad.gz', compression='gzip')
    
if __name__ == "__main__":
    main()