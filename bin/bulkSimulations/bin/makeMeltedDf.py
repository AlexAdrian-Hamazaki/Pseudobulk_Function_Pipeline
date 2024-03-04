#!/usr/bin/env python3
import pandas as pd
import sys

def main():
    tissue = sys.argv[1] # This will be a string
    lo_EGAD_paths = sys.argv[2:] # this will be a lsit of EGADS of a variety of variances
    
    # Isolate the variance level of each dataframe
    lo_variances = [get_variance_int(path) for path in lo_EGAD_paths] # checked
    
    # First lets load all of the EGADs
    lo_EGADs = [pd.read_csv(path, index_col=0) for path in lo_EGAD_paths]
    
    # Process the column names of all the AUCs of the lo_EGAds
    df_EGAD = process_auc_colnames(lo_dfs=lo_EGADs, lo_variance_levels=lo_variances) # Now we have a df of all the EGADs

    # add columns of interest
    df_EGAD['tissue'] = tissue
    
    # Melt
    df_melted = df_EGAD.reset_index().melt(id_vars = ['index', 'tissue'],
                                            var_name='variance',
                                            value_name='auc')
    
    df_melted.to_csv(f"{tissue}_melted_df.csv.gz", compression='gzip')    

    
    
def process_auc_colnames(lo_dfs:list, lo_variance_levels:list) -> pd.DataFrame:
    lo_ser_pro = [] # list hold AUC vectors wwhere names are the variance level
    for i, df in enumerate(lo_dfs):
        var_int = lo_variance_levels[i]
        # Renaming the 'auc' column to 'var'
        df_pro = df.rename(columns={'auc': str(var_int)})
        auc_ser = df_pro.loc[:,str(var_int)]
        lo_ser_pro.append(auc_ser)
    concat_df = pd.concat(lo_ser_pro, axis = 1)
    return concat_df 
    

def get_variance_int(path:str) -> str:
    return str(path.split('_')[-2])

if __name__ == "__main__":
    main()