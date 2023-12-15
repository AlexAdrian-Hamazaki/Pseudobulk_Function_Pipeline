#!/usr/bin/env python3
import os
import sys
import pandas as pd
import anndata as ad
import argparse
import ast



def main():
    
    parser = argparse.ArgumentParser(description="Make a melted dataframe for each organism part that has EGAD values for each bootstrap")
    
    parser.add_argument('--lo_bootstraps', type=ast.literal_eval, required = True, help='The first list as a Python expression')
    parser.add_argument('--lo_EGAD_paths', type=str, nargs='+', required = True, help='The second list as a Python expression')
    parser.add_argument('--organism_part', type=str, required=True, help='The organism part as a string')
    args = parser.parse_args()
    
    organism_part = args.organism_part # This will be a string
    lo_bootstrap = args.lo_bootstraps # this will be a string 
    lo_EGAD_paths = args.lo_EGAD_paths # this will be a lsit of EGADS of a variety of variances
    

    
    # fix the egad paths cuse nextflow is weirda
    lo_EGAD_paths = lo_EGAD_paths[0]
    lo_EGAD_paths = lo_EGAD_paths.replace(" ", "\",\"")
    lo_EGAD_paths = "[\"" + lo_EGAD_paths + "\"]"
 
 
    lo_EGAD_paths = ast.literal_eval(lo_EGAD_paths)

    # Load each dataframe for this organism part, and then add a column indicating what bootstrap it is
    lo_dfs = [pd.read_csv(path, index_col=0) for path in lo_EGAD_paths]
    
    # add bootstrap id to each bootstrap
    lo_dfs_with_boot = [df.assign(bootstrap=lo_bootstrap[n]) for n, df in enumerate(lo_dfs)]
    
    # melt each dataframe
    lo_melted_dfs = [df.reset_index().melt(id_vars = ['index', 'bootstrap'],value_vars = ['auc']) for df in lo_dfs_with_boot]

    # print(lo_melted_dfs[0].head())
    
    # Concat the dataframes
    merged_df = pd.concat(lo_melted_dfs, axis = 0)

    del merged_df['variable']
    merged_df.rename(columns = {'value':'auc'}, inplace = True)
    merged_df.to_csv(f"{organism_part}_melted_EGADs.csv.gz", compression='gzip')    
    
    

if __name__ == "__main__":
    main()