#!/usr/bin/env python3

import sys
import pandas as pd

def main():
    lo_expression_matrix_paths = sys.argv[1:]
    
    # Extract relevant metadata
    lo_tissues = [path.split("/")[-1].split("_")[0] for path in lo_expression_matrix_paths]
    lo_ct_pro =  [path.split("/")[-1].split("_")[2] for path in lo_expression_matrix_paths]
    lo_type = [path.split("/")[-1].split("_")[8] for path in lo_expression_matrix_paths]
    # load dfs
    lo_expression_matrices = [pd.read_csv(path, index_col=0) for path in lo_expression_matrix_paths]
    
    # Add metadata to dfs
    for i,df in enumerate(lo_expression_matrices):
        tissue = lo_tissues[i]
        ct = lo_ct_pro[i]
        type_m = lo_type[i]
        df['tissue'] = tissue
        df['reg_ct'] = ct
        df['type'] = type_m
    
    # concat
    df = pd.concat(lo_expression_matrices)
    
    df.to_csv("melted_EGAD.csv")
    
if __name__=="__main__":
    main()