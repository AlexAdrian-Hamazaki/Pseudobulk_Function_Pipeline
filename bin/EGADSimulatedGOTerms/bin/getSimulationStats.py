#!/usr/bin/env python3

import os
import sys
import pandas as pd
from functools import reduce
import re
import json
import plotly.express as px


def main():
    variance_level = sys.argv[1]
    op_name = sys.argv[2]
    
    print(op_name)
    lo_composition_file_paths = sys.argv[3:]
        
    # Load all the pandas dataframes
    lodfs = [pd.read_csv(file_path, index_col=0) for file_path in lo_composition_file_paths]
    
    # Get the num of each simulation
    lo_sim_ids = [path.split("/")[-1] for path in lo_composition_file_paths]
    
    # Concat all the dfs and change column names
    merged_df = pd.concat(lodfs, axis = 1)
    merged_df.columns = lo_sim_ids
    
    # Save this file as the merged proportions for all of our cells
    merged_df.to_csv(f"all_simulation_{variance_level}_{op_name}_proportions.csv")
    
    # GEt the avg cell count proportion and save that as well
    avg_df = merged_df.mean(axis = 1)
    avg_df.to_csv(f"average_simulation_{variance_level}_{op_name}_proportion.csv")
    
if __name__ == "__main__":
    main()