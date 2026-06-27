#!/usr/bin/env python3

import pandas as pd
import sys



def main():
   
    lo_EGAD_paths = sys.argv[1:] # this will be a lsit of EGADS of a variety of variances
    
    lo_melted_egads = []
    n = 0
    for path in lo_EGAD_paths:
        if n == 0:
            df = pd.read_csv(path, index_col=0)
            lo_melted_egads.append(df)
            n+=1
        else:
            df = pd.read_csv(path, index_col=0)
            # Remove column names by setting them to empty strings
            # df.columns = [""] * len(df.columns)
            lo_melted_egads.append(df)

    all_egads = pd.concat(lo_melted_egads, axis = 0)
    all_egads.to_csv("master_melted_df.csv.gz", compression='gzip')

if __name__ == "__main__":
	main()