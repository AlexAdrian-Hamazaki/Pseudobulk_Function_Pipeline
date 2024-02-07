#!/usr/bin/env python

import sys
import pandas as pd

def main():
    go_annot_path = sys.argv[1]
    go_annot_name = sys.argv[2]
    pc_map_path = sys.argv[3]

    # Load go annotations
    go_annot = pd.read_csv(go_annot_path, sep = ",", header = 0, index_col=None)
    
    # load PC map
    pc_map = pd.read_csv(pc_map_path, sep = ",", header = 0, index_col= None)
    
    # Merge the GO annotations onto
    merged = go_annot.merge(pc_map, left_on = "DB_Object_Symbol", right_on = "hgnc_symbol")
    
    # save merged
    print(f"{go_annot_name}_withGeneData.csv")
    merged.to_csv(f"{go_annot_name}_withGeneData.csv.gz", index = False, compression='gzip')
    
if __name__== "__main__":
    main()