#!/usr/bin/env python3

import sys
import scanpy as sc
import pandas as pd
import numpy as np

def main():
    path_to_h5ad = sys.argv[1]
    tissue = sys.argv[2]
    # Load adata object
    adata =sc.read_h5ad(path_to_h5ad)
    
    dic_markers = {}
    for ct in adata.obs.loc[:,'Cell type']:
        dic_markers[f'{ct}'] = np.random.choice(adata.var_names, 100)
    
    lo_dfs = []
    for k in dic_markers.keys():
        df = pd.DataFrame(dic_markers.get(k))
        with_hgnc = add_hgnc(df, adata)
        with_hgnc.columns = [f"{k}_ensg", f"{k}_hgnc"]
        lo_dfs.append(with_hgnc)
    
    concat = pd.concat(lo_dfs, axis = 1)
    # Save 
    concat.head(100).to_csv(f"{tissue}_random_markers.csv")
    
def add_hgnc(df, adata):
    merged = df.merge(adata.var, left_on=0, right_index=True)
    merged = merged.loc[:,[0, 'hgnc_symbol']]
    merged.columns = ['ensg', 'hgnc_symbol']
    return merged

    
if __name__=="__main__":
    main()