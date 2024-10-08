#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import sys



def symboltoENSG(adata:ad.AnnData, pc:pd.DataFrame):
    """filteres for PC genes by subseting adata by only hgnc_symbols that are PC

    Args:
        adata (ad.AnnData): a Bulk dataframe
        pc (pd.DataFrame): the pc master list

    Returns:
        adata: adata object with only PC genes in it. Varnames are hgnc_symbols
    """
    
    subset_adata = adata[:, adata.var_names.isin(pc.loc[:,"hgnc_symbol"])]
    
    return subset_adata
    

def main():
    pcPath = sys.argv[1]
    h5adPath = sys.argv[2]

    
    adata = ad.read_h5ad(h5adPath)    
    pc = pd.read_csv(pcPath, sep = ",")
    
    subsetadata = symboltoENSG(adata=adata, pc = pc)
    print(f'Subsetted Data shape {subsetadata}')
    
    subsetadata.write_h5ad('GtexMergedPC.h5ad')


if __name__ == "__main__":
    main()