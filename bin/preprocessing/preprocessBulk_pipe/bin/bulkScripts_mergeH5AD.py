#!/usr/bin/env python3

import sys
import anndata as ad
import os




def loadLOA(loaPaths:list):
    """Loads anndata objects for a list of paths

    Args:
        loaPaths (list): List of paths toanndata objects

    Raises:
        FileExistsError: if file is not found
    Returns:
        list: List of loaded anndata objects
    """
    #list of anndata objects
    loaObjects = []
    
    for file_path in loaPaths:
        if os.path.isfile(file_path):
            adata = ad.read_h5ad(file_path)
            loaObjects.append(adata)
        else:
            raise FileExistsError(f"File not found: {file_path}")
    return loaObjects
        
        
        
def main():
    loaPaths = sys.argv[1:]
    
    
    loaObjects = loadLOA(loaPaths=loaPaths)
    
    #Merged
    merged_adata = ad.concat(loaObjects, join = 'outer')

    #Make observation names unique. Hopefully this should yield no changes
    merged_adata.obs_names_make_unique()

    
    merged_adata.write_h5ad('bulkMerged.h5ad')

        
if __name__ == "__main__":
    main()
