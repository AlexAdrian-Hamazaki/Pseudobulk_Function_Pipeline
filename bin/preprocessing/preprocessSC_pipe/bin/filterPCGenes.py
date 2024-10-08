#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import sys




def main():
    h5ad_path = sys.argv[1]
    name = sys.argv[2]
    pc_map_path = sys.argv[3]
    
    # Load anndata
    adata = sc.read_h5ad(h5ad_path)
    pc_map = pd.read_csv(pc_map_path)
    
    # Create a boolean mask indicating whether a gene is in the gene_symbol_list
    gene_mask = [gene in pc_map.ensembl_gene_id.to_list() for gene in adata.var_names]
    
    # Subset the AnnData object using the gene mask
    adata = adata[:, gene_mask]
    
    adata.write_h5ad(f"{name}_pc.h5ad")




if __name__ == "__main__":
    main()