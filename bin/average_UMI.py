#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
import csv
import os
import sys
import anndata
import re


# In[189]:


def get_avg_gene_count_for_cell_type(adata, cell_type_column:str, layer_name:str):
    """
    @PARAM: adata is a h5df file

    Output: A file containing the average gene counts for each gene for each cell type in adata
    
    cell_type_column is the str name of the column in adata.obs indicating cell type differences
    """

    try:
        #Check if depth_normed is created
        assert adata.layers['depth_normalized'].shape[0] > 0
    except:
        return (print('DEPTH NORM DOES NOT EXIST'))

        
    avg_UMIs = []
    cell_types = adata.obs[cell_type_column].unique()

    for cell_type in cell_types:



        #Splices H5data to contain only cells of a given cell type
        h5df_grouped_cells = adata[adata.obs[cell_type_column] == cell_type]

        #print('Sending h5df.X to dense format from sparse')
        #denseX = h5df_grouped_cells.X.todense()


        UMI_mean = np.mean(h5df_grouped_cells.layers["depth_normalized"].todense(), axis = 0)
        UMI_mean = np.asarray(UMI_mean)[0]
        #avg = np.squeeze(np.asarray(avg))
        
        avg_UMIs.append(UMI_mean)
    return avg_UMIs


# In[139]:


def test_average_function(adata):
    cell0 =adata.layers['depth_normalized'][0].toarray()[0]
    cell1 =adata.layers['depth_normalized'][1].toarray()[0]
    cell2 =adata.layers['depth_normalized'][2].toarray()[0]
    means = np.mean([cell0, cell1, cell2], axis = 0)
    means = means[means!=0]
    
    matrix = np.mean(adata[0:3].layers["depth_normalized"].todense(), axis = 0)
    array = np.asarray(matrix)[0]
    array = array[array!=0]

    return all(array == means)






# In[173]:


#path = "data/subsampled_normalized/TS_Skin.h5ad"
path = sys.argv[1]
cell_type_column = sys.argv[2]
adata = anndata.read_h5ad(path)


# In[174]:


cell_types = adata.obs[cell_type_column].unique()
gene_symbols = np.array(adata.var['gene_symbol'])

# In[192]:


avg_gene_count = get_avg_gene_count_for_cell_type(adata = adata, cell_type_column =cell_type_column, layer_name = 'depth_normalized')


# In[193]:

# Identify File name
split = re.split('/', path)

comp = re.compile(".*\.h5ad")

file = list(filter(comp.match, split))[0]
print(f'file {file}')

organism_part = re.sub(".h5ad", "", file)

pd_avg_gene_count = pd.DataFrame(avg_gene_count, columns=gene_symbols, index=[organism_part+'_'+ct for ct in cell_types])



# In[ ]:


assert test_average_function(adata)


# In[195]:



outpath = 'data/average_counts'
if not os.path.exists(outpath):
    print('created dir')
    os.makedirs(outpath)

print( f"Wrote Avg Counts for {organism_part} and its cell types")

pd_avg_gene_count.to_csv(f"{outpath}/{organism_part}.csv",
                         header = True,
                         index = True,
                         index_label = 'Cell_Type')

