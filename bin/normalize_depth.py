#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import anndata
import os
import sys
import scipy
import re


# # Normalize Single Cell Data By Depth

# In[27]:


#path = '../data/h5ad/subsampled/TS_Skin_Subsampled.h5ad'
path = sys.argv[1]


# In[28]:


adata = anndata.read_h5ad(path)


# In[29]:


######## Global Variables ############

### Gene Names
gene_symbols = adata.var.gene_symbol


### Num Cells


### Cell Types
#Larger compartments
#cell_types = adata.obs.compartment.unique()
cell_types = adata.obs.cell_ontology_class.unique()


### Organism_Part
organism_part = adata.obs.organ_tissue.unique()[0]



# # Normalize for UMI Depthdef sum_UMIs_per_cell(adata):
# 

# ### Step 1. Get Sum of UMIs for Each Cell

# In[30]:


def sum_UMIs_per_cell(adata):
    X_counts_UMIs = [np.sum(adata[cell_index].X) for cell_index in range(adata.shape[0])]
    return np.array(X_counts_UMIs)

X_counts_UMIs = sum_UMIs_per_cell(adata)


# In[31]:



adata.obs['X_counts_UMIs'] = X_counts_UMIs


# In[32]:


assert len(X_counts_UMIs) == adata.shape[0]
assert np.sum(adata[0].X) == X_counts_UMIs[0]
assert np.sum(adata[adata.shape[0]-1].X) == X_counts_UMIs[adata.shape[0]-1]


# ### Step 2. Add the Sum to the adata.obs

# In[33]:

adata.obs['X_counts_UMIs'] = X_counts_UMIs


# ### Step 3. Devide each UMI By its sum count

# In[34]:


def normalize_UMIs(adata):
    """
    Requires adata.obs to contain X_counts_UMIs
    
    Normalize adata.X by deviding it by the indexed value of X_counts_UMIs

    Adds depth_normed layer to adata and saves to disk
    """
    depth_normed_layer = []
    for i, cell in enumerate(adata):
        depth_normed_cell = adata.X[i] / adata.obs.X_counts_UMIs[i]

        depth_normed_layer.append(depth_normed_cell)

    sparse_normed_layer = scipy.sparse.vstack(blocks = depth_normed_layer)
    
    adata.layers['depth_normalized'] = sparse_normed_layer

    return None
depth_normed = normalize_UMIs(adata)


# In[35]:


assert adata.layers['depth_normalized'].shape[0]>0


# In[36]:


assert len(adata.layers['depth_normalized'][0].toarray()[0]) == adata.shape[1]


# #adata.layers['depth_normalized'][0].toarray()[0] [adata.layers['depth_normalized'][0].toarray()[0] !=0] [0]

# ###  Write to Disk

# In[17]:


# Make Path

outpath = "../../../../pipeline42/datasets/TabulaSapiens/normalized"
if os.path.exists(outpath) == False:
    os.mkdir(outpath)


# In[ ]:


# Identify File name
split = re.split('/', path)

comp = re.compile(".*\.h5ad")

file = list(filter(comp.match, split))[0]

print(file)


# In[43]:

save_path = f'{outpath}/{file}'
adata.write(save_path)


# In[44]:


adata.file.close()

