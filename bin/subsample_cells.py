#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import anndata
import re
import sys
import os


# ### Purpose
# 
# So Ishka can run things quickly, I'll remove cells and genes. 
# 
# Also, saves the adata.X dense object
# 

# In[2]:




# In[2]:


#path = '../data/h5ad/'
path = sys.argv[1]
save_path = 'data/subsampled'


# In[ ]:


#os.chdir("..")
# files = os.listdir(path)


# In[11]:


# true_files = []
# for file in files:
#     s = re.search('.h5ad', file)
#     if s != None:
#         true_files.append(file)


# In[13]:


def main(path:str, from_cli = True):
    
    if from_cli == True:
        adata = anndata.read_h5ad(path)
        adata = adata[0:int(adata.shape[0]/4), 0:int(adata.shape[1]/5)]

        split = re.split('/', path)

        comp = re.compile(".*\.h5ad")

        file = list(filter(comp.match, split))[0]

        print(file)

        outpath = 'data/subsampled'
        if os.path.exists(outpath) == False:
            os.mkdir(outpath)
            
        final_save = f'{outpath}/{file}'

        adata.write(final_save)
        


# In[5]:


if __name__ == "__main__":
    main(path)

