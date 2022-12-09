#!/usr/bin/env python
# coding: utf-8

# In[66]:


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import os


# In[7]:


path = 'data/EAGD/EGAD.csv'


# In[62]:


save_path = 'data/figs'


# In[8]:


EGAD = pd.read_csv(path)


# In[42]:


EGAD


# In[14]:


EGAD_auc = EGAD.loc[:,"auc"]


# In[22]:


EGAD_auc = EGAD_auc[~np.isnan(EGAD_auc)]


# In[45]:


np.absolute(EGAD_auc).max()


# In[57]:


plt.hist(np.log(np.absolute(EGAD_auc+1)), bins = 10, range = (0, 1))
plt.xlabel('log AUC')
plt.ylabel('Count')


# In[67]:

if not os.path.isdir(save_path):
    os.mkdir(save_path)
plt.savefig(f'{save_path}/auc.jpeg')

