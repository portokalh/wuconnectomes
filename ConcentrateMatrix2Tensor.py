
# coding: utf-8

# In[1]:


from __future__ import division
import numpy as np
import pandas as pd
import scipy.io as scio


# In[2]:


l=['N54717','N54718','N54719','N54720','N54722','N54759','N54760','N54761','N54762','N54763','N54764','N54765','N54766','N54770','N54771','N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806','N54807','N54818','N54824','N54825','N54826','N54837','N54838','N54843','N54844','N54856','N54857','N54858','N54859','N54860','N54861','N54873','N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892','N54893','N54897','N54898','N54899','N54900','N54915','N54916','N54917']


# In[3]:


mypath = '/Users/wenlin_wu/Downloads/dipyconnectomes55_10p/'
outpath = mypath


# In[4]:


M_all = np.zeros((332,332))
for i in range(55):
    runno = l[i]
    M = pd.read_csv(mypath+runno+'_connectivityCSAbm.csv',delimiter=',',header=None,dtype=np.int)
    #M = (M-M.values.min())/(M.values.max()-M.values.min())
    #print(np.max(np.diag(M)))
    M_all = np.dstack((M_all,M))
M_all = M_all[:,:,1:56]


# In[25]:


scio.savemat(outpath+'connectivity_all332DipyNor.mat', {'connectivity332Nor': M_all})


# In[100]:


np.sum(M_all,1).shape


# In[13]:


l[3]


# In[12]:


l[19]


# In[9]:


np.max(M_all[:,:,45])


# In[8]:


M_all[:,:,51]

