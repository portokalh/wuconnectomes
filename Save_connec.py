
# coding: utf-8

# In[7]:


import numpy as np
from dipy.io.streamline import load_trk
from dipy.tracking import utils
from dipy.tracking.streamline import Streamlines
from dipy.io.image import load_nifti
import copy
import matplotlib.pyplot as plt


# In[2]:


#set path
mypath = '/Users/wenlin_wu/Downloads/'
outpath = '/Users/wenlin_wu/Desktop/results/'
labels_path = '/Users/wenlin_wu/Desktop/Summer Research/Summer Research Data/N54900/fa_labels_warp_N54900_RAS.nii.gz'
runno = 'N54900'


# In[3]:


#load label information
labels, affine_labels = load_nifti(labels_path)

#transform labels
labels_ = copy.copy(labels)
nonz = np.nonzero(labels)
for i in range(len(nonz[0])):
    if labels_[nonz[0][i], nonz[1][i], nonz[2][i]] >= 1000:
        labels_[nonz[0][i], nonz[1][i], nonz[2][i]] -= 1000
        labels_[nonz[0][i], nonz[1][i], nonz[2][i]] += 166
print runno + ' label transformation finished'


# In[4]:


#load streamlines
streamlines, header = load_trk(mypath+runno+'_bmCSA_detr_small.trk')


#remove the small streamlines
streamlines_cut = lambda: (sl for sl in streamlines if len(sl)>1)
streamlines = Streamlines(streamlines_cut())

print runno + ' streamlines paparation finished'


# In[5]:


#build the connectivity matrix
M= utils.connectivity_matrix(streamlines, labels_, affine=affine_labels,
                                        return_mapping=False,
                                        mapping_as_streamlines=False)
#fill dialgonal with 0
np.fill_diagonal(M,0)

print runno+ ' connectivity matrix building finished'


# In[20]:


M=M[1:,1:]


# In[22]:


plt.imshow(M)


# In[11]:


np.diag(M)


# In[19]:


np.savetxt('/Users/wenlin_wu/Desktop/t1.csv', M, fmt='%.0e', delimiter=',')


# In[23]:


import scipy.io as scio


# In[ ]:


scio.loadmat('')

