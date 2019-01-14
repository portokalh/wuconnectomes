
# coding: utf-8

# In[1]:


import multiprocessing
from time import time
import numpy as np
from dipy.tracking import utils
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import (ConstrainedSphericalDeconvModel,
                                   auto_response)
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.streamline import Streamlines
from dipy.tracking.local import LocalTracking, BinaryTissueClassifier
from dipy.reconst import peaks
from dipy.io.streamline import save_trk
import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes


# In[2]:


#Animal List
l=['N54717','N54718','N54719','N54720','N54722','N54759','N54760',
   'N54761','N54762','N54763','N54764','N54765','N54766','N54770','N54771',
   'N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806',
   'N54807','N54818','N54824','N54825','N54826','N54837','N54838','N54843',
   'N54844','N54856','N54857','N54858','N54859','N54860','N54861','N54873',
   'N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892',
   'N54893','N54897','N54898','N54899','N54900','N54915','N54916','N54917']


# In[3]:


def save_trk_heavy_duty(fname, streamlines, affine, vox_size=None, shape=None, header=None):
    """ Saves tractogram files (*.trk)

    Parameters
    ----------
    fname : str
        output trk filename
    streamlines : list of 2D arrays, generator or ArraySequence
        Each 2D array represents a sequence of 3D points (points, 3).
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline points.
    vox_size : array_like (3,), optional
        The sizes of the voxels in the reference image (default: None)
    shape : array, shape (dim,), optional
        The shape of the reference image (default: None)
    header : dict, optional
        Metadata associated to the tractogram file(*.trk). (default: None)
    """
    if vox_size is not None and shape is not None:
        if not isinstance(header, dict):
            header = {}
        header[Field.VOXEL_TO_RASMM] = affine.copy()
        header[Field.VOXEL_SIZES] = vox_size
        header[Field.DIMENSIONS] = shape
        header[Field.VOXEL_ORDER] = "".join(aff2axcodes(affine))

    tractogram = nib.streamlines.LazyTractogram(streamlines)
    tractogram.affine_to_rasmm = affine
    trk_file = nib.streamlines.TrkFile(tractogram, header=header)
    nib.streamlines.save(trk_file, fname)


# In[4]:


# please set the parameter here

#set path
mypath = '/Users/alex/code/Wenlin/data/wenlin_data/'
outpath = '/Users/alex/code/Wenlin/preprocessing/results/'


# In[ ]:


tall = time()
for j in range(55):
    print(j+1)
    runno = l[j]
    print(runno)

    fdwi = mypath+'4Dnii/'+runno+'_nii4D_RAS.nii.gz'    
    ffalabels = mypath+'labels/'+'fa_labels_warp_'+runno+'_RAS.nii.gz'   
    fbvals = mypath+'4Dnii/'+runno+'_RAS_ecc_bvals.txt'
    fbvecs = mypath+'4Dnii/'+runno+'_RAS_ecc_bvecs.txt'

    labels, affine_labels = load_nifti(ffalabels)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    # Correct flipping issue
    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[: ,2]]
    gtab = gradient_table(bvals, bvecs)
    data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm
    # Compute odfs in Brain Mask
    t2 = time()
    radius = 20
    response, ratio, nvl = auto_response(gtab, data, roi_radius=radius, return_number_of_voxels=True)
    print('We use the roi_radius={},\nand the response is {},\nthe ratio is {},\nusing {} of voxels'
    .format(radius, response, ratio, nvl))
    csd_model = ConstrainedSphericalDeconvModel(gtab, response, sh_order=6)
    csd_peaks = peaks_from_model(csd_model, data, sphere=peaks.default_sphere,
                                 relative_peak_threshold=0.5, min_separation_angle=60,mask=mask, 
                                 return_odf=True,return_sh=False, parallel=True,  
                                 nbr_processes=multiprocessing.cpu_count())

    duration2 = time() - t2
    print(runno+' CSD duration %.3f' % (duration2,))
    #mask and classifier
    seeds = utils.seeds_from_mask(mask, density=1,affine=np.eye(4))  
    classifier = BinaryTissueClassifier(bm)
    #local tracking
    streamlines_generator = LocalTracking(csd_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=.5)
    # save a smaller part by only keeping one in 10 streamlines
    t3 = time()
    sg_small = lambda: (s for i, s in enumerate(streamlines_generator) if i % 10 == 0)
    save_trk_heavy_duty(outpath+runno+"_bmCSD_detr_small60odf.trk", streamlines=sg_small,
                        affine=affine,
                        shape=mask.shape, vox_size=vox_size)
    duration3 = time() - t3
    print(runno+' Tracking duration %.3f' % (duration3,))

duration_all = time() - tall
print('All animals CSD det tracking finished, running time is {}'.format(duration_all))


# In[48]:





# In[53]:





# In[56]:





# In[ ]:


# save everything - will generate a 20+ GBytes of data - hard to manipulate
#     sg = lambda: (s for s in streamlines_generator)
#     save_trk_heavy_duty(outpath+runno+"_bmCSD_det.trk", streamlines=sg, affine=affine,
#                        shape=mask.shape,
#                        vox_size=vox_size)

