#!/usr/bin/env python3
#$ -l h_vmem=50000M,vf=50000M
#$ -M ${USER}@duke.edu 
#$ -m ea 
#$ -o /mnt/munin6/Badea/Lab/mouse/sinha_sbatch/slurm-$JOB_ID.out 
#$ -e /mnt/munin6/Badea/Lab/mouse/sinha_sbatch/slurm-$JOB_ID.out
#$ -N ${1}_LPCA_denoising

import os.path
from os import path
import sys

#runno=sys.argv[1] # switching to more generic "id"
id=sys.argv[1]
fdwi=sys.argv[2]
#bval_folder=sys.argv[3] # Need to have gtab handling that 1: takes in a single value, & 2: supports DSI Studio btables
# However, for now we'll stick with the bval/bvec pair to keep things moving forward.
bval_or_bvec_or_btable=sys.argv[3]
outpath=sys.argv[4]
#find other pointers@
#https://github.com/nipy/nipype/blob/fbf2c35f533b7805ca93c742006472e0809d8d03/nipype/workflows/dmri/mrtrix/diffusion.py
#to do: coreg/eddy correction/bias field (if not part of denoising already)
#https://github.com/nipy/dipy/blob/349e6b181ac89f333f07146751a2295b732b5c04/scratch/very_scratch/registration_example.py

from time import time
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import CsaOdfModel
import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from dipy.workflows.denoise import NLMeansFlow
from dipy.denoise.denspeed import nlmeans_3d
from dipy.denoise.localpca import localpca
from dipy.denoise.pca_noise_estimate import pca_noise_estimate
from dipy.denoise.non_local_means import non_local_means
from dipy.denoise.adaptive_soft_matching import adaptive_soft_matching
from dipy.denoise.nlmeans import nlmeans
from dipy.denoise.noise_estimate import estimate_sigma
from dipy.segment.mask import median_otsu
import matplotlib.pyplot as plt
from time import time
from dipy.segment.mask import median_otsu
from dipy.tracking.utils import random_seeds_from_mask
from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.streamline import Streamlines

from dipy.io.streamline import save_trk


np.seterr(divide='ignore', invalid='ignore')

print(id)

t1 = time()

#fbval=bval_folder + runno + '_fsl_bvals.txt' 
#fbvec=bval_folder + runno + '_fsl_bvecs.txt'

# The following shenanigans should allow us to specify EITHER the bval OR the bvec file.
# It is assuming that the both have the same prefix, with the exeception of ending in bvecs.txt or bvals.txt

# 
fbval=bval_or_bvec_or_btable.replace("bvecs.txt","bvals.txt")
fbvec=fbval.replace("bvals.txt","bvecs.txt")

bvals, bvecs = read_bvals_bvecs(fbval, fbvec)

    
print(f'b values {bvals}')
print(f'b vecs {bvecs}')
    

gtab = gradient_table(bvals, bvecs)

# Currently need to run over pre-masked data.
no_masking=1
if no_masking:
   # masked_path=outpath + id + '_nii4D_masked.nii.gz'
    masked_path=fdwi
    if path.exists(masked_path):
        maskdata, affine, vox_size = load_nifti(masked_path, return_voxsize=True)
else:
    data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)
    idx = bvals > 10
    print(f'non zero b values {bvals[idx]}')
    #mask1 = data[..., 0] > 0
    #maskdata, mask = median_otsu(data, 3, 1, True, dilate=2)
    maskdata, mask = median_otsu(data, idx,  3, 1, False, dilate=2)
    
    #mask_img = nib.Nifti1Image(maskdata.astype(np.int16), affine)
    #nib.save(mask_img, 'tensor_fa.nii.gz')
    # masked_path=outpath+runno+'_nii4D_masked.nii.gz'
    save_nifti(masked_path, maskdata, affine)
    mask_path=outpath + id +'_mask.nii.gz'
    save_nifti(mask_path, mask.astype(np.ubyte), affine)
        
data=maskdata

from dipy.reconst.dti import TensorModel

duration1 = time() - t1
        
# print('BIAC006'+' DTI duration %.3f' % (duration1,))

lpca_path=outpath+'/LPCA_' + id + '_nii4D.nii.gz'
if path.exists(lpca_path):
    print('File already exists; Skipping LPCA denoising (path: ' + lpca_path + ')' )
else:
    print('Beginning LPCA denoising for: '+ id + '.  (Expected result: ' + lpca_path + ')' )
    t = time()
        
    print(data.shape)
    data2=data
    sigma1 = pca_noise_estimate(data2, gtab, correct_bias=True, smooth=1)
    print("Sigma estimation time", time() - t)
    
    #lpca
    t = time()
    denoised_arr = localpca(data2, sigma=sigma1, patch_radius=2,pca_method='svd', tau_factor=2.3)
    save_nifti(lpca_path, denoised_arr, affine)
    print("Time taken for local PCA denoising", -t + time())

