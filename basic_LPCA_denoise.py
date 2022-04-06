#!/usr/bin/env python3
#$ -l h_vmem=50000M,vf=50000M
#$ -M ${USER}@duke.edu 
#$ -m ea 
#$ -o /mnt/munin6/Badea/Lab/mouse/sinha_sbatch/slurm-$JOB_ID.out 
#$ -e /mnt/munin6/Badea/Lab/mouse/sinha_sbatch/slurm-$JOB_ID.out
#$ -N ${1}_LPCA_denoising


"""
#runno=sys.argv[1] # switching to more generic "id"
#bval_folder=sys.argv[3] # Need to have gtab handling that 1: takes in a single value, & 2: supports DSI Studio btables
# However, for now we'll stick with the bval/bvec pair to keep things moving forward.

#id = 58214
#fdwi="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale//diffusion_prep_58214/nii4D_N58214.nii.gz"
#bval_or_bvec_or_btable="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale//diffusion_prep_58214/58214_bvecs.txt"
#outpath="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale//diffusion_prep_58214"


#find other pointers@
#https://github.com/nipy/nipype/blob/fbf2c35f533b7805ca93c742006472e0809d8d03/nipype/workflows/dmri/mrtrix/diffusion.py
#to do: coreg/eddy correction/bias field (if not part of denoising already)
#https://github.com/nipy/dipy/blob/349e6b181ac89f333f07146751a2295b732b5c04/scratch/very_scratch/registration_example.py
"""

from os import path
import sys
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table

#from dipy.denoise.localpca import localpca
from denoise_processes import localpca, mppca
from dipy.denoise.pca_noise_estimate import pca_noise_estimate
import os
from time import time
from dipy.segment.mask import median_otsu
"""
from dipy.denoise.non_local_means import non_local_means
from dipy.denoise.adaptive_soft_matching import adaptive_soft_matching
from dipy.denoise.nlmeans import nlmeans
from dipy.denoise.noise_estimate import estimate_sigma
import matplotlib.pyplot as plt
from dipy.tracking.utils import random_seeds_from_mask
from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.streamline import Streamlines
from dipy.reconst.shm import CsaOdfModel
import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from dipy.workflows.denoise import NLMeansFlow
from dipy.denoise.denspeed import nlmeans_3d
from dipy.io.streamline import save_trk
import os.path
"""

"""
id=sys.argv[1]
fdwi=sys.argv[2]
bval_or_bvec_or_btable=sys.argv[3]
outpath=sys.argv[4]
"""

def basic_LPCA_denoise_func(id,fdwi,bval_or_bvec_or_btable,outpath, processes=1, denoise="LPCA", verbose=False):

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
    masked_path = fdwi
    if no_masking:
        # masked_path=outpath + id + '_nii4D_masked.nii.gz'
        if path.exists(masked_path):
            maskdata, affine, vox_size = load_nifti(masked_path, return_voxsize=True)
        elif path.exists(masked_path.replace(".nii.gz",".nii")):
            masked_path=masked_path.replace(".nii.gz",".nii")
            maskdata, affine, vox_size = load_nifti(masked_path, return_voxsize=True)
        else:
            txt = f'Could not find masked path at {masked_path}'
            raise Exception(txt)
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
        outdir = os.path.dirname(outpath)
        mask_path=os.path.join(outdir,id +'_mask.nii.gz')
        save_nifti(mask_path, mask.astype(np.ubyte), affine)

    data=maskdata

    from dipy.reconst.dti import TensorModel

    duration1 = time() - t1

    # print('BIAC006'+' DTI duration %.3f' % (duration1,))

    if denoise.lower() == 'lpca':
        id=str(id)
        if path.exists(outpath):
            print('File already exists; Skipping LPCA denoising (path: ' + outpath + ')' )
        else:
            print('Beginning LPCA denoising for: '+ id + '. (Expected result: ' + outpath + ')' )
            t = time()

            print(data.shape)
            data2=data
            print(masked_path)
            print(gtab)
            sigma1 = pca_noise_estimate(data2, gtab, correct_bias=True, smooth=1)
            print("Sigma estimation time", time() - t)

            #lpca
            t = time()
            denoised_arr = localpca(data2, sigma=sigma1, patch_radius=2,pca_method='svd', tau_factor=2.3,
                                    processes=processes, verbose=verbose)
            save_nifti(outpath, denoised_arr, affine)
            print("Time taken for local PCA denoising", -t + time())
    elif denoise.lower() == 'mpca':
        id=str(id)
        #mpca_path=outpath+'/MPCA_' + id + '_nii4D.nii.gz'
        if path.exists(outpath):
            print('File already exists; Skipping LPCA denoising (path: ' + outpath + ')' )
        else:
            print('Beginning MPCA denoising for: '+ id + '. (Expected result: ' + outpath + ')' )
            t = time()

            print(data.shape)
            data2=data
            print(masked_path)
            print(gtab)
            #sigma1 = pca_noise_estimate(data2, gtab, correct_bias=True, smooth=1)
            #print("Sigma estimation time", time() - t)
            #lpca
            t = time()
            denoised_arr = mppca(data2, patch_radius=2, return_sigma=False, processes=processes, verbose=verbose)
            save_nifti(outpath, denoised_arr, affine)
            print("Time taken for Marcenko-Pastur denoising", -t + time())
