
"""
Created by Jacques Stout
Functions for applying and creating masks from MRI acquisitions
"""


import nibabel as nib
from dipy.io.image import load_nifti, save_nifti
from dipy.segment.mask import median_otsu
import numpy as np

def applymask_samespace(file, mask, outpath=None):
    #note: there should definitely be a bet option which would probably end up faster, but for some insane reason there is no
    #obvious documentation on it, so this will do for now -_-
    img_nii= nib.load(file)
    mask_nii = nib.load(mask)
    img_data = img_nii.get_fdata()
    img_data_new = img_data
    mask_data = mask_nii.get_fdata()
    img_shape = img_nii.shape
    dims = np.size(img_shape)
    if (dims == 3 or dims == 4) and mask_nii.shape[0:3] == img_nii.shape[0:3]:
        for i in np.arange(img_shape[0]):
            for j in np.arange(img_shape[1]):
                for k in np.arange(img_shape[2]):
                    if mask_data[i, j, k] == 1:
                        if dims == 3:
                            img_data_new[i, j, k] = img_data[i, j, k]
                        else:
                            for l in range(img_shape[3]):
                                img_data_new[i, j, k, l] = img_data[i, j, k, l]
                    else:
                        if dims == 3:
                            img_data_new[i, j, k] = 0
                        else:
                            for l in range(img_shape[3]):
                                img_data_new[i, j, k, l] = 0
    elif dims not in [3,4]:
        raise TypeError("Could not interpret the dimensions of the entering image")
    elif mask_nii.shape[0:3] != img_nii.shape[0:3]:
        raise TypeError("The shape of the mask and the image are not equal, therefore readjustments are needed")
    img_nii_new = nib.Nifti1Image(img_data_new, img_nii.affine, img_nii.header)
    if outpath is None:
        outpath = file
    nib.save(img_nii_new, outpath)


def median_mask_make(inpath, outpath, outpathmask=None, median_radius=4, numpass=4,binary_dilation=None):

    if outpathmask is None:
        outpathmask=outpath.replace(".nii","_mask.nii")
    data, affine = load_nifti(inpath)
    data = np.squeeze(data)
    data_masked, mask = median_otsu(data, median_radius=median_radius, numpass=numpass, dilate=binary_dilation)
    save_nifti(outpath, data_masked.astype(np.float32), affine)
    save_nifti(outpathmask, mask.astype(np.float32), affine)
