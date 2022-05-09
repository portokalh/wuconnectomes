
"""
Created by Jacques Stout
Part of the DTC pipeline, mostly handles dwi files before calculating trk.
Tries to create masks, determines the parameters of a denoising request, handles fa files, etc
"""

import matplotlib.pyplot as plt
from dipy.core.histeq import histeq
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.segment.mask import median_otsu
import os
from denoise_processes import mppca
from dipy.denoise.gibbs import gibbs_removal
from time import time
from figures_handler import denoise_fig
import glob
from mask_handler import applymask_array
import nibabel as nib
from computer_nav import checkfile_exists_remote, load_nifti_remote, save_nifti_remote

#from os.path import join as pjoin
#from dipy.data import get_fnames
#from nifti_handler import getdwidata

def string_inclusion(string_option,allowed_strings,option_name):
    "checks if option string is part of the allowed list of strings for that option"
    try:
        string_option=string_option.lower()
    except AttributeError:
        if string_option is None:
            #raise Warning(option_name + " stated as None value, option will not be implemented")
            print(option_name + " stated as None value, option will not be implemented")
            return
        else:
            raise AttributeError('Unrecognized value for ' + option_name)

    if not any(string_option == x for x in allowed_strings):
        raise ValueError(string_option + " is an unrecognized string, please check your input for " + option_name)
    if string_option == "none":
        print(option_name + " stated as None value, option will not be implemented")


def dwi_to_mask(data, subject, affine, outpath, masking='median', makefig=False, vol_idx=None, median_radius = 5,
                numpass=6, dilate = 2, forcestart = False, header = None, verbose = False, sftp=None):

    data = np.squeeze(data)
    binarymask_path = os.path.join(outpath, subject + '_dwi_binary_mask.nii.gz')
    maskeddwi_path = os.path.join(outpath, subject + '_dwi_mask.nii.gz')
    if not checkfile_exists_remote(binarymask_path,sftp) and not checkfile_exists_remote(maskeddwi_path, sftp) \
            and not forcestart:
        if masking == 'median':
            b0_mask, mask = median_otsu(data, vol_idx=vol_idx, median_radius=median_radius, numpass=numpass,
                                        dilate=dilate)
        if masking == 'extract':
            if np.size(np.shape(data)) == 3:
                mask=data>0
            if np.size(np.shape(data)) == 4:
                mask=data[:,:,:,0]>0
            mask = mask.astype(np.float32)
            b0_mask = applymask_array(data,mask)

        if verbose:
            txt = f"Creating binarymask at {binarymask_path} and masked data at {maskeddwi_path}"
            print(txt)
        if header is None:
            save_nifti(binarymask_path, mask, affine)
            save_nifti(maskeddwi_path, b0_mask.astype(np.float32), affine)
        else:
            binarymask_nii = nib.Nifti1Image(mask, affine, header)
            save_nifti_remote(binarymask_nii, binarymask_path, sftp=sftp)
            #maskeddwi_nii = nib.Nifti1Image(b0_mask, affine, header)
            #save_nifti_remote(maskeddwi_nii, maskeddwi_path, sftp=sftp)
    else:
        mask = load_nifti_remote(binarymask_path,sftp=sftp)
        mask = mask[0]
        b0_mask = load_nifti_remote(maskeddwi_path,sftp=sftp)
        b0_mask = b0_mask[0]

    if makefig:
        sli = data.shape[2] // 2
        if len(b0_mask.shape) ==4:
            b0_mask_2 = b0_mask[:,:,:,0]
        else:
            b0_mask_2 = b0_mask
        if len(data.shape) ==4:
            data = data[:,:,:,0]
        plt.figure('Brain segmentation')
        plt.subplot(1, 2, 1).set_axis_off()
        plt.imshow(histeq(data[:, :, sli].astype('float')).T,
                   cmap='gray', origin='lower')

        plt.subplot(1, 2, 2).set_axis_off()
        plt.imshow(histeq(b0_mask_2[:, :, sli].astype('float')).T,
                   cmap='gray', origin='lower')
        plt.savefig(outpath + 'median_otsu.png')

    return(mask.astype(np.float32), b0_mask.astype(np.float32))

def dwi_to_mask_old(data, subject, affine, outpath, masking='median', makefig=False, vol_idx=None, median_radius = 5,
                numpass=6, dilate = 2, forcestart = False, header = None, verbose = False, sftp=None):

    data = np.squeeze(data)
    binarymask_path = os.path.join(outpath, subject + '_dwi_binary_mask.nii.gz')
    maskeddwi_path = os.path.join(outpath, subject + '_dwi_mask.nii.gz')
    if not os.path.exists(binarymask_path) and not os.path.exists(maskeddwi_path) and not forcestart:
        if masking == 'median':
            b0_mask, mask = median_otsu(data, vol_idx=vol_idx, median_radius=median_radius, numpass=numpass,
                                        dilate=dilate)
        if masking == 'extract':
            if np.size(np.shape(data)) == 3:
                mask=data>0
            if np.size(np.shape(data)) == 4:
                mask=data[:,:,:,0]>0
            mask = mask.astype(np.float32)
            b0_mask = applymask_array(data,mask)

        if verbose:
            txt = f"Creating binarymask at {binarymask_path} and masked data at {maskeddwi_path}"
            print(txt)
        if header is None:
            save_nifti(binarymask_path, mask, affine)
            save_nifti(maskeddwi_path, b0_mask.astype(np.float32), affine)
        else:
            binarymask_nii = nib.Nifti1Image(mask, affine, header)
            nib.save(binarymask_nii, binarymask_path)
            maskeddwi_nii = nib.Nifti1Image(b0_mask, affine, header)
            nib.save(maskeddwi_nii, maskeddwi_path)
    else:
        mask = load_nifti(binarymask)
        mask = mask[0]
        b0_mask = load_nifti(maskeddwi)
        b0_mask = b0_mask[0]

    if makefig:
        sli = data.shape[2] // 2
        if len(b0_mask.shape) ==4:
            b0_mask_2 = b0_mask[:,:,:,0]
        else:
            b0_mask_2 = b0_mask
        if len(data.shape) ==4:
            data = data[:,:,:,0]
        plt.figure('Brain segmentation')
        plt.subplot(1, 2, 1).set_axis_off()
        plt.imshow(histeq(data[:, :, sli].astype('float')).T,
                   cmap='gray', origin='lower')

        plt.subplot(1, 2, 2).set_axis_off()
        plt.imshow(histeq(b0_mask_2[:, :, sli].astype('float')).T,
                   cmap='gray', origin='lower')
        plt.savefig(outpath + 'median_otsu.png')

    return(mask.astype(np.float32), b0_mask.astype(np.float32))


def check_for_fa(outpath, subject, getdata=False):
    #Checks for fa files ('bmfa') in specified outpath folder. Returns with the path
    #whether it exists or not, and the fa nifti if specified to do so
    if os.path.isdir(outpath):
        outpathglobfa = os.path.join(outpath, subject + '_*fa.nii.gz')
    elif os.path.isfile(outpath):
        outpathglobfa = os.path.join(os.path.dirname(outpath), subject + '_*fa.nii.gz')
    outpathfa = glob.glob(outpathglobfa)
    if np.size(outpathfa) == 1:
        outpathfa=outpathfa[0]
        if getdata is True:
            fa = load_nifti(outpathfa)
            return outpathfa, True, fa
        else:
            return outpathfa, True, None
    elif np.size(outpathfa) == 0:
        outpathglobfa.replace("*fa","bmfa")
        return outpathfa, False, None

def make_tensorfit(data,mask,gtab,affine,subject,outpath, overwrite=False, forcestart = False, verbose=None):
    #Given dwi data, a mask, and other relevant information, creates the fa and saves it to outpath, unless
    #if it already exists, in which case it simply returns the fa

    from dipy.reconst.dti import TensorModel
    outpathbmfa, exists, _ = check_for_fa(outpath, subject, getdata = False)
    if exists and not forcestart:
        fa = load_nifti(outpathbmfa)
        fa_array = fa[0]
        if verbose:
            txt = "FA already computed at " + outpathbmfa
            print(txt)
        return outpathbmfa, fa_array
    else:
        if verbose:
            print('Calculating the tensor model from bval/bvec values of ', subject)
        tensor_model = TensorModel(gtab)

        t1 = time()
        if len(mask.shape) == 4:
            mask = mask[...,0]
        tensor_fit = tensor_model.fit(data, mask)

        duration1 = time() - t1
        if verbose:
            print(subject + ' DTI duration %.3f' % (duration1,))

        save_nifti(outpathbmfa, tensor_fit.fa, affine)
        if verbose:
            print('Saving subject'+ subject+ ' at ' + outpathbmfa)

        return outpathbmfa, tensor_fit.fa


def denoise_pick(data, affine, hdr, outpath, mask, type_denoise='macenko', processes=1, savedenoise=True, verbose=False,
                 forcestart=False, datareturn=False, display=None):
    allowed_strings = ['mpca', 'yes', 'all', 'gibbs', 'none', 'macenko']
    string_inclusion(type_denoise, allowed_strings, "type_denoise")

    if type_denoise == 'macenko' or type_denoise == 'mpca' or type_denoise == 'yes' or type_denoise == 'all':
        type_denoise = '_mpca_'
    if type_denoise == 'gibbs':
        type_denoise = "_gibbs"
    if type_denoise is None or type_denoise == 'none':
        type_denoise = "_"

    outpath_denoise = outpath + type_denoise + 'dwi.nii.gz'
    if os.path.exists(outpath_denoise) and not forcestart:
        if verbose:
            txt = "Denoising already done at " + outpath_denoise
            print(txt)
        if datareturn:
            data = load_nifti(outpath_denoise)
    else:
        if type_denoise == '_mpca_':
            # data, snr = marcenko_denoise(data, False, verbose=verbose)
            t = time()
            denoised_arr, sigma = mppca(data, patch_radius=2, return_sigma=True, processes=processes, verbose=verbose)
            save_nifti(outpath_denoise, denoised_arr, affine, hdr=hdr)
            if verbose:
                txt = ("Saved image at " + outpath_denoise)
                print(txt)
            mask = np.array(mask, dtype=bool)
            mean_sigma = np.mean(sigma[mask])
            b0 = denoised_arr[..., 0]

            mean_signal = np.mean(b0[mask])

            snr = mean_signal / mean_sigma
            if verbose:
                print("Time taken for local MP-PCA ", -t +
                      time())
                print("The SNR of the b0 image appears to be at " + str(snr))
            if display:
                denoise_fig(data, denoised_arr, type='macenko')
            data = denoised_arr

        if type_denoise == 'gibbs':
            outpath_denoise = outpath + '_gibbs.nii.gz'
            if os.path.exists(outpath_denoise) and not forcestart:
                if verbose:
                    txt = "Denoising already done at " + outpath_denoise
                    print(txt)
                if datareturn:
                    data = load_nifti(outpath_denoise)
            t = time()
            data_corrected = gibbs_removal(data, slice_axis=2)
            save_nifti(outpath_denoise, denoised_arr, affine, hdr=hdr)
            if verbose:
                print("Time taken for the gibbs removal ", - t + time())
            if display:
                denoise_fig(data, data_corrected, type='gibbs')

            data = data_corrected

        if type_denoise == "_":
            print('No denoising was done')
            save_nifti(outpath_denoise, data, affine, hdr=hdr)

    return data, outpath_denoise
