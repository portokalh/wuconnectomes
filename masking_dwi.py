from os.path import join as pjoin
import numpy as np
from dipy.data import get_fnames
from dipy.io.image import load_nifti, save_nifti
from dipy.segment.mask import median_otsu
from nifti_handler import getdwidata

BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse"
dwipath = BIGGUS_DISKUS + "/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
subject = "H26637"
bvec_orient=[1,2,-3]


data, affine, gtab, mask, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, bvec_orient)
#data_fnames = get_fnames('scil_b0')
#fdwi_data, affine = load_nifti(data_fnames[1])
#data = np.squeeze(data[:,:,:,0])
data = np.squeeze(data)


b0_mask, mask = median_otsu(data, vol_idx = [0,1,2,3], median_radius=5, numpass=6, dilate = 2)

fname = "/Users/alex/whiston_figures_2/H22637_medrad4_sum_numpass6_dilate2"
save_nifti(fname + '_binary_mask.nii.gz', mask.astype(np.float32), affine)
save_nifti(fname + '_mask.nii.gz', b0_mask.astype(np.float32), affine)

import matplotlib.pyplot as plt
from dipy.core.histeq import histeq

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

plt.savefig(fname + 'median_otsu.png')
