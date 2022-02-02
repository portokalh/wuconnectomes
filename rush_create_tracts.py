import numpy as np
from dif_to_trk import QCSA_tractmake
import os
from BIAC_tools import send_mail, isempty
from nifti_handler import extract_nii_info, getgtab
import nibabel as nib

input_diff_file = ''
output_trk_file = ''
subject = ''
overwrite = False
verbose = True
bvec_orient = [1,2,-3]
mask_file = ''
mask = nib.load(mask_file).dataobj
masktype = "subjspace"
ratio = 100
step_size = 2
peak_processes=1
get_params=False
doprune = False
#mask, _ = getmask(os.path.dirname(input_diff_file), subject, masktype, verbose)

#check_dif_ratio(outpath, subject, strproperty, ratio)

if os.path.exists(output_trk_file) and overwrite is False:
    print("The tract creation of subject " + subject + " is already done")

if verbose:
    print('Running the ' + subject + ' file')

diff_data, affine, vox_size, header, ref_info = extract_nii_info(input_diff_file, verbose)
gtab = getgtab(os.path.dirname(input_diff_file), subject, bvec_orient)

if np.size(np.shape(mask)) == 1:
    mask = mask[0]
if np.size(np.shape(mask)) == 4:
    mask = mask[:, :, :, 0]
print("Mask shape is " + str(np.shape(mask)))

#if classifier == "FA":
#    outpathbmfa, mask = make_tensorfit(diff_data ,mask ,gtab ,affine ,subject ,outpath=diffpath ,verbose=verbose)

print(verbose)
if verbose:
    txt = ("The QCSA Tractmake is ready to launch for subject " + subject)
    print(txt)
    send_mail(txt ,subject="QCSA main function start")
    print("email sent")

outpathtrk, trkstreamlines, params = QCSA_tractmake(diff_data, affine, vox_size, gtab, mask, masktype, ref_info,
                                                    step_size, peak_processes, output_trk_file, subject, ratio,
                                                    overwrite, get_params, doprune, figspath=None,
                                                    verbose=verbose)
