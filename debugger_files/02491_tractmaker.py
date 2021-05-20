import os
from nifti_handler import getfa, getdwidata, getlabelmask, move_bvals, getmask
import numpy as np
from tract_manager import create_tracts
from dipy.segment.mask import median_otsu
from dipy.io.image import load_nifti, save_nifti
from diff_preprocessing import dwi_to_mask, denoise_pick
from dif_to_trk import make_tensorfit, QCSA_tractmake
from Daemonprocess import MyPool
import multiprocessing as mp

def orient_to_str(bvec_orient):
    mystr="_"
    for i in np.arange(3):
        if np.abs(bvec_orient[i]) == 1:
            if bvec_orient[i]<0:
                mystr = mystr+"mx"
            else:
                mystr = mystr+"px"
        if np.abs(bvec_orient[i]) == 2:
            if bvec_orient[i] < 0:
                mystr = mystr + "my"
            else:
                mystr = mystr + "py"
        if np.abs(bvec_orient[i])==3:
            if bvec_orient[i]<0:
                mystr = mystr+"mz"
            else:
                mystr = mystr+"pz"
    return mystr


masktype = "FA"
masktype = "T1"
masktype = "dwi"
subject_processes = 1
function_processes = 10
ratio = 10
get_params = None
doprune = True
labelslist = []
figspath = "/Volumes/Data/Badea/ADdecode.01/Analysis/"
outpath = "/Volumes/Data/Badea/ADdecode.01/Analysis/"
dwipath = "/Volumes/Data/Badea/ADdecode.01/Data/Anat/20210216_02491/bia6_02491_003.nii.gz"
stepsize = 0.5
vol_b0 = [0,1,2]
subject = "02491"
strproperty = "_FA"
verbose = True
overwrite = False
outpathtrk = os.path.join(outpath,subject + strproperty + '_pruned.trk')
subject_processes = 1

outpathmask = os.path.join(outpath, subject)
#data, affine, gtab, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, None)
#mask, _ = dwi_to_mask(data, affine, outpathmask, makefig=False, vol_idx=vol_b0, median_radius=5, numpass=6,
#                      dilate=2)

#if masktype == "FA":
#    data, affine, gtab, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, None)
#    outpathbmfa, mask = make_tensorfit(data, mask, gtab, affine, subject, outpath=dwipath, verbose=verbose)
if masktype == "dwi":
    outpathmask = os.path.join(outpath, subject)
    data, affine, gtab, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, None)
    mask, _ = dwi_to_mask(data, affine, outpathmask, makefig=False, vol_idx=vol_b0, median_radius=5, numpass=6,
                          dilate=2)
elif masktype == "T1":
    #bet bia6_02491_40006.nii.gz 02491.nii.gz -m -o -f 0.4
    #mv 02491_mask.nii.gz 02491_T1_binary_mask.nii.gz
    mask, affinemask = getmask(outpath,subject,"T1",verbose)

import itertools
bvec_orient1 = (np.array(list(itertools.permutations([1, 2, 3]))))
bvec_orient2 = [elm*[-1, 1, 1] for elm in bvec_orient1]
bvec_orient3 = [elm*[1, -1, 1] for elm in bvec_orient1]
bvec_orient4 = [elm*[1, 1, -1] for elm in bvec_orient1]

bvec_orient_list = np.concatenate((bvec_orient1, bvec_orient2, bvec_orient3, bvec_orient4))

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(create_tracts, [(dwipath, outpath, subject, figspath, stepsize, function_processes,
                                                        orient_to_str(bvec_orient), ratio, masktype, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, verbose) for bvec_orient in bvec_orient_list]).get()
    pool.close()
else:
    txtfile = os.path.join(outpath, subject + "_params.txt")
    for bvec_orient in bvec_orient_list:
        tract_results = []
        print(bvec_orient)
        str_identifier = orient_to_str(bvec_orient)
        str_identifier = strproperty + str_identifier
        tract_results.append(create_tracts(dwipath, outpath, subject, figspath, stepsize, function_processes,
                          str_identifier, ratio, masktype, 'FA', labelslist, bvec_orient, doprune, overwrite, get_params,
                          verbose))
        print(tract_results)

    """
    with open(txtfile, 'a') as f:
        for item in tract_results:
            f.write("Subject %s with %s %s %s \n" % (item[0],str(bvec_orient[0]),str(bvec_orient[1]),str(bvec_orient[2])))
            f.write("Num tracts: %s \n" % item[2][0])
            f.write("Min tract length: %s \n" % item[2][1])
            f.write("Max tract length: %s \n" % item[2][2])
            f.write("Average tract length: %s \n" % item[2][3])
            f.write("Standard deviancy tract length: %s \n" % item[2][4])
    """
