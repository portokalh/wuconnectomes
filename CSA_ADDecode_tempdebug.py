import os
from nifti_handler import getfa, getdwidata, getlabelmask, move_bvals, getmask
import numpy as np
from tract_manager import create_tracts, dwi_preprocessing, tract_connectome_analysis, get_diffusionattributes
from dipy.segment.mask import median_otsu
from dipy.io.image import load_nifti, save_nifti
from diff_preprocessing import dwi_to_mask, denoise_pick
from dif_to_trk import make_tensorfit, QCSA_tractmake
from bvec_handler import checkbxh
from Daemonprocess import MyPool
import multiprocessing as mp
import glob
import os
from time import time

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

def mkcdir(folderpaths):
    for folderpath in folderpaths:
        if not os.path.exists(folderpath):
            os.mkdir(folderpath)

datapath = "/Volumes/Data/Badea/ADdecode.01/Data/Anat/"
subjects = ["01912", "02110", "02224", "02227", "02230", "02231", "02266", "02289", "02320", "02361", "02363", "02373", "02386", "02390", "02402", "02410", "02421", "02424", "02446", "02451", "02469", "02473", "02485", "02490", "02491", "02506"]
#subjects = ["02227"]
#subjects = ["02402", "02410", "02421", "02424", "02446", "02451", "02469", "02473", "02485", "02490", "02491", "02506"]
#subjects = ["02402", "02410", "02421"]
#subjects = ["02424", "02446", "02451", "02469", "02473", "02485", "02490", "02491", "02506"]
subjects = ["02231", "02266", "02289"]
subjects = ["02490", "02491", "02506"]
subjects = ["01912", "02110", "02224", "02227", "02231", "02266", "02289", "02320", "02361", "02363", "02373", "02386", "02390", "02402", "02410", "02421", "02424", "02446", "02451", "02469", "02473", "02485", "02491", "02506"]

subjects = ["01912", "02110", "02224", "02227", "02231", "02266"]
#"02230" "02490" these subjects are strange, to investigate

"""
subjfolder = glob.glob(os.path.join(datapath, "*" + identifier + "*"))[0]
subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
for bxhfile in subjbxh:
    bxhtype = checkbxh(bxhfile, False)
    if bxhtype == "dwi":
        dwipath = bxhfile.replace(".bxh", ".nii.gz")
        break
"""

outpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/AD_Decode_temp/"
figspath = os.path.join(outpath, "Figures")
dwi_preprocessed = os.path.join(outpath, "DWI")
trkpath = os.path.join(outpath, "TRK")

mkcdir([outpath, figspath, dwi_preprocessed, trkpath])
masktype = "FA"
masktype = "T1"
masktype = "dwi"


max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

subject_processes = np.size(subjects)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
function_processes = np.int(max_processors / subject_processes)

"""
if masktype == "dwi":
    outpathmask = os.path.join(outpath, subject)
    data, affine, gtab, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, None)
    mask, _ = dwi_to_mask(data, affine, outpathmask, makefig=False, vol_idx=vol_b0, median_radius=5, numpass=6,
                          dilate=2)
elif masktype == "T1":
    #bet bia6_02491_40006.nii.gz 02491.nii.gz -m -o -f 0.4
    #mv 02491_mask.nii.gz 02491_T1_binary_mask.nii.gz
    mask, affinemask = getmask(outpath,subject,"T1",verbose)
"""

stepsize = 2


ratio = 1
if ratio == 1:
    saved_streamlines = "_all"
else:
    saved_streamlines = "_ratio_" + str(ratio)

trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
#str_identifier = '_stepsize_' + str(stepsize) + saved_streamlines+ roistring
str_identifier = '_stepsize_' + str(stepsize).replace(".","_") + saved_streamlines + roistring

duration1=time()
overwrite = False
get_params = False
forcestart = False
if forcestart:
    print("WARNING: FORCESTART EMPLOYED. THIS WILL COPY OVER PREVIOUS DATA")
picklesave = True
verbose = True
get_params = None
doprune = True
classifier = "FA"
labelslist = []
bvec_orient = [1,2,-3]
vol_b0 = [0,1,2]

dwi_results = []
donelist = []
notdonelist = []
createmask = masktype
inclusive = True
denoise = "mpca"
savefa = True
make_connectomes = True

classifiertype = "FA"
classifiertype = "binary"
brainmask = "dwi"

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"


#atlas_legends = None
#atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
atlas_legends = outpath + "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"


if make_connectomes:
    for subject in subjects:
        picklepath_connect = figspath + subject + str_identifier + '_connectomes.p'
        excel_path = figspath + subject + str_identifier + "_connectomes.xlsx"
        if os.path.exists(picklepath_connect) and os.path.exists(excel_path):
            print("The writing of pickle and excel of " + str(subject) + " is already done")
            donelist.append(subject)
        else:
            notdonelist.append(subject)

dwi_results = []
tract_results = []


if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    dwi_results = pool.starmap_async(dwi_preprocessing, [(datapath, dwi_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
                                     createmask, vol_b0, verbose) for subject in subjects]).get()
    tract_results = pool.starmap_async(create_tracts, [(dwi_preprocessed, trkpath, subject, figspath, stepsize, function_processes,
                                                        str_identifier, ratio, masktype, classifier, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, verbose) for subject in subjects]).get()
    if make_connectomes:
        tract_results = pool.starmap_async(tract_connectome_analysis, [(dwi_preprocessed, trkpath, str_identifier, figspath,
                                                                       subject, atlas_legends, bvec_orient, inclusive,
                                                                       function_processes, forcestart, picklesave, verbose)
                                                                     for subject in subjects]).get()
    pool.close()
else:
    for subject in subjects:
        #dwi_results.append(dwi_preprocessing(datapath, dwi_preprocessed, subject, bvec_orient, denoise, savefa,
        #                                     function_processes, createmask, vol_b0, verbose))
        #tract_results.append(
        #    create_tracts(dwi_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
        #                  ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params,
        #                  verbose))
        #get_diffusionattributes(dwi_preprocessed, dwi_preprocessed, subject, str_identifier, vol_b0, ratio, bvec_orient,
        #                        createmask, overwrite, verbose)
        if make_connectomes:
            tract_results.append(tract_connectome_analysis(dwi_preprocessed, trkpath, str_identifier, figspath, subject,
                                                           atlas_legends, bvec_orient,  brainmask, inclusive,
                                                           function_processes, forcestart, picklesave, verbose))
    print(tract_results)
