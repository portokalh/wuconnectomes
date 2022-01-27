
import numpy as np
from tract_manager import create_tracts, diff_preprocessing, tract_connectome_analysis, get_diffusionattributes
from Daemonprocess import MyPool
import multiprocessing as mp
import os
from file_tools import mkcdir
from time import time
from argument_tools import parse_arguments
import sys
import glob

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

subjects = ["S02666","S02670","S02686","S02654", "S02686", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781",
            "S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02938", "S02939", "S02967", "S02987", "S02987",
            "S03010", "S03033", "S03034", "S03045", "S02524","S02535","S02690","S02715","S02771","S02804","S02812",
            "S02817", "S02840","S02871","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028",
            "S03048", "S03069"]
#removed S02666, S02670, S02686 for now
#"S02230" "S02690" "S02804" these subjects are strange, to investigatei
#02490 has been since confirmed to be incomplete and can be safely ignored, as it is the same as 02491 subject wise
#02804 seems to have streamlines in negative voxel space, thought I fixed that through pruning, skip FOR NOW
#02926, 02954 03017 03048 has an out of bounds error, will skip FOR NOW
"""
subjfolder = glob.glob(os.path.join(datapath, "*" + identifier + "*"))[0]
subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
for bxhfile in subjbxh:
    bxhtype = checkbxh(bxhfile, False)
    if bxhtype == "dwi":
        dwipath = bxhfile.replace(".bxh", ".nii.gz")
        break
"""

outpath = "/Volumes/Data/Badea/ADdecode.01/Analysis/"
outpath="/mnt/paros_MRI/jacques/AD_Decode/Analysis/"
outpath = '/Volumes/Data/Badea/Lab/jacques/ADDecode_test_connectomes/'
diff_preprocessed = os.path.join(outpath, "DWI")
trkpath = os.path.join(outpath, "TRK_MPCA_fixed")
trkpath = os.path.join(outpath, "TRK_MPCA_100")

mkcdir([outpath, diff_preprocessed, trkpath])
masktype = "FA"
masktype = "T1"
masktype = "subjspace"
subjects_all = glob.glob(os.path.join(diff_preprocessed,'*subjspace_dwi*.nii.gz'))

for subject in subjects_all:
    subject_name = os.path.basename(subject)
    subjects.append(subject_name[:6])
print(subjects)
removed_list = ['S02821']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

subjects = ['S01912']
subject_processes, function_processes = parse_arguments(sys.argv,subjects)

stepsize = 2


ratio = 100
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
#classifier = ["FA", "binary"]
classifier = "binary"
labelslist = []
bvec_orient = [1,2,-3]
vol_b0 = [0,1,2]

dwi_results = []
donelist = []
notdonelist = []
createmask = masktype
symmetric = True
inclusive = True

denoise = "coreg"
savefa = True
make_connectomes = True

classifiertype = "FA"
classifiertype = "binary"
brainmask = "subjspace"
labeltype='lrordered'

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"


#atlas_legends = None
atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
#atlas_legends = outpath + "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

figspath = os.path.join(outpath, "Figures_MPCA"+inclusive_str+symmetric_str+saved_streamlines)
mkcdir(figspath)

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

    tract_results = pool.starmap_async(create_tracts, [(diff_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                          verbose) for subject
                                                       in subjects]).get()
    if make_connectomes:
        tract_results = pool.starmap_async(tract_connectome_analysis, [(diff_preprocessed, trkpath, str_identifier, figspath,
                                                                       subject, atlas_legends, bvec_orient, inclusive,
                                                                       function_processes, forcestart, picklesave, labeltype, symmetric, verbose)
                                                                     for subject in subjects]).get()
    pool.close()
else:
    for subject in subjects:
        tract_results.append(
            create_tracts(diff_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                          verbose))
        #get_diffusionattributes(diff_preprocessed, diff_preprocessed, subject, str_identifier, vol_b0, ratio, bvec_orient,
        #                        createmask, overwrite, verbose)
        if make_connectomes:
            tract_results.append(tract_connectome_analysis(diff_preprocessed, trkpath, str_identifier, figspath, subject,
                                                           atlas_legends, bvec_orient,  brainmask, inclusive,
                                                           function_processes, forcestart, picklesave, labeltype, symmetric, verbose))
    print(tract_results)

# dwi_results.append(diff_preprocessing(datapath, diff_preprocessed, subject, bvec_orient, denoise, savefa,
#                                     function_processes, createmask, vol_b0, verbose)) ##Unnecessary, replaced by SAMBA_prep
#dwi_results = pool.starmap_async(diff_preprocessing, [(datapath, diff_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
#                                 createmask, vol_b0, verbose) for subject in subjects]).get()
