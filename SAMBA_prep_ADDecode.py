import numpy as np
from tract_manager import create_tracts, dwi_preprocessing, tract_connectome_analysis, get_diffusionattributes
from dipy.segment.mask import median_otsu
from dipy.io.image import load_nifti, save_nifti
from diff_preprocessing import dwi_to_mask, denoise_pick, make_tensorfit
from dif_to_trk import QCSA_tractmake
from bvec_handler import checkbxh
from Daemonprocess import MyPool
import multiprocessing as mp
import glob
import os
from bvec_handler import extractbvals, rewrite_subject_bvalues, fix_bvals_bvecs
from time import time
import shutil
from mouse_diffusion_preprocessing import launch_preprocessing
from file_tools import mkcdir, largerfile
import shutil

"""
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
"""

gunniespath = "/Users/alex/bass/gitfolder/gunnies/"
mainpath="/Volumes/Data/Badea/ADdecode.01/"
dwipath = mainpath + "Data/Anat"
subject = "58214"
#outpath = None
#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/Volumes/Data/Badea/Lab/human/AD_Decode/diffusion_prep_locale/"
shortcutpath = "/Volumes/Data/Badea/Lab/mouse/ADDeccode_symlink_pool/"
mkcdir(outpath)
#subjects = ["02690", "02720", "02737", "02745", "02753", "02765", "02771"]
subjects = ["02753", "02765", "02771"]
weirdones = ["02745", "02771"]
subjects = ["02654"]
proc_name ="diffusion_prep_"

makebtables=False
if makebtables:
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        fbvals, fbvecs = extractbvals(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
        #fbvals, fbvecs = rewrite_subject_bvalues(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)

quickfix = True
if quickfix:
    bval_file="/Volumes/Data/Badea/ADdecode.01/Analysis/DWI/01912_bvals.txt"
    bvec_file="/Volumes/Data/Badea/ADdecode.01/Analysis/DWI/01912_bvec.txt"
    bval_file, bvec_file = fix_bvals_bvecs(bval_file, bvec_file, outpath= outpath, identifier="", format="dsi")
    results=[]
    copybtables=True
    if copybtables:
        for subject in subjects:
            subjectpath = os.path.join(outpath, proc_name + subject)
            mkcdir(subjectpath)
            #subjectpath = glob.glob(os.path.join(os.path.join(outpath, "*" + subject + "*")))
            #subjectpath = subjectpath[0]
            new_bval_file=os.path.join(subjectpath, subject+"_bvals.txt")
            new_bvec_file=os.path.join(subjectpath, subject+"_bvecs.txt")
            if not os.path.exists(new_bval_file):
                shutil.copyfile(bval_file,new_bval_file)
            if not os.path.exists(new_bvec_file):
                shutil.copyfile(bvec_file,new_bvec_file)

max_processors = 50
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()
subject_processes = np.size(subjects)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)
results=[]
if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    results = pool.starmap_async(launch_preprocessing, [(subject, largerfile(glob.glob(os.path.join(os.path.join(dwipath, "*" + subject + "*")))[0]), outpath) for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        subjectpath = glob.glob(os.path.join(os.path.join(dwipath, "*" + subject + "*")))[0]
        max_file=largerfile(subjectpath)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        launch_preprocessing(subject, max_file, outpath, nominal_bval=1000, shortcutpath=shortcutpath)
        #results.append(launch_preprocessing(subject, max_file, outpath))

"""
subjectlist = ["58215","58216","58217","58218","58219","58221","58222","58223","58224","58225","58226","58228","58229","58230","58231","58232","58633","58634","58635","58636","58649","58650","58651","58653","58654"]
for subj in subjectlist:
    fbvals_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    fbvecs_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
"""
