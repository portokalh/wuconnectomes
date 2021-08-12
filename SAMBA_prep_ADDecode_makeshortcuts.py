import numpy as np
from tract_manager import create_tracts
import multiprocessing as mp
from Daemonprocess import MyPool
import glob
import os
from bvec_handler import extractbvals, extractbvals_research, rewrite_subject_bvalues, fix_bvals_bvecs
from time import time
import shutil
from diffusion_preprocessing import launch_preprocessing
from file_tools import mkcdir, largerfile
from img_transform_exec import get_transpose
import shutil

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

gunniespath = "~/gunnies/"
mainpath="/mnt/munin6/Badea/ADdecode.01/"

#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale/"
bonusshortcutfolder = "/mnt/munin6/Badea/Lab/mouse/ADDeccode_symlink_pool/"
diffpath = mainpath + "Data/Anat"

subjects = ["02654", "02690", "02720", "02737", "02753", "02765", "02781", "02802", "02804", "02813", "02817", "02840", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987", "02987", "03010", "03017", "03033", "03034", "03045", "03048"]

atlas = None

overwrite=False
cleanup = True
atlas = None
gettranspose=False
if gettranspose:
    transpose = get_transpose(atlas)

proc_subjn="S"
proc_name ="diffusion_prep_"+proc_subjn
denoise = "lpca"
masking = "bet"
overwrite=False
cleanup = True
atlas = None
recenter=0
gettranspose=False
verbose=True
nominal_bval=1000
if gettranspose:
    transpose = get_transpose(atlas)
ref = "md"

btables="None"

if btables == "extract":
    for subject in subjects:
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite_b=False
        proc_name = "diffusion_prep_"  # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline
        outpath_subj = os.path.join(outpath,proc_name+subject)
        mkcdir(outpath_subj)
        fbvals, fbvecs = extractbvals_research(diffpath, subject, outpath=outpath_subj, fix=False, writeformat=writeformat, overwrite=overwrite_b)

max_processors = 1
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()
subject_processes = np.size(subjects)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them
nominal_bval=1000
verbose=True
function_processes = np.int(max_processors/subject_processes)
results=[]
if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    results = pool.starmap_async(launch_preprocessing, [(proc_subjn+subject,
                                                         largerfile(glob.glob(os.path.join(os.path.join(diffpath, "*" + subject + "*")))[0]),
                                                         outpath, cleanup, nominal_bval, bonusshortcutfolder,
                                                         gunniespath, function_processes, masking, ref, transpose,
                                                         overwrite, denoise, recenter, verbose)
                                                        for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        subjectpath = glob.glob(os.path.join(os.path.join(outpath, "diffusion*"+subject+"*")))[0]
        max_file=largerfile(subjectpath)
        max_file= os.path.join(subjectpath, "nii4D_"+subject+".nii.gz")
        print(max_file)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        launch_preprocessing(subject, max_file, outpath, cleanup, nominal_bval, bonusshortcutfolder,
                                                         gunniespath, function_processes, atlas, transpose,
                             overwrite, denoise, verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))

