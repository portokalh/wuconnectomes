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
from transform_handler import get_transpose
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


#gunniespath = "/mnt/clustertmp/common/rja20_dev/gunnies/"
gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
diffpath = "/Volumes/dusom_civm-atlas/20.abb.15/research/"
outpath = "/Volumes/Data/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"
#outpath = "/mnt/munin6/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"
bonusshortcutfolder = "/Volumes/Data/Badea/Lab/jacques/APOE_series/19abb14/"
#bonusshortcutfolder = "/mnt/munin6/Badea/Lab/19abb14/"
bonusshortcutfolder = None

subjects = ["N58214","N58215","N58216","N58217","N58218","N58219","N58221","N58222","N58223","N58224","N58225","N58226","N58228",
            "N58229","N58230","N58231","N58232","N58633","N58634","N58635","N58636","N58649","N58650","N58651","N58653","N58654",
            'N58408', 'N58610', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
            'N58612', 'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
            'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
            'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
            'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
            'N58613', 'N58732', 'N58516', 'N58813', 'N58402']

subjects = ["N58214","N58215","N58216","N58217","N58218","N58219","N58221","N58222","N58223","N58224","N58225","N58226","N58228",
            "N58229","N58230","N58231","N58232","N58633","N58634","N58635","N58636","N58649","N58650","N58651","N58653","N58654"]
#subjects = ['N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302', 'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608', 'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305', 'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396', 'N58613', 'N58732', 'N58516', 'N58402']

atlas = "/mnt/munin6/Badea/Lab/atlases/chass_symmetric3/chass_symmetric3_DWI.nii.gz"

proc_subjn=""
denoise="None"
recenter=0
proc_name ="diffusion_prep_"+proc_subjn
cleanup = True
masking = "median"
makebtables = False
gettranspose=False
copybtables = True
verbose=True
transpose=None
overwrite=False
ref="coreg"
#btables=["extract","copy","None"]
btables="none"

if btables == "extract":
    for subject in subjects:
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite_b=False
        proc_name = "diffusion_prep_"  # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline
        outpath_subj = os.path.join(outpath,proc_name+subject)
        mkcdir(outpath_subj)
        fbvals, fbvecs = extractbvals_research(dwipath, subject, outpath=outpath_subj, fix=False, writeformat=writeformat, overwrite=overwrite_b)

max_processors = 1
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()
subject_processes = np.size(subjects)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them
nominal_bval=4000
verbose=True
function_processes = np.int(max_processors/subject_processes)
results=[]
if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    results = pool.starmap_async(launch_preprocessing, [(proc_subjn+subject,
                                                         largerfile(glob.glob(os.path.join(os.path.join(diffpath, "diffusion*" + subject + "*")))[0]),
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
        launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, bonusshortcutfolder,
                             gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter,
                             verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))

