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


gunniespath = "/Users/alex/bass/gitfolder/gunnies/"
dwipath = "/Volumes/dusom_civm-atlas/18.abb.11/research/"
#dwipath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/"
subject = "58214"
#outpath = None
#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/"
outpath = "/Volumes/Data/Badea/Lab/jacques/APOE_series/diffusion_prep_locale/"

bonusshortcutfolder = "/Volumes/Data/Badea/Lab/19abb14/"

subjects = ['N58408', 'N58610', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302', 'N58612', 'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608', 'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305', 'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396', 'N58613', 'N58732', 'N58516', 'N58813', 'N58402']

cleanup = True

makebtables = False
copybtables = True
if makebtables:
    for subject in subjects:
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        proc_name = "diffusion_prep_"  # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline
        outpath_subj = os.path.join(outpath,proc_name+subject)
        mkcdir(outpath_subj)
        fbvals, fbvecs = extractbvals_research(dwipath, subject, outpath=outpath_subj, writeformat=writeformat, overwrite=overwrite)


bval_file="/Volumes/Data/Badea/Lab/jacques/APOE_series/diffusion_prep_locale/diffusion_prep_N58408/N58408_bvals.txt"
bvec_file="/Volumes/Data/Badea/Lab/jacques/APOE_series/diffusion_prep_locale/diffusion_prep_N58408/N58408_bvecs.txt"
proc_name ="diffusion_prep_"
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

max_processors = 10
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

    results = pool.starmap_async(launch_preprocessing, [(subject,
                                                         largerfile(glob.glob(os.path.join(os.path.join(dwipath, "diffusion*"+subject+"*")))[0]),
                                                         outpath, cleanup, nominal_bval,
                                                         bonusshortcutfolder, gunniespath, function_processes, verbose)
                                                        for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        subjectpath = glob.glob(os.path.join(os.path.join(dwipath, "diffusion*"+subject+"*")))[0]
        max_file=largerfile(subjectpath)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        launch_preprocessing(subject, max_file, outpath, bonusshortcutfolder=bonusshortcutfolder, gunniespath=gunniespath, processes=function_processes, cleanup=cleanup)
        #results.append(launch_preprocessing(subject, max_file, outpath))

"""
subjectlist = ["58215","58216","58217","58218","58219","58221","58222","58223","58224","58225","58226","58228","58229","58230","58231","58232","58633","58634","58635","58636","58649","58650","58651","58653","58654"]
for subj in subjectlist:
    fbvals_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    fbvecs_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)


makebtables=False
if makebtables:
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/diffusion_prep_"+subject
        mkcdir(outpath)
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        fbvals, fbvecs = extractbvals(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
        #fbvals, fbvecs = rewrite_subject_bvalues(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)

quickfix = False
if quickfix:
    bval_file="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/diffusionN58302dsi_studio/N58302_bvals.txt"
    bvec_file="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/diffusionN58302dsi_studio/N58302_bvecs.txt"
    bval_file, bvec_file = fix_bvals_bvecs(bval_file, bvec_file, outpath= outpath, identifier="", format="dsi")

bval_file="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/N58302_bvals.txt"
bvec_file="/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/diffusionN58302dsi_studio/N58302_bvecs.txt"
proc_name ="diffusion_prep_"

results=[]
copybtables=False
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

"""
