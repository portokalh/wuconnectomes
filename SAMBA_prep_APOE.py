import numpy as np
from tract_manager import create_tracts
import multiprocessing as mp
from Daemonprocess import MyPool
import glob
import os
from bvec_handler import extractbvals, rewrite_subject_bvalues, fix_bvals_bvecs
from time import time
import shutil
from mouse_diffusion_preprocessing import launch_preprocessing
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
dwipath = "/Volumes/dusom_civm-atlas/20.abb.15/research/"
dwipath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/"
subject = "58214"
#outpath = None
#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/"

subjects = ["58302", "58303", "58305", "58309", "58310", "58344","58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58303", "58305", "58309", "58310", "58344","58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58303"]
subjects = ["58305", "58309", "58310", "58344","58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
#subjects = ["58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
#subjects = ["58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58612","58613","58706","58708","58712"]
subjects = ["58606","58608","58610","58611"]
subjects = ["58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510"]
subjects = ["58611"]
subjects = ["58612","58613","58706","58708","58712"]
subjects = ["58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611"]
subjects = ["58613","58706"]
subjects = ["58302", "58303", "58305", "58309", "58310", "58344","58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
subjects = ["58613","58706","58708","58712"]
subjects = ["58612"]

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
        launch_preprocessing(subject, max_file, outpath)
        #results.append(launch_preprocessing(subject, max_file, outpath))

"""
subjectlist = ["58215","58216","58217","58218","58219","58221","58222","58223","58224","58225","58226","58228","58229","58230","58231","58232","58633","58634","58635","58636","58649","58650","58651","58653","58654"]
for subj in subjectlist:
    fbvals_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    fbvecs_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
"""
