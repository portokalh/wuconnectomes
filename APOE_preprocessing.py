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
from bvec_handler import extractbvals, rewrite_subject_bvalues
from time import time
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

def mkcdir(folderpaths):
    for folderpath in folderpaths:
        if not os.path.exists(folderpath):
            os.mkdir(folderpath)

dwipath = "/Volumes/dusom_civm-atlas/20.abb.15/research/"
dwipath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/"
subject = "58214"
subjects = ["58302", "58303", "58305", "58309", "58310", "58344","58346","58350","58355","58359","58361","58394","58396","58398","58400","58402","58404","58406","58408","58477","58500","58510","58512","58514","58516","58604","58606","58608","58610","58611","58612","58613","58706","58708","58712"]
for subject in subjects:
    outpath = "/Users/alex/jacques/APOE_temp"
    outpath = None
    outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
    outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/diffusion_prep_"+subject
    mkcdir(outpath)
    writeformat="tab"
    writeformat="dsi_format"
    overwrite=True
    #fbvals, fbvecs = extractbvals(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
    fbvals, fbvecs = rewrite_subject_bvalues(dwipath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)

"""
subjectlist = ["58215","58216","58217","58218","58219","58221","58222","58223","58224","58225","58226","58228","58229","58230","58231","58232","58633","58634","58635","58636","58649","58650","58651","58653","58654"]
for subj in subjectlist:
    fbvals_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    fbvecs_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
"""