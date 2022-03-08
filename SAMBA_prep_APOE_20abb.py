import numpy as np
from tract_manager import create_tracts
import multiprocessing as mp
from Daemonprocess import MyPool
import glob
import os, sys
from bvec_handler import writebfiles, extractbvals, extractbvals_research, rewrite_subject_bvalues, fix_bvals_bvecs
from time import time
import shutil
from diffusion_preprocessing import launch_preprocessing
from file_tools import mkcdir, largerfile
import shutil
from argument_tools import parse_arguments

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

bonusshortcutfolder = "/Volumes/Data/Badea/Lab/19abb14/"
bonusshortcutfolder = None

gunniespath = "/Users/alex/bass/gitfolder/gunnies/"
diffpath = "/Volumes/dusom_civm-atlas/20.abb.15/research/"
outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/"
outpath = "/Volumes/Data/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"

#bonusshortcutfolder = "/Volumes/Data/Badea/Lab/jacques/APOE_series/19abb14/"

subjects = ["N58214","N58215","N58216","N58217","N58218","N58219","N58221","N58222","N58223","N58224","N58225","N58226","N58228",
            "N58229","N58230","N58231","N58232","N58633","N58634","N58635","N58636","N58649","N58650","N58651","N58653","N58654"]
#atlas = "/Volumes/Data/Badea/Lab/atlases/chass_symmetric3/chass_symmetric3_DWI.nii.gz"

subject_processes, function_processes = parse_arguments(sys.argv, subjects)

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
btables=""


#Neither copy nor extract are fully functioning right now, for now the bvec extractor from extractdiffdirs works
#go back to this if ANY issue with bvals/bvecs
#extract is as the name implies here to extract the bvals/bvecs from the files around subject data
#copy takes one known good file for bval and bvec and copies it over to all subjects
if btables=="extract":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        fbvals, fbvecs = extractbvals(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite) #extractbvals_research
        #fbvals, fbvecs = rewrite_subject_bvalues(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
elif btables=="copy":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = os.path.join(outpath,proc_name+subject)
        outpathbval= os.path.join(outpathsubj, proc_subjn + subject+"_bvals.txt")
        outpathbvec= os.path.join(outpathsubj, proc_subjn + subject+"_bvecs.txt")
        outpathrelative = os.path.join(outpath, "relative_orientation.txt")
        newoutpathrelative= os.path.join(outpathsubj, "relative_orientation.txt")
        shutil.copy(outpathrelative, newoutpathrelative)
        if not os.path.exists(outpathbval) or not os.path.exists(outpathbvec) or overwrite:
            mkcdir(outpathsubj)
            writeformat="tab"
            writeformat="dsi"
            overwrite=True
            bvals = glob.glob(os.path.join(outpath, "*N58214*bvals*.txt"))
            bvecs = glob.glob(os.path.join(outpath, "*N58214*bvec*.txt"))
            if np.size(bvals)>0 and np.size(bvecs)>0:
                shutil.copy(bvals[0], outpathbval)
                shutil.copy(bvecs[0], outpathbvec)

# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them
nominal_bval=4000
verbose=True
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
        subjectpath = glob.glob(os.path.join(os.path.join(diffpath, "diffusion*"+subject+"*")))[0]
        max_file=largerfile(subjectpath)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        launch_preprocessing(proc_subjn+subject, max_file, outpath, cleanup, nominal_bval, bonusshortcutfolder,
         gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter, verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))
