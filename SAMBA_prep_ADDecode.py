import numpy as np
from tract_manager import create_tracts
from Daemonprocess import MyPool
import multiprocessing as mp
import glob
import os
import sys
from bvec_handler import extractbvals, rewrite_subject_bvalues, fix_bvals_bvecs
from diffusion_preprocessing import launch_preprocessing
from file_tools import mkcdir, largerfile
from transform_handler import get_transpose
import shutil
from argument_tools import parse_arguments

munin=False
if munin:
    gunniespath = "~/wuconnectomes/gunnies"
    mainpath = "/mnt/munin6/Badea/ADdecode.01/"
    outpath = "/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale_mpca/"
    bonusshortcutfolder = "/mnt/munin6/Badea/Lab/mouse/ADDeccode_symlink_pool/"
else:
    gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
    mainpath="/Volumes/Data/Badea/ADdecode.01/"
    #outpath = "/Users/alex/jacques/APOE_temp"
    outpath = "/Volumes/Data/Badea/Lab/human/AD_Decode/diffusion_prep_locale_mpca/"
    bonusshortcutfolder = "/Volumes/Data/Badea/Lab/mouse/ADDeccode_symlink_pool/"
    bonusshortcutfolder = None

diffpath = os.path.join(mainpath, "Data","Anat")


#gunniespath = "~/gunnies/"
#mainpath = "/mnt/munin6/Badea/ADdecode.01/"
#outpath = "/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale/"
#bonusshortcutfolder = "/mnt/munin6/Badea/Lab/mouse/ADDeccode_symlink_pool/"

mkcdir(outpath)
subjects = ["02654", "02690", "02720", "02737", "02745", "02753", "02765", "02771", "02781", "02802", "02804", "02812", "02813", "02817", "02840", "02842", "02871", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987", "02987", "03010", "03017", "03028", "03033", "03034", "03045", "03048"]

subjects = ["02654", "02690", "02720", "02737", "02753", "02765", "02781", "02802", "02804", "02813", "02817", "02840", "02842", "02871", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987", "02987", "03010", "03017", "03028", "03033", "03034", "03045", "03048"]
subjects = ["02654", "02666", "02670", "02686", "02690", "02695", "02715", "02720", "02737", "02753", "02765", "02771", "02781", "02802", "02804", "02813", "02817", "02840", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987", "02987", "03010", "03017", "03033", "03034", "03045", "03048"]
#subjects = ["02666"]
#subjects = ["02871"]
#subjects = ["02842", "02812", "02871", "02715", "02771","03069"]
subjects = ["02654"]

#subjects = ["02871", "02877", "02898", "02926", "02938", "02939", "02954", "02967", "02987", "02987", "03010", "03017", "03028", "03033", "03034", "03045", "03048"]
#02745 was not fully done, discount
#02771 has 21 images in 4D space even though there should be 23?
#"02812", 02871 is a strange subject, to investigate
#02842, 03028 has apparently a 92 stack ? to investigate

subject_processes, function_processes = parse_arguments(sys.argv,subjects)

proc_subjn="S"
proc_name ="diffusion_prep_"+proc_subjn
denoise = "mpca"
masking = "bet"
overwrite=False
cleanup = True
atlas = None
gettranspose=False
verbose=True
nominal_bval=1000
if gettranspose:
    transpose = get_transpose(atlas)
ref = "md"
recenter=0

transpose=None

#btables=["extract","copy","None"]
btables="copy"
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
        fbvals, fbvecs = extractbvals(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=True)
        #fbvals, fbvecs = rewrite_subject_bvalues(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
elif btables=="copy":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = os.path.join(outpath,proc_name+subject)
        outpathbval= os.path.join(outpathsubj, proc_subjn + subject+"_bvals.txt")
        outpathbvec= os.path.join(outpathsubj, proc_subjn + subject+"_bvecs.txt")
        if not os.path.exists(outpathbval) or not os.path.exists(outpathbvec) or overwrite:
            mkcdir(outpathsubj)
            writeformat="tab"
            writeformat="dsi"
            bvals = glob.glob(os.path.join(outpath, "*bvals*.txt"))
            bvecs = glob.glob(os.path.join(outpath, "*bvec*.txt"))
            if np.size(bvals)>0 and np.size(bvecs)>0:
                shutil.copy(bvals[0], outpathbval)
                shutil.copy(bvecs[0], outpathbvec)
#quickfix was here

# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

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
        subjectpath = glob.glob(os.path.join(os.path.join(diffpath, "*" + subject + "*")))[0]
        max_file=largerfile(subjectpath)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        #max_file="/Volumes/Data/Badea/ADdecode.01/Data/Anat/20210522_02842/bia6_02842_003.nii.gz"
        #launch_preprocessing(subject, max_file, outpath, nominal_bval=1000, shortcutpath=shortcutpath, bonusshortcutfolder = bonusshortcutfolder, gunniespath="/Users/alex/bass/gitfolder/gunnies/")
        #max_file = '/Volumes/Data/Badea/ADdecode.01/Data/Anat/20210522_02842/bia6_02842_003.nii.gz'
        launch_preprocessing(proc_subjn+subject, max_file, outpath, cleanup, nominal_bval, bonusshortcutfolder,
         gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter, verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))

"""
subjectlist = ["58215","58216","58217","58218","58219","58221","58222","58223","58224","58225","58226","58228","58229","58230","58231","58232","58633","58634","58635","58636","58649","58650","58651","58653","58654"]
for subj in subjectlist:
    fbvals_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    fbvecs_new = fbvals.replace("58214", subj)
    shutil.copyfile(fbvals, fbvals_new)
    
    
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
