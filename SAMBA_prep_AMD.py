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
    mainpath = "/mnt/munin/Whitson/BrainChAMD.01/"
    outpath = "/mnt/munin6/Badea/Lab/human/AMD/diffusion_prep_locale/"
    bonusshortcutfolder = "/mnt/munin6/Badea/Lab/mouse/whiston_symlink_pool_test/"
else:
    gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
    mainpath="/Volumes/Data/Whitson/BrainChAMD.01/"
    outpath = "/Volumes/Data/Badea/Lab/human/AMD/diffusion_prep_locale/"
    bonusshortcutfolder = "/Volumes/Data/Badea/Lab/mouse/whiston_symlink_pool_test/"
    bonusshortcutfolder = None

diffpath = os.path.join(mainpath, "Data","Anat")

mkcdir(outpath)
subjects = ["29056", "26578", "29060", "26637", "29264", "26765", "29225", "26660", "29304", "26890", "29556", "26862", "29410", "26966", "29403", "26841", "21593", "27126", "29618", "27111", "29627", "27164", "29502", "27100", "27381", "21836", "27391", "21850", "27495", "21729", "27488", "21915", "27682", "21956", "27686", "22331", "28208", "21990", "28955", "29878", "27719", "22102", "27841", "22101", "27842", "22228", "28029", "22140", "27852", "22276", "27999", "22369", "28115", "22644", "28308", "22574", "28377", "22368", "28325", "22320", "28182", "22898", "28748", "22683", "28373", "22536", "28433", "22825", "28662", "22864", "28698", "23143", "28861", "23157", "28820", "23028", "29002", "23210", "29020", "23309", "29161", "26841", "26862", "26949", "26966", "27100", "27126", "27163", "27246", "27488", "27682", "27686", "27719", "27841", "27842", "27852", "27869", "27999", "28029", "28068", "28208", "28262", "28325", "28820", "28856", "28869", "28955", "29002", "29044", "29089", "29127", "29161", "29242", "29254", "26578", "26637", "26660", "26745", "26765", "26850", "26880", "26890", "26958", "26974", "27017", "27111", "27164", "27381", "27391", "27495", "27610", "27640", "27680", "27778", "27982", "28115", "28308", "28338", "28373", "28377", "28433", "28437", "28463", "28532", "28662", "28698", "28748", "28809", "28857", "28861", "29013", "29020", "29025"]
subjects = ["27999"]

subject_processes, function_processes = parse_arguments(sys.argv,subjects)

proc_subjn="H"
proc_name ="diffusion_prep_"+proc_subjn
denoise = "mpca"
masking = "bet"
overwrite=True
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
        launch_preprocessing(proc_subjn+subject, max_file, outpath, cleanup, nominal_bval, bonusshortcutfolder,
         gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter, verbose)
