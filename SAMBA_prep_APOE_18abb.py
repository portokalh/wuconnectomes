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
from bvec_handler import orient_to_str


gunniespath = "/Users/jas/bass/gitfolder/gunnies/"
diffpath = "/Volumes/dusom_civm-atlas/18.abb.11/research/"
#diffpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/research/"
#outpath = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_locale/"
outpath = "/Volumes/Data/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"

#bonusshortcutfolder = "/Volumes/Data/Badea/Lab/19abb14/"

SAMBA_inputs_folder = "/Volumes/Data/Badea/Lab/19abb14/"
shortcuts_all_folder = "/Volumes/Data/Badea/Lab/mouse/APOE_symlink_pool_allfiles/"
shortcuts_all_folder = None

"""
subjects = ['N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302', 'N58612'
    , 'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608', 'N58779', 'N58500',
            'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305', 'N58514', 'N58794',
            'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303',
            'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396', 'N58613', 'N58732', 'N58516', 'N58813', 'N58402']
            
subjects = ['N58408', 'N59072', 'N58398', 'N58935', 'N58714', 'N58740', 'N58477', 'N59003', 'N58734', 'N58309', 'N58792'
, 'N58819', 'N58302', 'N59078', 'N59116', 'N58909', 'N58784', 'N58919', 'N58706', 'N58889', 'N58361', 'N58355', 'N59066'
, 'N58712', 'N58790', 'N59010', 'N58859', 'N58917', 'N58606', 'N58815', 'N59118', 'N58997', 'N58350', 'N59022', 'N58999'
, 'N58881', 'N59026', 'N58608', 'N58853', 'N58779', 'N58995', 'N58500', 'N58604', 'N58749', 'N58877', 'N58883', 'N59109'
, 'N59120', 'N58510', 'N58885', 'N58906', 'N59065', 'N58394', 'N58821', 'N58855', 'N58346', 'N58861', 'N58344', 'N59099'
, 'N58857', 'N58788', 'N58305', 'N58514', 'N58851', 'N59076', 'N59097', 'N58794', 'N58733', 'N58655', 'N58887', 'N58735'
, 'N58310', 'N59035', 'N58879', 'N58400', 'N59041', 'N58952', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303', 'N58404'
, 'N58751', 'N58611', 'N58829', 'N58913', 'N58745', 'N58831', 'N58406', 'N58359', 'N58742', 'N58396', 'N58941', 'N59033'
, 'N58516', 'N59080', 'N58813', 'N59039', 'N58402']


subjects = ['N58305', 'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303',
            'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396', 'N58613', 'N58732', 'N58516', 'N58813', 'N58402']

subjects = ['N58712', 'N58790', 'N58606', 'N58350', 'N58608', 'N58779', 'N58500',
            'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788']
"""


subjects = ['N58302', 'N58612'
    , 'N58784', 'N58706', 'N58361', 'N58355']
subjects = ['N57452','N57456','N57458','N57462']


subjects_folders = glob.glob(os.path.join(diffpath,'diffusion*/'))
subjects = []
for subject_folder in subjects_folders:
    subjects.append(subject_folder.split('diffusion')[1][:6])



subjects = sorted(subjects)
subjects = subjects[76:]
#subjects.reverse()
#subjects = subjects[:8]


removed_list = ['N58610', 'N58612', 'N58613','N58732']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

print(subjects)

#subjects = ['N58610', 'N58612', 'N58613']
subjects = ['N58610']
subject_processes, function_processes = parse_arguments(sys.argv, subjects)

#subject 'N58610' retired, weird? to investigate
proc_subjn=""
denoise="None"
recenter=0
proc_name ="diffusion_prep_"+proc_subjn
cleanup = True
masking = "median"
makebtables = False
gettranspose=False
ref = "coreg"
copybtables = True
verbose=True
transpose=None
overwrite=False

#btables=["extract","copy","None"]
btables="extract"
#Neither copy nor extract are fully functioning right now, for now the bvec extractor from extractdiffdirs works
#go back to this if ANY issue with bvals/bvecs
#extract is as the name implies here to extract the bvals/bvecs from the files around subject data
#copy takes one known good file for bval and bvec and copies it over to all subjects
if btables=="extract":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = os.path.join(outpath,proc_name+subject)
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        fbvals, fbvecs = extractbvals(diffpath, subject, outpath=outpathsubj, writeformat=writeformat, overwrite=overwrite) #extractbvals_research
        #fbvals, fbvecs = rewrite_subject_bvalues(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
elif btables=="copy":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = os.path.join(outpath,proc_name+subject)
        outpathbval= os.path.join(outpathsubj, proc_subjn + subject+"_bvals.txt")
        outpathbvec= os.path.join(outpathsubj, proc_subjn + subject+"_bvecs.txt")
        outpathrelative = os.path.join(outpath, "relative_orientation.txt")
        newoutpathrelative= os.path.join(outpathsubj, "relative_orientation.txt")
        mkcdir(outpathsubj)
        shutil.copy(outpathrelative, newoutpathrelative)
        if not os.path.exists(outpathbval) or not os.path.exists(outpathbvec) or overwrite:
            mkcdir(outpathsubj)
            writeformat="tab"
            writeformat="dsi"
            overwrite=True
            bvals = glob.glob(os.path.join(outpath, "*N58408*bvals*.txt"))
            bvecs = glob.glob(os.path.join(outpath, "*N58408*bvec*.txt"))
            if np.size(bvals)>0 and np.size(bvecs)>0:
                shutil.copy(bvals[0], outpathbval)
                shutil.copy(bvecs[0], outpathbvec)


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
        #writeformat="dsi"
        #fix_bvals_bvecs(bval_file, fbvecs=bvec_file, outpath=subjectpath, writeformat=writeformat)
        if not os.path.exists(new_bval_file):
            shutil.copyfile(bval_file,new_bval_file)
        if not os.path.exists(new_bvec_file):
            shutil.copyfile(bvec_file,new_bvec_file)

max_processors = 20
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

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
    results = pool.starmap_async(launch_preprocessing, [launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                                 shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter,
                             recenter, verbose) for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        subjectpath = glob.glob(os.path.join(os.path.join(diffpath, "diffusion*"+subject+"*")))[0]
        subject_outpath = os.path.join(outpath, 'diffusion_prep_' + proc_subjn + subject)
        max_file=largerfile(subjectpath)
        if os.path.exists(os.path.join(subject_outpath, f'{subject}_subjspace_fa.nii.gz')) and not overwrite:
            print(f'already did subject {subject}')
        elif os.path.exists(os.path.join('/Volumes/Badea/Lab/APOE_symlink_pool/', f'{subject}_subjspace_coreg.nii.gz')) and not overwrite:
            print(f'Could not find subject {subject} in main diffusion folder but result was found in SAMBA prep folder')
        elif os.path.exists(os.path.join('/Volumes/Data/Badea/Lab/mouse/VBM_20APOE01_chass_symmetric3_allAPOE-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n32_i5/reg_images/',f'{subject}_rd_to_MDT.nii.gz')) and not overwrite:
            print(f'Could not find subject {subject} in main diff folder OR samba init but was in results of SAMBA')
        else:
            launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                                 shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose,
                                 overwrite, denoise, recenter, verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))

"""

(subj, raw_nii, outpath, cleanup=False, nominal_bval=4000, SAMBA_inputs_folder=None,
                         shortcuts_all_folder = None, gunniespath="~/gunnies/", processes=1, masking="bet", ref=None,
                         transpose=None, overwrite=False, denoise='None', recenter=0, verbose=False)

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
