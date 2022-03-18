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


gunniespath = "~/gunnies/"
mainpath="/Volumes/Data/Badea/Lab/"
diffpath = "/Volumes/Data/Badea/Lab/RaulChavezValdez/"

#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/mnt/munin6/Badea/Lab/mouse/Chavez_series/diffusion_prep_locale/"
SAMBA_inputs_folder = "/mnt/munin6/Badea/Lab/Chavez_prep/"
shortcuts_all_folder = "/mnt/munin6/Badea/Lab/mouse/Chavez_symlink_pool_allfiles/"

if SAMBA_inputs_folder is not None:
    mkcdir(SAMBA_inputs_folder)
if shortcuts_all_folder is not None:
    mkcdir(shortcuts_all_folder)
mkcdir(outpath)
subjects = ['C_20220124_001', 'C_20220124_002', 'C_20220124_003', 'C_20220124_004', 'C_20220124_005', 'C_20220124_006', 'C_20220124_007']

"""
subjects_folders = glob.glob(os.path.join(diffpath,'diffusion*/'))
subjects = []
for subject_folder in subjects_folders:
    subjects.append(subject_folder.split('diffusion')[1][:6])
"""
removed_list = []
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

print(subjects)

# subjects = ['N58610', 'N58612', 'N58613']

subject_processes, function_processes = parse_arguments(sys.argv, subjects)

# subject 'N58610' retired, weird? to investigate
proc_subjn = ""
denoise = "None"
#denoise = "mpca"
recenter = 0
proc_name = "diffusion_prep_" + proc_subjn
cleanup = True
masking = "median_5"
makebtables = False
gettranspose = False
ref = "coreg"
copybtables = True
verbose = True
transpose = None
overwrite = False

# btables=["extract","copy","None"]
btables = "extract"
# Neither copy nor extract are fully functioning right now, for now the bvec extractor from extractdiffdirs works
# go back to this if ANY issue with bvals/bvecs
# extract is as the name implies here to extract the bvals/bvecs from the files around subject data
# copy takes one known good file for bval and bvec and copies it over to all subjects
if btables == "extract":
    for subject in subjects:
        # outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = outpath + "_" + subject
        writeformat = "tab"
        writeformat = "dsi"
        outpathsubj = os.path.join(outpath, proc_name + subject)
        mkcdir(outpathsubj)
        outpathrelative = os.path.join(outpath, "relative_orientation.txt")
        newoutpathrelative = os.path.join(outpathsubj, "relative_orientation.txt")
        shutil.copy(outpathrelative, newoutpathrelative)
        fbvals, fbvecs = extractbvals(diffpath, subject, outpath=outpathsubj, writeformat=writeformat,
                                      overwrite=False)  # extractbvals_research
        # fbvals, fbvecs = rewrite_subject_bvalues(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
elif btables == "copy":
    for subject in subjects:
        # outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = os.path.join(outpath, proc_name + subject)
        outpathbval = os.path.join(outpathsubj, proc_subjn + subject + "_bvals.txt")
        outpathbvec = os.path.join(outpathsubj, proc_subjn + subject + "_bvecs.txt")
        outpathrelative = os.path.join(outpath, "relative_orientation.txt")
        newoutpathrelative = os.path.join(outpathsubj, "relative_orientation.txt")
        mkcdir(outpathsubj)
        shutil.copy(outpathrelative, newoutpathrelative)
        if not os.path.exists(outpathbval) or not os.path.exists(outpathbvec) or overwrite:
            mkcdir(outpathsubj)
            writeformat = "tab"
            writeformat = "dsi"
            bvals = glob.glob(os.path.join(outpath, "*N58408*bvals*.txt"))
            bvecs = glob.glob(os.path.join(outpath, "*N58408*bvec*.txt"))
            if np.size(bvals) > 0 and np.size(bvecs) > 0:
                shutil.copy(bvals[0], outpathbval)
                shutil.copy(bvecs[0], outpathbvec)

proc_name = "diffusion_prep_"

max_processors = 20
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them
nominal_bval = 2401
verbose = True
function_processes = np.int(max_processors / subject_processes)
results = []
if subject_processes > 1:
    if function_processes > 1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)
    results = pool.starmap_async(launch_preprocessing, [
        launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                             shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose, overwrite,
                             denoise, recenter,
                             recenter, verbose) for subject in subjects]).get()
else:
    for subject in subjects:
        max_size = 0
        subjectpath = glob.glob(os.path.join(os.path.join(diffpath,subject + "*")))[0]
        subject_outpath = os.path.join(outpath, 'diffusion_prep_' + proc_subjn + subject)
        max_file = largerfile(subjectpath,identifier=".nii")
        if os.path.exists(os.path.join(subject_outpath, f'{subject}_subjspace_fa.nii.gz')) and not overwrite:
            print(f'already did subject {subject}')
        elif os.path.exists(os.path.join('/Volumes/Badea/Lab/APOE_symlink_pool/',
                                         f'{subject}_subjspace_coreg.nii.gz')) and not overwrite:
            print(
                f'Could not find subject {subject} in main diffusion folder but result was found in SAMBA prep folder')
        #elif os.path.exists(os.path.join(
        #        '/Volumes/Data/Badea/Lab/mouse/VBM_20APOE01_chass_symmetric3_allAPOE-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n32_i5/reg_images/',
        #        f'{subject}_rd_to_MDT.nii.gz')) and not overwrite:
        #    print(f'Could not find subject {subject} in main diff folder OR samba init but was in results of SAMBA')
        else:
            launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                                 shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose,
                                 overwrite, denoise, recenter, verbose)
        # results.append(launch_preprocessing(subject, max_file, outpath))


