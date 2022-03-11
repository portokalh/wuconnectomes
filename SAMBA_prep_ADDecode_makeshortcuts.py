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
from bvec_handler import orient_to_str

gunniespath = "~/gunnies/"
mainpath="/mnt/munin6/Badea/ADdecode.01/"

#outpath = "/Users/alex/jacques/APOE_temp"
outpath = "/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale/"
SAMBA_inputs_folder = "/mnt/munin6/Badea/Lab/mouse/ADDeccode_symlink_pool/"
shortcuts_all_folder = "/mnt/munin6/Badea/Lab/human/ADDeccode_symlink_pool_allfiles/"
diffpath = mainpath + "Data/Anat"
mkcdir([SAMBA_inputs_folder, shortcuts_all_folder])

subjects = ['01912', '02110', '02224', '02227', '02230', '02231', '02266', '02289', '02320', '02361', '02363', '02373', '02386', '02390', '02402', '02410', '02421', '02424', '02446', '02451', '02469', '02473', '02485', '02491', '02490', '02506', '02523', '02524', '02535', '02654', '02666', '02670', '02686', '02690', '02695', '02715', '02720', '02737', '02745', '02753', '02765', '02771', '02781', '02802', '02804', '02813', '02812', '02817', '02840', '02842', '02871', '02877', '02898', '02926', '02938', '02939', '02954', '02967', '02987', '03010', '03017', '03028', '03033', '03034', '03045', '03048', '03069', '03225', '03265', '03293', '03308', '03321', '03343','03350','03378','03391','03394']

subjects = ['02524', '02535', '02654', '02666', '02670', '02686', '02690', '02695', '02715', '02720', '02737', '02745', '02753', '02765', '02771', '02781', '02802', '02804', '02813', '02812', '02817', '02840', '02842', '02871', '02877', '02898', '02926', '02938', '02939', '02954', '02967', '02987', '03010', '03017', '03028', '03033', '03034', '03045', '03048', '03069', '03225', '03265', '03293', '03308', '03321', '03343','03350','03378','03391','03394']
subjects = ['02227']
subjects = ["02695","02686", "03010", "02670", "02666", "02654"]
removed_list = ['02230','02231','02490','02523','02745','02266','02289','02320','02361','02363','02373','02386','02390','S02402']
subjects = ['01912']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

atlas = None

overwrite=False
cleanup = True
atlas = None
gettranspose=False
if gettranspose:
    transpose = get_transpose(atlas)

proc_subjn="S"
proc_name ="diffusion_prep_"+proc_subjn
denoise = "mpca"
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
transpose=None
btables="None"
if btables=="extract":
    for subject in subjects:
        #outpathsubj = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/diffusion_prep_58214/"
        outpathsubj = outpath + "_" + subject
        writeformat="tab"
        writeformat="dsi"
        overwrite=True
        fbvals, fbvecs = extractbvals(diffpath, subject, outpath=outpath, writeformat=writeformat, overwrite=overwrite)
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
            overwrite=True
            bvals = glob.glob(os.path.join(outpath, "*bvals*.txt"))
            bvecs = glob.glob(os.path.join(outpath, "*bvec*.txt"))
            if np.size(bvals)>0 and np.size(bvecs)>0:
                shutil.copy(bvals[0], outpathbval)
                shutil.copy(bvecs[0], outpathbvec)

overwrite=False
#quickfix was here
max_processors = 1
if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()
subject_processes = np.size(subjects)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them
nominal_bval=1000
masking = 'bet'
denoise = 'mpca'
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
                                                         outpath, cleanup, nominal_bval, SAMBA_inputs_folder,shortcuts_all_folder,
                                                         gunniespath, function_processes, masking, ref, transpose,
                                                         overwrite, denoise, recenter, verbose)
                                                        for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        subjectpath = glob.glob(os.path.join(os.path.join(outpath, "diffusion*"+subject+"*")))[0]
        print(subjectpath)
        max_file=largerfile(subjectpath)
        max_file= os.path.join(subjectpath, "nii4D_"+proc_subjn+subject+".nii.gz")
        print(max_file)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        subject_f = proc_subjn + subject

        if os.path.exists(os.path.join('/mnt/munin6/Badea/Lab/mouse/ADDeccode_symlink_pool/', f'{subject_f}_subjspace_coreg.nii.gz')):
            print(f'Could not find subject {subject_f} in main diffusion folder but result was found in SAMBA prep folder')
        #elif os.path.exists(os.path.join('/mnt/munin6/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/reg_images/',f'{subject_f}_rd_to_MDT.nii.gz')):
        #    print(f'Could not find subject {subject_f} in main diff folder OR samba init but was in results of SAMBA')
        else:
            launch_preprocessing(proc_subjn+subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,shortcuts_all_folder,
                                                         gunniespath, function_processes, masking, ref, transpose,
                             overwrite, denoise, verbose)
        #results.append(launch_preprocessing(subject, max_file, outpath))

