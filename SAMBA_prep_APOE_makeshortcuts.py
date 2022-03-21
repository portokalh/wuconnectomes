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

#gunniespath = "/mnt/clustertmp/common/rja20_dev/gunnies/"
#gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
#diffpath = "/Volumes/dusom_civm-atlas/20.abb.15/research/"
#outpath = "/Volumes/Data/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"

gunniespath = ""
outpath = "/mnt/munin6/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/"
#bonusshortcutfolder = "/Volumes/Data/Badea/Lab/jacques/APOE_series/19abb14/"
#bonusshortcutfolder = "/mnt/munin6/Badea/Lab/APOE_symlink_pool/"
#bonusshortcutfolder = None
SAMBA_inputs_folder = "/mnt/munin6/Badea/Lab/19abb14/"
shortcuts_all_folder = "/mnt/munin6/Badea/Lab/mouse/APOE_symlink_pool_allfiles/"

mkcdir([SAMBA_inputs_folder, shortcuts_all_folder])


subjects = ["N58214","N58215","N58216","N58217","N58218","N58219","N58221","N58222","N58223","N58224","N58225","N58226","N58228",
            "N58229","N58230","N58231","N58232","N58633","N58634","N58635","N58636","N58649","N58650","N58651","N58653","N58654",
            'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
            'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
            'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
            'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
            'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
            'N58613', 'N58732', 'N58516', 'N58813', 'N58402']
#58610, N58612 removed
subjects = ['N58408', 'N58398', 'N58935', 'N58714', 'N58740', 'N58477', 'N59003', 'N58734', 'N58309', 'N58792', 'N58819', 'N58302', 'N58909', 'N58784', 'N58919', 'N58706', 'N58889', 'N58361', 'N58355', 'N58712', 'N58790', 'N59010', 'N58859', 'N58917', 'N58606', 'N58815', 'N58997', 'N58350', 'N58999', 'N58881', 'N58608', 'N58853', 'N58779', 'N58995', 'N58500', 'N58604', 'N58749', 'N58877', 'N58883', 'N58510', 'N58885', 'N58906', 'N58394', 'N58821', 'N58855', 'N58346', 'N58861', 'N58344', 'N58857', 'N58788', 'N58305', 'N58514', 'N58851', 'N58794', 'N58733', 'N58655', 'N58887', 'N58735', 'N58310', 'N58879', 'N58400', 'N58708', 'N58780', 'N58512', 'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58829', 'N58913', 'N58745', 'N58831', 'N58406', 'N58359', 'N58742', 'N58396', 'N58941', 'N58516', 'N58813', 'N58402']
removed_list = []

subjects_fpath = glob.glob("/mnt/munin6/Badea/Lab/mouse/APOE_series/diffusion_prep_locale/diffusion_prep*")
subjects = []
for subject in subjects_fpath:
    subject_fname = os.path.basename(subject)
    subjects.append(subject_fname.split('diffusion_prep_')[1])

subjects = ['N58952', 'N58995', 'N58997', 'N58999', 'N59003', 'N59010', 'N59022', 'N59026', 'N59033', 'N59035',
            'N59039', 'N59041', 'N59065', 'N59066', 'N59072', 'N59076', 'N59078', 'N59080', 'N59097', 'N59099',
            'N59109', 'N59116', 'N59118', 'N59120']

removed_list = []
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

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

    results = pool.starmap_async(launch_preprocessing, [launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                                 shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose, overwrite, denoise, recenter,
                              verbose) for subject in subjects]).get()
else:
    for subject in subjects:
        max_size=0
        print(os.path.join(os.path.join(outpath, "diffusion*"+subject+"*")))
        subjectpath = glob.glob(os.path.join(os.path.join(outpath, "diffusion*"+subject+"*")))[0]
        max_file=largerfile(subjectpath)
        max_file= os.path.join(subjectpath, "nii4D_"+subject+".nii.gz")
        print(max_file)
        #command = gunniespath + "mouse_diffusion_preprocessing.bash"+ f" {subject} {max_file} {outpath}"
        if os.path.exists(os.path.join(shortcuts_all_folder,f'{proc_subjn + subject}_fa.nii.gz')) and os.path.exists(os.path.join(SAMBA_inputs_folder, f'{proc_subjn + subject}_fa.nii.gz')):
            print(f'already did subject {proc_subjn + subject}')
        else:
            #print('notyet')
            launch_preprocessing(proc_subjn + subject, max_file, outpath, cleanup, nominal_bval, SAMBA_inputs_folder,
                                 shortcuts_all_folder, gunniespath, function_processes, masking, ref, transpose, overwrite, denoise,
                             recenter, verbose)
