
import numpy as np
from tract_manager import create_tracts, diff_preprocessing, tract_connectome_analysis, get_diffusionattributes
from Daemonprocess import MyPool
import multiprocessing as mp
import os
from file_tools import mkcdir
from time import time
from argument_tools import parse_arguments
import sys
import glob
from bvec_handler import orient_to_str


mainpath = "/Volumes/Data/Badea/ADdecode.01/Analysis/"
mainpath="/mnt/paros_MRI/jacques/AD_Decode/Analysis/"
diff_preprocessed = os.path.join(mainpath, "DWI")

mkcdir([mainpath, diff_preprocessed])

"""
subjects = ["S02666","S02670","S02686","S02654", "S02686", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781",
            "S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02938", "S02939", "S02967", "S02987", "S02987",
            "S03010", "S03033", "S03034", "S03045", "S02524","S02535","S02690","S02715","S02771","S02804","S02812",
            "S02817", "S02840","S02871","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028",
            "S03048", "S03069"]
"""
"""
subjfolder = glob.glob(os.path.join(datapath, "*" + identifier + "*"))[0]
subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
for bxhfile in subjbxh:
    bxhtype = checkbxh(bxhfile, False)
    if bxhtype == "dwi":
        dwipath = bxhfile.replace(".bxh", ".nii.gz")
        break
"""
subjects = ["S02666","S02670","S02686","S02654", "S02686", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781",
            "S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02938", "S02939", "S02967", "S02987", "S02987",
            "S03010", "S03033", "S03034", "S03045", "S02524","S02535","S02690","S02715","S02771","S02804","S02812",
            "S02817", "S02840","S02871","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028",
            "S03048", "S03069"]


subjects_all = glob.glob(os.path.join(diff_preprocessed,'*subjspace_dwi*.nii.gz'))
subjects = []
for subject in subjects_all:
    subject_name = os.path.basename(subject)
    subjects.append(subject_name[:6])
print(subjects)

subjects = ["S02666","S02670","S02686","S02654", "S02686", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781",
            "S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02938", "S02939", "S02967", "S02987", "S02987",
            "S03010", "S03033", "S03034", "S03045", "S02524","S02535","S02690","S02715","S02771","S02804","S02812",
            "S02817", "S02840","S02871","S02877","S02898","S02926","S02938","S02939","S02954", "S03017", "S03028",
            "S03048", "S03069","S03225", "S03293", "S03308", "S02842"]

removed_list = ['S02771',"S03343", "S03350", "S03378", "S03391", "S03394","S03225", "S03293", "S03308", "S02842", "S02804", "S02771"]

#removed_list = []
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

subjects.reverse()
subject_processes, function_processes = parse_arguments(sys.argv,subjects)

#mask types => ['FA', 'T1', 'subjspace']
masktype = "subjspace"
stepsize = 2
overwrite = True
get_params = False
forcestart = False
picklesave = True
verbose = True
get_params = None
doprune = True
bvec_orient = [1,2,3]
vol_b0 = [0,1,2]
classifier = "binary"
symmetric = False
inclusive = False
denoise = "coreg"
savefa = True

make_tracts = True
make_connectomes = False

#classifier types => ["FA", "binary"]
classifiertype = "binary"
brainmask = "subjspace"
labeltype='lrordered'
ratio = 1


if ratio == 1:
    saved_streamlines = "_all"
    trk_folder_name = ""
else:
    saved_streamlines = "_ratio_" + str(ratio)
    trk_folder_name = "_" + str(ratio)

trkpath = os.path.join(mainpath, "TRK_MPCA_fixed")
trkpath = os.path.join(mainpath, "TRK_MPCA_100")
trkpath = os.path.join(mainpath, "TRK_MPCA_fixed"+trk_folder_name)
#trkpath = os.path.join(mainpath, "TRK_MPCA_neworient"+trk_folder_name)
mkcdir(trkpath)


trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
str_identifier = '_stepsize_' + str(stepsize).replace(".","_") + saved_streamlines + roistring

duration1=time()

if forcestart:
    print("WARNING: FORCESTART EMPLOYED. THIS WILL COPY OVER PREVIOUS DATA")



labelslist = []
dwi_results = []
donelist = []
notdonelist = []

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"


#atlas_legends = None
#atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
atlas_legends = mainpath + "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

figspath = os.path.join(mainpath, "Figures_testzone","Figures_MPCA"+inclusive_str+symmetric_str+saved_streamlines)
mkcdir(figspath)

if make_connectomes:
    for subject in subjects:
        picklepath_connect = figspath + subject + str_identifier + '_connectomes.p'
        excel_path = figspath + subject + str_identifier + "_connectomes.xlsx"
        if os.path.exists(picklepath_connect) and os.path.exists(excel_path):
            print("The writing of pickle and excel of " + str(subject) + " is already done")
            donelist.append(subject)
        else:
            notdonelist.append(subject)

dwi_results = []
tract_results = []

print(f'Overwrite is {overwrite}')

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)
    if make_tracts:
        tract_results = pool.starmap_async(create_tracts, [(diff_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                          verbose) for subject
                                                       in subjects]).get()
    if make_connectomes:
        tract_results = pool.starmap_async(tract_connectome_analysis, [(diff_preprocessed, trkpath, str_identifier, figspath,
                                                                       subject, atlas_legends, bvec_orient, inclusive,
                                                                       function_processes, forcestart, picklesave, labeltype, symmetric, overwrite, verbose)
                                                                     for subject in subjects]).get()
    pool.close()
else:
    for subject in subjects:
        if make_tracts:
            tract_results.append(
            create_tracts(diff_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                          verbose))
        #get_diffusionattributes(diff_preprocessed, diff_preprocessed, subject, str_identifier, vol_b0, ratio, bvec_orient,
        #                        masktype, overwrite, verbose)
        if make_connectomes:
            tract_results.append(tract_connectome_analysis(diff_preprocessed, trkpath, str_identifier, figspath, subject,
                                                           atlas_legends, bvec_orient,  brainmask, inclusive,
                                                           function_processes, forcestart, picklesave, labeltype, symmetric, overwrite, verbose))
    print(tract_results)

# dwi_results.append(diff_preprocessing(datapath, diff_preprocessed, subject, bvec_orient, denoise, savefa,
#                                     function_processes, masktype, vol_b0, verbose)) ##Unnecessary, replaced by SAMBA_prep
#dwi_results = pool.starmap_async(diff_preprocessing, [(datapath, diff_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
#                                 masktype, vol_b0, verbose) for subject in subjects]).get()
