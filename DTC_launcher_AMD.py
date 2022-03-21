
import numpy as np
from tract_manager import create_tracts, diff_preprocessing, tract_connectome_analysis, get_diffusionattributes, get_str_identifier
from Daemonprocess import MyPool
import multiprocessing as mp
import os, sys
from file_tools import mkcdir
from time import time
from argument_tools import parse_arguments
from bvec_handler import orient_to_str
import socket

computer_name = socket.gethostname()

if 'samos' in computer_name:
    mainpath = '/mnt/paros_MRI/jacques'
    atlas_legends = mainpath + "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in computer_name or 'hydra' in computer_name:
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    atlas_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    raise Exception('No other computer name yet')

project = "AMD"

if project == 'AD_Decode':
    mainpath = os.path.join(mainpath, 'AD_Decode', 'Analysis')
else:
    mainpath = os.path.join(mainpath, project)

diff_preprocessed = os.path.join(mainpath, "DWI")

mkcdir([mainpath, diff_preprocessed])

subjects = ["H29056", "H26578", "H29060", "H26637", "H29264", "H26765", "H29225", "H26660", "H29304", "H26890",
            "H29556", "H26862", "H29410", "H26966", "H29403", "H26841", "H21593", "H27126", "H29618", "H27111", "H29627",
            "H27164", "H29502", "H27100", "H27381", "H21836", "H27391", "H21850", "H27495", "H21729", "H27488", "H21915",
            "H27682", "H21956", "H27686", "H22331", "H28208", "H21990", "H28955", "H29878", "H27719", "H22102", "H27841",
            "H22101", "H27842", "H22228", "H28029", "H22140", "H27852", "H22276", "H27999", "H22369", "H28115", "H22644",
            "H28308", "H22574", "H28377", "H22368", "H28325", "H22320", "H28182", "H22898", "H28748", "H22683", "H28373",
            "H22536", "H28433", "H22825", "H28662", "H22864", "H28698", "H23143", "H28861", "H23157", "H28820", "H23028",
            "H29002", "H23210", "H29020", "H23309", "H29161", "H26841", "H26862", "H26949", "H26966", "H27100", "H27126",
            "H27163", "H27246", "H27488", "H27682", "H27686", "H27719", "H27841", "H27842", "H27852", "H27869", "H27999",
            "H28029", "H28068", "H28208", "H28262", "H28325", "H28820", "H28856", "H28869", "H28955", "H29002", "H29044",
            "H29089", "H29127", "H29161", "H29242", "H29254", "H26578", "H26637", "H26660", "H26745", "H26765", "H26850",
            "H26880", "H26890", "H26958", "H26974", "H27017", "H27111", "H27164", "H27381", "H27391", "H27495", "H27610",
            "H27640", "H27680", "H27778", "H27982", "H28115", "H28308", "H28338", "H28373", "H28377", "H28433", "H28437",
            "H28463", "H28532", "H28662", "H28698", "H28748", "H28809", "H28857", "H28861", "H29013", "H29020", "H29025"]


subject_processes, function_processes = parse_arguments(sys.argv,subjects)

#"S02230" "S02490" these subjects are strange, to investigate

"""
subjfolder = glob.glob(os.path.join(datapath, "*" + identifier + "*"))[0]
subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
for bxhfile in subjbxh:
    bxhtype = checkbxh(bxhfile, False)
    if bxhtype == "dwi":
        dwipath = bxhfile.replace(".bxh", ".nii.gz")
        break
"""

#outpath = "/Volumes/Data/Badea/ADdecode.01/Analysis/"
#outpath = "/Users/alex/jacques/AMD_TRK_testing/"
diff_preprocessed = os.path.join(mainpath, "DWI")

mkcdir([mainpath, figspath, diff_preprocessed])
masktype = "FA"
masktype = "T1"
masktype = "dwi"

stepsize = 2
ratio = 1
trkroi="wholebrain"

str_identifier = get_str_identifier(stepsize, ratio, trkroi)
if ratio == 1:
    saved_streamlines = "_all"
else:
    saved_streamlines = "_ratio_" + str(ratio)

trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
#str_identifier = '_stepsize_' + str(stepsize) + saved_streamlines+ roistring
str_identifier = '_stepsize_' + str(stepsize).replace(".","_") + saved_streamlines + roistring

duration1=time()
overwrite = False
get_params = False
forcestart = False
if forcestart:
    print("WARNING: FORCESTART EMPLOYED. THIS WILL COPY OVER PREVIOUS DATA")
picklesave = True
verbose = True
get_params = None
doprune = True
#classifier = ["FA", "binary"]
classifier = "binary"
labelslist = []
bvec_orient = [1,2,-3]
vol_b0 = [0,1,2]
dwi_results = []
donelist = []
notdonelist = []
createmask = masktype
inclusive = False
denoise = "coreg"
savefa = True
make_connectomes = False

stepsize = 2
symmetric = True

reference_weighting = 'fa'
volume_weighting = True
make_tracts = False
make_connectomes = True

classifiertype = "FA"
classifiertype = "binary"
brainmask = "dwi"
labeltype='lrordered'

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"

if ratio == 1:
    saved_streamlines = "_all"
    trk_folder_name = ""
else:
    saved_streamlines = "_ratio_" + str(ratio)
    trk_folder_name = "_" + str(ratio)
trkpath = os.path.join(mainpath, "TRK_MPCA_fixed"+trk_folder_name)

#atlas_legends = None
#atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

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

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

figspath = os.path.join(mainpath,"Figures_MPCA"+inclusive_str+symmetric_str+saved_streamlines)

mkcdir(figspath, trkpath)

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
                                                                       subject, atlas_legends, bvec_orient, brainmask, inclusive,
                                                                       function_processes, overwrite, picklesave, labeltype, symmetric, reference_weighting, volume_weighting, verbose)
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
                                                           function_processes, overwrite, picklesave, labeltype, symmetric, reference_weighting, volume_weighting, verbose))
    print(tract_results)