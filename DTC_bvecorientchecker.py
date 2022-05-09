
import numpy as np
from tract_manager import create_tracts, diff_preprocessing, tract_connectome_analysis, get_diffusionattributes, get_str_identifier
from Daemonprocess import MyPool
import multiprocessing as mp
import os, sys
from file_tools import mkcdir, getfromfile
from time import time
from argument_tools import parse_arguments
from bvec_handler import orient_to_str
import socket
from computer_nav import get_mainpaths, get_atlas
import random

remote=False
project='AMD'
if remote:
    username, passwd = getfromfile('/Users/jas/samos_connect.rtf')
else:
    username = None
    passwd = None
inpath, outpath, atlas_folder, sftp = get_mainpaths(remote,project = project, username=username,password=passwd)
atlas_legends = get_atlas(atlas_folder, 'IIT')

inpath = '/Volumes/Data/Badea/Lab/jacques/AMD_compare/'
#diff_preprocessed = os.path.join(inpath, "DWI_v2")
diff_preprocessed = os.path.join(inpath, 'DWI_orig')

txtpath = os.path.join(inpath, "Parameters")


if not remote:
    mkcdir([outpath, diff_preprocessed, txtpath])
else:
    mkcdir([outpath, diff_preprocessed, txtpath], sftp)

subjects = ["H28029"]

subjects = sorted(subjects)
print(subjects)
random.shuffle(subjects)

subject_processes, function_processes = parse_arguments(sys.argv,subjects)

#"S02230" "S02490" these subjects are strange, to investigate


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
bvec_orient = [1,2,3]
vol_b0 = [0,1,2]
dwi_results = []
donelist = []
notdonelist = []
createmask = masktype
inclusive = False
denoise = "coreg"
savefa = True

stepsize = 2
symmetric = True

reference_weighting = 'fa'
volume_weighting = True
make_tracts = True
make_connectomes = False

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
trkpath = os.path.join(inpath, "TRK"+trk_folder_name)

#atlas_legends = None
#atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

ratio = 1

if ratio == 1:
    saved_streamlines = "_all"
    trk_folder_name = ""
else:
    saved_streamlines = "_ratio_" + str(ratio)
    trk_folder_name = "_" + str(ratio)

fixed=False
if fixed:
    fixed_str = '_fixed'
else:
    fixed_str = ''

trkpath = os.path.join(inpath, "TRK_MPCA_fixed")
trkpath = os.path.join(inpath, "TRK_MPCA_100")
trkpath = os.path.join(inpath, f"TRK_MPCA_v2{fixed_str}"+trk_folder_name)

trkpath = os.path.join(inpath, f"TRK_MPCA_coregdiffrun{fixed_str}"+trk_folder_name)

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

figspath = os.path.join(outpath,"Figures_MPCA"+inclusive_str+symmetric_str+saved_streamlines)

if not remote:
    mkcdir([figspath, trkpath])
else:
    mkcdir([figspath, trkpath], sftp)

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

import itertools
bvec_orient1 = (np.array(list(itertools.permutations([1, 2, 3]))))
bvec_orient2 = [elm*[-1, 1, 1] for elm in bvec_orient1]
bvec_orient3 = [elm*[1, -1, 1] for elm in bvec_orient1]
bvec_orient4 = [elm*[1, 1, -1] for elm in bvec_orient1]

bvec_orient_list = np.concatenate((bvec_orient4, bvec_orient1, bvec_orient2, bvec_orient3))

get_params = False
print(bvec_orient_list)

print(f'Overwrite is {overwrite}')
for subject in subjects:
    txtfile = os.path.join(txtpath, subject + "_" + "params.txt")
    if get_params:
        with open(txtfile, 'a') as fi:
            fi.write("Parameters for subject %s \n" % subject)
    for bvec_orient in bvec_orient_list:
        tract_results = []
        print(bvec_orient)
        strproperty = orient_to_str(bvec_orient)
        print(f'this is the strproperty {strproperty}')
        tract_results.append(
            create_tracts(diff_preprocessed, trkpath, subject, figspath, stepsize, function_processes, strproperty,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params,
                          denoise,
                          verbose, sftp))
        if get_params:
            with open(txtfile, 'a') as f:
                for item in tract_results:
                    f.write("Subject %s with %s %s %s \n" % (
                    item[0], str(bvec_orient[0]), str(bvec_orient[1]), str(bvec_orient[2])))
                    f.write("Num tracts: %s \n" % item[2][0])
                    f.write("Min tract length: %s \n" % item[2][1])
                    f.write("Max tract length: %s \n" % item[2][2])
                    f.write("Average tract length: %s \n" % item[2][3])
                    f.write("Standard deviancy tract length: %s \n" % item[2][4])