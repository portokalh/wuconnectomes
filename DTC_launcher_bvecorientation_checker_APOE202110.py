#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eleftherios and Serge

Wenlin make some changes to track on the whole brain
Wenlin add for loop to run all the animals 2018-20-25
"""

from time import time
import numpy as np
import os

import multiprocessing as mp
import pickle

from tract_manager import create_tracts
from BIAC_tools import send_mail
from Daemonprocess import MyPool
from argument_tools import parse_arguments
from file_tools import mkcdir
import sys

def orient_to_str(bvec_orient):
    mystr=""
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

l = ['N57500']
l = ['N58651']

subject_processes, function_processes = parse_arguments(sys.argv, l)

#main_folder = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/"
#dwipath = main_folder + "/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"

main_folder = "/Users/alex/jacques/APOE_temp/"
atlas_legends = "/Volumes/Data/Badea/Lab/atlases/CHASSSYMM3AtlasLegends.xlsx"
samos = True
if samos:
    main_folder = "/mnt/paros_MRI/jacques/APOE/"
    atlas_legends = "/mnt/paros_MRI/jacques/atlases/CHASSSYMM3AtlasLegends.xlsx"

dwipath_preprocessed = os.path.join(main_folder,"diff_whiston_preprocessed")
dwipath = os.path.join(main_folder,'DWI_allsubj')
outtrkpath = os.path.join(main_folder,'TRK_bvecs')
figspath = os.path.join(main_folder,"Figures_RAS_allsubj_lr")
txtpath = os.path.join(main_folder, "Parameters")

mkcdir([outtrkpath,figspath,txtpath])

outpathpickle = figspath

stepsize = 2

# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

targetrois = ["Cerebellum"]
ratio = 1
if ratio == 1:
    saved_streamlines = "_all_"
else:
    saved_streamlines = "_ratio_" + str(ratio) + '_'

savefa="yes"
verbose=True
denoise='none'
denoise='coreg'
savedenoise=True
display=False
savefig=False
doprune = True
inclusive=False
allsave=True
labelslist=None

brainmask = "dwi"
classifiertype = "FA"
classifiertype = "binary"
brainmask = "dwi"
labeltype = "lrordered"

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"

trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] #+ "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring #+ "_"
#str_identifier = '_stepsize_' + str(stepsize) + saved_streamlines+ roistring
str_identifier = '_stepsize_' + str(stepsize) + classifiertype + roistring + saved_streamlines #to be reimplemented if full calc, disabled for now
str_identifier = roistring + saved_streamlines + 'stepsize_' + str(stepsize)

tall = time()
tract_results = []


if verbose:
    txt=("Process running with % d max processes available on % d subjects with % d subjects in parallel each using % d processes"
      % (mp.cpu_count(), np.size(l), subject_processes, function_processes))
    print(txt)
    send_mail(txt,subject="Main process start msg ")

duration1=time()
overwrite = False
get_params = False

if overwrite:
    print("WARNING: FORCESTART EMPLOYED. THIS WILL COPY OVER PREVIOUS DATA")
picklesave = True

donelist = []
notdonelist = []
for subject in l:
    picklepath_connect = figspath + subject + str_identifier + '_connectomes.p'
    excel_path = figspath + subject + str_identifier + "_connectomes.xlsx"
    if os.path.exists(picklepath_connect) and os.path.exists(excel_path):
        print("The writing of pickle and excel of " + str(subject) + " is already done")
        donelist.append(subject)
    else:
        notdonelist.append(subject)

#str_identifier='_wholebrain_small_stepsize_2'
createmask = True

dwi_results = []
vol_b0 = [0,1,2,3]

labeltype = 'lrordered'


labelslist = []
bvec_orient = [1, 2, 3]
# ---------------------------------------------------------
tall = time()
tract_results = []

import itertools
bvec_orient1 = (np.array(list(itertools.permutations([1, 2, 3]))))
bvec_orient2 = [elm*[-1, 1, 1] for elm in bvec_orient1]
bvec_orient3 = [elm*[1, -1, 1] for elm in bvec_orient1]
bvec_orient4 = [elm*[1, 1, -1] for elm in bvec_orient1]

bvec_orient_list = np.concatenate((bvec_orient1, bvec_orient2, bvec_orient3, bvec_orient4))

if verbose:
    txt = ("Process running with % d max processes available on % d subjects with % d subjects in parallel each using % d processes"
      % (mp.cpu_count(), np.size(l), subject_processes, function_processes))
    print(txt)
    send_mail(txt, subject="Main process start msg ")

duration1 = time()
txtfile = "/Users/alex/bass/testdata/"

get_params = True
print(bvec_orient_list)

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, figspath, stepsize, function_processes,
                                                        str_identifier, ratio, brainmask, classifiertype, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, denoise, verbose) for subject in l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                         figspath, function_processes, doprune, display, verbose)
#                                                        for subject in l]).get()
    pool.close()
else:
    for subject in l:
        txtfile = os.path.join(txtpath, subject+"_" + "params.txt")
        with open(txtfile, 'a') as fi:
            fi.write("Parameters for subject %s \n" % subject)
        for bvec_orient in bvec_orient_list:
            tract_results = []
            print(bvec_orient)
            strproperty = orient_to_str(bvec_orient)
            print(f'this is the strproperty {strproperty}')
            tract_results.append(create_tracts(dwipath, outtrkpath, subject, figspath, stepsize, function_processes, strproperty,
                                              ratio, brainmask, classifiertype, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                                           verbose))
            print(tract_results)
            with open(txtfile, 'a') as f:
                for item in tract_results:
                    f.write("Subject %s with %s %s %s \n" % (item[0],str(bvec_orient[0]),str(bvec_orient[1]),str(bvec_orient[2])))
                    f.write("Num tracts: %s \n" % item[2][0])
                    f.write("Min tract length: %s \n" % item[2][1])
                    f.write("Max tract length: %s \n" % item[2][2])
                    f.write("Average tract length: %s \n" % item[2][3])
                    f.write("Standard deviancy tract length: %s \n" % item[2][4])
