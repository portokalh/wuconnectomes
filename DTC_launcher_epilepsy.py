#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eleftherios and Serge

Wenlin make some changes to track on the whole brain
Wenlin add for loop to run all the animals 2018-20-25
"""


import numpy as np
from tract_manager import create_tracts, dwi_preprocessing, tract_connectome_analysis
from Daemonprocess import MyPool
import multiprocessing as mp
import os
from file_tools import mkcdir
from time import time

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

import sys, getopt

subjects = ["00393", "00490", "00560", "00680", "00699", "00795","01952","02263"]
weird_subjects = ["02432"]
atlas_legends = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"

outpath = "/Volumes/Data/Badea/Lab/human/Sinha_epilepsy/"
datapath = os.path.join(outpath, "DWI")
figspath = os.path.join(outpath, "Figures")
dwi_preprocessed = os.path.join(outpath, "DWI_temp")
trkpath = os.path.join(outpath, "TRK")

mkcdir([outpath, figspath, dwi_preprocessed, trkpath])
masktype = "FA"
masktype = "T1"
masktype = "dwi"


max_processors = 10

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

subject_processes = np.size(subjects)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
function_processes = np.int(max_processors / subject_processes)

"""
if masktype == "dwi":
    outpathmask = os.path.join(outpath, subject)
    data, affine, gtab, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, None)
    mask, _ = dwi_to_mask(data, affine, outpathmask, makefig=False, vol_idx=vol_b0, median_radius=5, numpass=6,
                          dilate=2)
elif masktype == "T1":
    #bet bia6_02491_40006.nii.gz 02491.nii.gz -m -o -f 0.4
    #mv 02491_mask.nii.gz 02491_T1_binary_mask.nii.gz
    mask, affinemask = getmask(outpath,subject,"T1",verbose)
"""

stepsize = 2


ratio = 1
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
classifier = "FA"
labelslist = []
bvec_orient = [1,2,-3]
vol_b0 = [0,1,2]

dwi_results = []
donelist = []
notdonelist = []
createmask = masktype
inclusive = True
denoise = None
savefa = True
make_connectomes = False

classifiertype = "FA"
classifiertype = "binary"
brainmask = "dwi"

if classifiertype == "FA":
    classifiertype = "_fa"
else:
    classifiertype = "_binary"


#atlas_legends = None
#atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
atlas_legends = outpath + "/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

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


if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    dwi_results = pool.starmap_async(dwi_preprocessing, [(datapath, dwi_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
                                     createmask, vol_b0, verbose) for subject in subjects]).get()
    tract_results = pool.starmap_async(create_tracts, [(dwi_preprocessed, trkpath, subject, figspath, stepsize, function_processes,
                                                        str_identifier, ratio, masktype, classifier, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, verbose) for subject in subjects]).get()
    if make_connectomes:
        tract_results = pool.starmap_async(tract_connectome_analysis, [(dwi_preprocessed, trkpath, str_identifier, figspath,
                                                                       subject, atlas_legends, bvec_orient, inclusive,
                                                                       function_processes, forcestart, picklesave, verbose)
                                                                     for subject in subjects]).get()
    pool.close()
else:
    for subject in subjects:
        dwi_results.append(dwi_preprocessing(datapath, dwi_preprocessed, subject, bvec_orient, denoise, savefa,
                                             function_processes, createmask, vol_b0, verbose))
        tract_results.append(
            create_tracts(dwi_preprocessed, trkpath, subject, figspath, stepsize, function_processes, str_identifier,
                          ratio, brainmask, classifier, labelslist, bvec_orient, doprune, overwrite, get_params,
                          verbose))
        #get_diffusionattributes(dwi_preprocessed, dwi_preprocessed, subject, str_identifier, vol_b0, ratio, bvec_orient,
        #                        createmask, overwrite, verbose)
        if make_connectomes:
            tract_results.append(tract_connectome_analysis(dwi_preprocessed, trkpath, str_identifier, figspath, subject,
                                                           atlas_legends, bvec_orient,  brainmask, inclusive,
                                                           function_processes, forcestart, picklesave, verbose))
    print(tract_results)
