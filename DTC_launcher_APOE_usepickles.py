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
from tract_manager import create_tracts, tract_connectome_analysis_pickle, dwi_preprocessing
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

import sys, getopt


l = ['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498']
#l = ['N57446', 'N57447','N57449','N57451','N57496','N57498']
#l=['N57447']
#l=['N57496']
l = ['N57498']
l = ['N57692']

#l = ["H29410", "H29060"]
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hb:e:", ["first=", "last="])
except getopt.GetoptError:
    print('test.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('test.py -b first -s last')
        sys.exit()
    elif opt in ("-b", "--first"):
        start = arg
    elif opt in ("-e", "--last"):
        end = arg

if 'start' in locals():
    del(start, end)
if 'start' in locals():
    start = int(start)
    if 'end' in locals():
        l = l[int(start):int(end)+1]
    else:
        l = l[start:]
if 'start' not in locals():
    if 'end' not in locals():
        l = l
    else:
        l = l[0:end]
print("Will go from subject "+ l[0] + " to subject "+l[-1])
max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

main_folder = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/"
main_folder = "/Users/alex/jacques/APOE_temp/"
dwipath = main_folder + "/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
dwipath_preprocessed = main_folder + "/C57_JS/diff_whiston_preprocessed/"
dwipath = main_folder + '/DWI_allsubj/'
outtrkpath = main_folder + '/TRK_allsubj/'
figspath = main_folder + "/Figures_RAS_40subj_lr/"

outpathpickle = figspath

atlas_legends = main_folder + "/../atlases/IITmean_RPI/IITmean_RPI_lookup.xlsx"
atlas_legends = main_folder + "/../atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
atlas_legends = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"

stepsize = 2
subject_processes = np.size(l)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

targetrois = ["Cerebellum"]
ratio = 100
if ratio == 1:
    saved_streamlines = "_all_"
else:
    saved_streamlines = "_ratio_" + str(ratio)

savefa="yes"
verbose=True
denoise='none'
#denoise=None
savedenoise=True
display=False
savefig=False
doprune = True
inclusive=False
allsave=True

classifiertype = "FA"
classifiertype = "binary"
brainmask = "dwi"

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
str_identifier = '_stepsize_' + str(stepsize) + classifiertype + roistring + saved_streamlines

labelslist=[]
if targetrois and (targetrois[0]!="wholebrain" or len(targetrois) > 1):
    df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()
    for roi in targetrois:
        rslt_df = df.loc[df['Structure'] == roi.lower()]
        if roi.lower() == "wholebrain" or roi.lower() == "brain":
            labelslist=None
        else:
            labelslist=np.concatenate((labelslist, np.array(rslt_df.index)))
print(labelslist)
if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
    txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
    print(txt)

bvec_orient=[1,2,-3]
bvec_orient=[-2,1,3]

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
forcestart = False
if forcestart:
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

str_identifier='_wholebrain_small_stepsize_2'
createmask = True

dwi_results = []
vol_b0 = [0,1,2,3]

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    #dwi_results = pool.starmap_async(dwi_preprocessing, [(dwipath, dwipath_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
    #                                 createmask, vol_b0, verbose) for subject in l]).get()
    #tract_results = pool.starmap_async(create_tracts, [(dwipath_preprocessed, outtrkpath, subject, figspath, stepsize, function_processes,
    #                                                    str_identifier, ratio, classifiertype, labelslist, bvec_orient, doprune,
    #                                                    overwrite, get_params, verbose) for subject in l]).get()
    tract_results = pool.starmap_async(tract_connectome_analysis_pickle, [(dwipath_preprocessed, outtrkpath, str_identifier, figspath,
                                                                   subject, atlas_legends, bvec_orient, brainmask,
                                                                    inclusive,function_processes, forcestart,
                                                                    picklesave, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in l:
       #dwi_results.append(dwi_preprocessing(dwipath, dwipath_preprocessed, subject, bvec_orient, denoise, savefa,
       #                                  function_processes, createmask, vol_b0, verbose))
       #tract_results.append(create_tracts(dwipath_preprocessed, outtrkpath, subject, figspath, stepsize, function_processes, str_identifier,
       #                                       ratio, classifiertype, labelslist, bvec_orient, doprune, overwrite, get_params,
       #                                    verbose))
       tract_results.append(tract_connectome_analysis_pickle(dwipath, outtrkpath, str_identifier, figspath, subject,
                                                     atlas_legends, bvec_orient, brainmask, inclusive, function_processes,
                                                     forcestart, picklesave, verbose))


subject=l[0]