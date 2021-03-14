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
from tract_manager import create_tracts, tract_connectome_analysis, dwi_preprocessing, copylabels
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

import sys, getopt

l = ["H26637", "H26966"] #"H29410", "H29060"
l = ["H26966"]
l = ["H29410"]
#l = ["H29060"]
#l = ["H26637"]
l = ["H26637", "H26966", "H29410", "H29060"]
l = ["H26637"]


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hb:e:",["first=","last="])
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
"""
l = ["H21593", "H21729"]
l = ["H21729"]
l = ['H22102', 'H27841', 'H22101',
 'H27842', 'H22228', 'H28029', 'H22140', 'H27852', 'H22276', 'H27999', 'H22369', 'H28115', 'H22644', 'H28308',
 'H22574', 'H28377', 'H22368', 'H28325', 'H22320', 'H28182', 'H22898', 'H28748', 'H22683', 'H28373', 'H22536',
 'H28433', 'H22825', 'H28662', 'H22864', 'H28698', 'H23143', 'H28861', 'H23157', 'H28820', 'H23028', 'H29002',
 'H23210', 'H29020', 'H23309', 'H29161', 'H26949', 'H27163', 'H27246', 'H27869', 'H28068', 'H28262', 'H28856',
 'H28869', 'H29044', 'H29089', 'H29127', 'H29242', 'H29254', 'H26745', 'H26850', 'H26880', 'H26958', 'H26974',
 'H27017', 'H27610', 'H27640', 'H27680', 'H27778', 'H27982', 'H28338', 'H28437', 'H28463', 'H28532', 'H28809',
 'H28857', 'H29013', 'H29025']

l = ['H29056', 'H26578', 'H29060', 'H26637', 'H29264', 'H26765', 'H29225', 'H26660', 'H29304', 'H26890', 'H29556',
     'H26862', 'H29410', 'H26966', 'H29403', 'H26841', 'H21593', 'H27126', 'H29618', 'H27111', 'H29627', 'H27164',
     'H29502', 'H27100', 'H27381', 'H21836', 'H27391', 'H21850', 'H27495', 'H21729', 'H27488', 'H21915', 'H27682',
     'H21956', 'H27686', 'H22331', 'H28208', 'H21990', 'H28955', 'H29878', 'H27719', 'H22102', 'H27841', 'H22101',
     'H27842', 'H22228', 'H28029', 'H22140', 'H27852', 'H22276', 'H27999', 'H22369', 'H28115', 'H22644', 'H28308',
     'H22574', 'H28377', 'H22368', 'H28325', 'H22320', 'H28182', 'H22898', 'H28748', 'H22683', 'H28373', 'H22536',
     'H28433', 'H22825', 'H28662', 'H22864', 'H28698', 'H23143', 'H28861', 'H23157', 'H28820', 'H23028', 'H29002',
     'H23210', 'H29020', 'H23309', 'H29161', 'H26949', 'H27163', 'H27246', 'H27869', 'H28068', 'H28262', 'H28856',
     'H28869', 'H29044', 'H29089', 'H29127', 'H29242', 'H29254', 'H26745', 'H26850', 'H26880', 'H26958', 'H26974',
     'H27017', 'H27610', 'H27640', 'H27680', 'H27778', 'H27982', 'H28338', 'H28437', 'H28463', 'H28532', 'H28809',
     'H28857', 'H29013', 'H29025']
"""

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

max_processors = 20

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse"
#BIGGUS_DISKUS = "/mnt/munin6/Badea/Lab/mouse"
#BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
#BIGGUS_DISKUS = "/mnt/munin6/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
dwipath = BIGGUS_DISKUS + "/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
dwipath_preprocessed = BIGGUS_DISKUS + "/C57_JS/diff_whiston_preprocessed/"
#outtrkpath = '/mnt/munin6/Badea/Lab/mouse/C57_JS/VBM_whistson_QA/'
outtrkpath = BIGGUS_DISKUS + '/C57_JS/VBM_whiston_QA_new/'
#outtrkpath = "/Volumes/dusom_dibs_ad_decode/all_staff/VBM_whiston_QA/"

figspath = BIGGUS_DISKUS + '/C57_JS/VBM_whiston_Figs_inclusive_new/'

outpathpickle = figspath

atlas_legends = BIGGUS_DISKUS + "/../atlases/IITmean_RPI/IITmean_RPI_lookup.xlsx"
atlas_legends = BIGGUS_DISKUS + "/../atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

stepsize = 2
subject_processes = np.size(l)

if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

targetrois = ["Cerebellum"]
ratio = 10
if ratio == 1:
    saved_streamlines = "_all"
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
inclusive=True
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
str_identifier = '_stepsize_' + str(stepsize) + maskuse + roistring + saved_streamlines

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
# ---------------------------------------------------------
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
forcestart = True
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

createmask = True

dwi_results = []
vol_b0 = [0,1,2,3]

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    dwi_results = pool.starmap_async(dwi_preprocessing, [(dwipath, dwipath_preprocessed, subject, bvec_orient, denoise, savefa, function_processes,
                                     createmask, vol_b0, verbose) for subject in l]).get()
    tract_results = pool.starmap_async(create_tracts, [(dwipath_preprocessed, outtrkpath, subject, figspath, stepsize, function_processes,
                                                        str_identifier, ratio, classifiertype, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, verbose) for subject in l]).get()
    pool.starmap_async = copylabels(dwipath, dwipath_preprocessed, subject, verbose)
    tract_results = pool.starmap_async(tract_connectome_analysis, [(dwipath_preprocessed, outtrkpath, str_identifier, figspath,
                                                                   subject, atlas_legends, bvec_orient, inclusive,
                                                                   function_processes, forcestart, picklesave, verbose)
                                                                 for subject in l]).get()
    pool.close()
else:
    for subject in l:
       #dwi_results.append(dwi_preprocessing(dwipath, dwipath_preprocessed, subject, bvec_orient, denoise, savefa,
       #                                  function_processes, createmask, vol_b0, verbose))
       #tract_results.append(create_tracts(dwipath_preprocessed, outtrkpath, subject, figspath, stepsize, function_processes, str_identifier,
       #                                       ratio, classifiertype, labelslist, bvec_orient, doprune, overwrite, get_params,
       #                                    verbose))
       copylabels(dwipath, dwipath_preprocessed, subject, verbose)
       tract_results.append(tract_connectome_analysis(dwipath_preprocessed, outtrkpath, str_identifier, figspath, subject,
                                                     atlas_legends, bvec_orient, inclusive, function_processes,
                                                     forcestart, picklesave, verbose))


subject=l[0]
