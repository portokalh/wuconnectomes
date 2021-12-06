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
from tract_manager import create_tracts, tract_connectome_analysis, diff_preprocessing
from BIAC_tools import send_mail
from Daemonprocess import MyPool
import sys, getopt
from file_tools import mkcdir
from argument_tools import parse_arguments

import pickle
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty
from bvec_handler import extractbvec_fromheader

#First set
#l=['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504', 'N57513',
# 'N57515','N57518','N57520','N57522','N57546','N57548','N57550','N57552','N57554','N57559','N57580','N57582','N57584',
# 'N57587','N57590','N57692','N57694','N57700','N57500','N57702','N57709']
#Second set, August update
#l = ["N58214", "N58215", "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
# "N58225", "N58226", "N58228", "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649",
# "N58650", "N58651", "N58653", "N58654", 'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309',
# 'N58792', 'N58302', 'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
# 'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
# 'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
# 'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
# 'N58613', 'N58732', 'N58516', 'N58402']


#'N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504', 'N57513',
# 'N57515','N57518','N57520','N57522','N57546','N57548','N57550','N57552','N57554','N57559','N57580','N57582',
# 'N57584','N57587','N57590','N57692','N57694','N57700','N57702','N57709', 
l = ["N58214", "N58215",
     "N58216", "N58217", "N58218", "N58219", "N58221", "N58222", "N58223", "N58224",
                "N58225", "N58226", "N58228",
                "N58229", "N58230", "N58231", "N58232", "N58633", "N58634", "N58635", "N58636", "N58649", "N58650",
                "N58651", "N58653", "N58654",
                'N58408', 'N58398', 'N58714', 'N58740', 'N58477', 'N58734', 'N58309', 'N58792', 'N58302',
                'N58784', 'N58706', 'N58361', 'N58355', 'N58712', 'N58790', 'N58606', 'N58350', 'N58608',
                'N58779', 'N58500', 'N58604', 'N58749', 'N58510', 'N58394', 'N58346', 'N58344', 'N58788', 'N58305',
                'N58514', 'N58794', 'N58733', 'N58655', 'N58735', 'N58310', 'N58400', 'N58708', 'N58780', 'N58512',
                'N58747', 'N58303', 'N58404', 'N58751', 'N58611', 'N58745', 'N58406', 'N58359', 'N58742', 'N58396',
                'N58613', 'N58732', 'N58516', 'N58402']
l = ["N58634"]
l = ['N57442']

print("Will go from subject "+ l[0] + " to subject "+l[-1])

subject_processes, function_processes = parse_arguments(sys.argv, l)

#main_folder = "/Volumes/dusom_dibs_ad_decode/all_staff/APOE_temp/"
#dwipath = main_folder + "/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"

main_folder = "/Users/alex/jacques/APOE_temp/"
atlas_legends = "/Volumes/Data/Badea/Lab/atlases/CHASSSYMM3AtlasLegends.xlsx"
samos = False
if samos:
    main_folder = "/mnt/paros_MRI/jacques/APOE/"
    atlas_legends = "/mnt/paros_MRI/jacques/atlases/CHASSSYMM3AtlasLegends.xlsx"

dwipath_preprocessed = os.path.join(main_folder,"diff_whiston_preprocessed")
dwipath = os.path.join(main_folder,'DWI_allsubj')
outtrkpath = os.path.join(main_folder + 'TRK_allsubj')
figspath = os.path.join(main_folder + "Figures_RAS_allsubj_lr")

mkcdir([outtrkpath,figspath])

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

if overwrite:x
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
#labelslist is only good if you want to build tracts pertaining to specific regions, not used currently

"""
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
"""
make_connectomes=True
if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)
    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, figspath, stepsize, function_processes,
                                                        str_identifier, ratio, brainmask, classifiertype, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, denoise, verbose) for subject in l]).get()
    if make_connectomes:
        tract_results = pool.starmap_async(tract_connectome_analysis, [(dwipath, outtrkpath, str_identifier, figspath,
                                                                   subject, atlas_legends, bvec_orient, brainmask,
                                                                    inclusive,function_processes, overwrite,
                                                                    picklesave, labeltype, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in l:
       tract_results.append(create_tracts(dwipath, outtrkpath, subject, figspath, stepsize, function_processes, str_identifier,
                                              ratio, brainmask, classifiertype, labelslist, bvec_orient, doprune, overwrite, get_params, denoise,
                                           verbose))
       if make_connectomes:
            tract_results.append(tract_connectome_analysis(dwipath, outtrkpath, str_identifier, figspath, subject,
                                                     atlas_legends, bvec_orient, brainmask, inclusive, function_processes,
                                                     overwrite, picklesave, labeltype, verbose))


subject=l[0]

# dwi_results.append(dwi_preprocessing(dwipath, dwipath_preprocessed, subject, bvec_orient, denoise, savefa,
#                                  function_processes, createmask, vol_b0, verbose))
