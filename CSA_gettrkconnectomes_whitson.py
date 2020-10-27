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

from tract_manager import tract_connectome_analysis
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

l = ['H29056', 'H26578', 'H29060', 'H26637', 'H29264', 'H26765', 'H29225', 'H26660', 'H29304', 'H26890', 'H29556', 'H26862', 'H29410', 'H26966', 'H29403', 'H26841', 'H21593', 'H27126', 'H29618', 'H27111', 'H29627', 'H27164', 'H29502', 'H27100', 'H27381', 'H21836', 'H27391', 'H21850', 'H27495', 'H21729', 'H27488', 'H21915', 'H27682', 'H21956', 'H27686', 'H22331', 'H28208', 'H21990', 'H28955', 'H29878', 'H27719', 'H22102', 'H27841', 'H22101', 'H27842', 'H22228', 'H28029', 'H22140', 'H27852', 'H22276', 'H27999', 'H22369', 'H28115', 'H22644', 'H28308', 'H22574', 'H28377', 'H22368', 'H28325', 'H22320', 'H28182', 'H22898', 'H28748', 'H22683', 'H28373', 'H22536', 'H28433', 'H22825', 'H28662', 'H22864', 'H28698', 'H23143', 'H28861', 'H23157', 'H28820', 'H23028', 'H29002', 'H23210', 'H29020', 'H23309', 'H29161', 'H26949', 'H27163', 'H27246', 'H27869', 'H28068', 'H28262', 'H28856', 'H28869', 'H29044', 'H29089', 'H29127', 'H29242', 'H29254', 'H26745', 'H26850', 'H26880', 'H26958', 'H26974', 'H27017', 'H27610', 'H27640', 'H27680', 'H27778', 'H27982', 'H28338', 'H28437', 'H28463', 'H28532', 'H28809', 'H28857', 'H29013', 'H29025']
l = ['H29056']

max_processors = 15

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

# please set the parameter here

BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"
BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
#BIGGUS_DISKUS = "/mnt/munin6/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
dwipath = BIGGUS_DISKUS

outtrkpath = '/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_QA/'

figspath = BIGGUS_DISKUS + "/C57_JS/Figures_RAS_40subj/"
#figspath = "/Users/alex/jacques/connectomes_testing/Fig_RAS_connectomes/"
#figspath = '/Users/alex/bass/testdata/lifetest/'

outpathpickle = "/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whiston_Figs/"

atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"
atlas_legends = "/Volumes/Data/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_lookup.xlsx"

df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()


stepsize = 2
subject_processes = np.size(l)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)



targetrois = ["Cerebellum"]
ratio = 1
saved_streamlines = "all"
overwrite = "yes"
savefa="only"
verbose=True
denoise='mpca'
savedenoise=True
display=False
savefig=False
doprune=True
#strproperty = "_pypxmz_wholebrain"
allsave=True

trkroi = ["wholebrain"]
if len(trkroi)==1:
    roistring = "_" + trkroi[0] + "_"
elif len(trkroi)>1:
    roistring="_"
    for roi in trkroi:
        roistring = roistring + roi[0:4]
    roistring = roistring + "_"
str_identifier = roistring + "_" + saved_streamlines + '__stepsize_' + str(stepsize)

if len(targetrois)==1:
    targetroistring = "_" + targetrois[0] + "_"
elif len(targetrois)>1:
    roistring="_"
    for roi in rois:
        targetoistring = roistring + roi[0:4]
    targetroistring = roistring + "_"
#labellist=[163,1163,120,1120]
#labelslist=[121,1121]#corpuscallosum

labelslist=[]
for roi in targetrois:
    rslt_df = df.loc[df['Structure'] == roi.lower()]
    if roi.lower() == "wholebrain" or roi.lower() == "brain":
        labelslist=None
    else:
        labelslist=np.concatenate((labelslist,np.array(rslt_df.index)))
print(labelslist)
if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
    txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
    print(txt)

#whitematter_label = []
#rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]
#whitematter_label = np.concatenate((whitematter_label,np.array(rslt_whitem.index2)))

#labelslist=None
bvec_orient = [1,2,-3]
# ---------------------------------------------------------
tall = time()
tract_results=[]


if verbose:
    txt=("Process running with % d max processes available on % d subjects with % d subjects in parallel each using % d processes"
      % (mp.cpu_count(), np.size(l), subject_processes, function_processes))
    print(txt)
    send_mail(txt,subject="Main process start msg ")

duration1=time()


if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)
    tract_results = pool.starmap_async(tract_connectome_analysis, [(dwipath, outtrkpath, str_identifier, figspath,
                                                                    subject, atlas_legends, bvec_orient, verbose)
                                                                   for subject in l]).get()
else:
    for subject in l:
        tract_results.append(tract_connectome_analysis(dwipath, outtrkpath, str_identifier, figspath, subject,
                                                       atlas_legends, bvec_orient, verbose))
