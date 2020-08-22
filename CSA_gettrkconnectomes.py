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

from tract_manager import create_tracts, evaluate_tracts, dwi_preprocessing, analysis_diffusion_figures, tract_connectome_analysis
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

l = ['N57433','N57434','N57435','N57436','N57437','N57440']
l = ['N57433']
l = ['N57437', 'N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
l = ['N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
#l = ['N57500','N57502','N57504','N57513','N57515','N57518', 'N57520','N57522','N57546','N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
#l = ['N57442', 'N57446', 'N57447','N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57520','N57522','N57546','N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
#l = ['N57520','N57522','N57546', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
#l = ['N57449','N57451','N57496','N57498','N57500','N57502','N57504','N57513','N57515','N57518','N57522','N57546','N57548']
#l = ['N57434','N57435','N57436','N57437','N57440']
#l = ['N57433']
l = ['N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']

max_processors = 5

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'
#dwipath = '/Users/alex/brain_data/19abb14/C57_RAS/'
#BIGGUS_DISKUS = os.environ.get('BIGGUS_DISKUS')
#BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse"
BIGGUS_DISKUS = "/Volumes/dusom_dibs_ad_decode/all_staff/munin3badea/mouse"
#BIGGUS_DISKUS = "/mnt/BIAC/munin3.dhe.duke.edu/Badea/Lab/mouse/"
dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_RAS_40subj/"
#dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_RAS"
#dwipath = "/Users/alex/jacques/connectomes_testing/DWI_RAS/"

#dwipath = "/Users/alex/jacques/tempDWI/"
#outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'

#outtrkpath = '/Users/alex/bass/testdata/lifetest/'
outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS_40subj/"
#outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS"
#outtrkpath = "/Users/alex/jacques/connectomes_testing/TRK_RAS/"

figspath = BIGGUS_DISKUS + "/C57_JS/Figures_RAS_40subj/"
#figspath = "/Users/alex/jacques/connectomes_testing/Fig_RAS_connectomes/"
#figspath = '/Users/alex/bass/testdata/lifetest/'

outpathpickle = BIGGUS_DISKUS + "/C57_JS/PicklesFig_RAS/"

atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"
atlas_legends = "/Users/alex/jacques/connectomes_testing//atlases/CHASSSYMM3AtlasLegends.xlsx"

df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()


stepsize = 2
subject_processes = np.size(l)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

rois = ["fimbria"]
ratio = 10
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
if len(rois)==1:
    strproperty = "_" + rois[0] + "_"
elif len(rois)>1:
    strproperty="_"
    for roi in rois:
        strproperty = strproperty + roi[0:4]
    strproperty = strproperty + "_"
#labellist=[163,1163,120,1120]
#labelslist=[121,1121]#corpuscallosum

labelslist=[]
for roi in rois:
    rslt_df = df.loc[df['Structure'] == roi.lower()]
    if roi.lower() == "wholebrain" or roi.lower() == "brain":
        labelslist=None
    else:
        labelslist=np.concatenate((labelslist,np.array(rslt_df.index2)))
print(labelslist)
if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
    txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
    print(txt)

whitematter_label = []
rslt_whitem = df.loc[df['Subdivisions_7'] == "7_whitematter"]
whitematter_label = np.concatenate((whitematter_label,np.array(rslt_whitem.index2)))

#labelslist=None
bvec_orient=[-2,1,3]
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
    tract_results = pool.starmap_async(tract_connectome_analysis, [(dwipath, outtrkpath, saved_streamlines, strproperty, stepsize, figspath, subject, whitematter_label, rois, labelslist, atlas_legends, bvec_orient, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          allsave, display, strproperty, ratio,
#                                                          verbose) for subject in l]).get()
#    pool.close()
else:
    for subject in l:
        tract_results.append(tract_connectome_analysis(dwipath, outtrkpath, saved_streamlines, strproperty, stepsize, figspath, subject, whitematter_label, rois, labelslist, atlas_legends, bvec_orient, verbose))
#        tract_results.append(evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, labelslist,
#                                             outpathpickle, figspath, function_processes, allsave, display, strproperty,
#                                             ratio, verbose))
