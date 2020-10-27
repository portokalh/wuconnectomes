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
from tract_manager import create_tracts, create_tracts_test
from tract_manager import create_tracts, evaluate_tracts, dwi_preprocessing
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

l = ["H21593", "H21729"]

max_processors = 2

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"
BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
#BIGGUS_DISKUS = "/mnt/munin6/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
dwipath = BIGGUS_DISKUS

#outtrkpath = '/mnt/munin6/Badea/Lab/mouse/C57_JS/VBM_whistson_QA/'
outtrkpath = '/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whistson_QA/'

#outpathpickle = BIGGUS_DISKUS + "/C57_JS/PicklesFig_RAS/"

stepsize = 2
subject_processes = np.size(l)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

rois = ["wholebrain"]
ratio = 10
saved_streamlines = "all"
savefa="no"
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

labelslist = []
if rois and (rois[0]!="wholebrain" or len(rois) > 1):
    atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"
    df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()
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

bvec_orient=[-2,1,3]
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


if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, stepsize, function_processes,
                                                        strproperty, ratio, savefa, labelslist, bvec_orient, doprune,
                                                        overwrite, get_params, verbose) for subject in l]).get()

    pool.close()
else:
    for subject in l:
        tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                              ratio, savefa, labelslist, bvec_orient, doprune, overwrite, get_params,
                                           verbose))


subject=l[0]
#dwip_results = dwi_preprocessing(dwipath,dwipath,subject,denoise,savedenoise=savedenoise, savefa=savefa, processes=function_processes, verbose=verbose)
#tract_results = dwi_create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes,
#                                          saved_streamlines, denoise, savefa, verbose)


#tracteval_results = evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, outpathfig=figspath,
#                                    processes=function_processes, doprune=True, display=display, verbose=verbose)

#        tract_results.append(evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, labelslist,
#                                             outpathpickle, figspath, function_processes, allsave, display, strproperty,
#                                             ratio, verbose))

#dwip_results = pool.starmap_async(dwi_preprocessing[(dwipath,outpath,subject,denoise,savefa,function_processes, verbose) for subject in l]).get()

#tract_results = pool.starmap_async(create_tracts,[(dwipath, outpath, subject, stepsize, function_processes,
#                                            saved_streamlines, denoise, savefa, verbose) for subject in l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          allsave, display, strproperty, ratio, verbose) for subject in l]).get()
"""
tall = time() - tall
if verbose:
    text = ("Process was finished after %.3f s" % (tall))
    print(text)
    send_mail(text, subject="End Process information")

#picklepath = '/Users/alex/jacques/allsubjects_test_eval.p'
#pickle.dump(tracteval_results, open(picklepath,"wb"))

# for j in range(np.size(l)):
#    print(j+1)
#    subject = l[j]
    #pool.starmap_async(create_tracts(mypath,outpath,subject,step_size,function_processes))



duration_all = time() -  tall
print('All animals tracking finished, running time is {}'.format(duration_all))


"""
#example to use parallel in parallel
"""def sleepawhile(t):
    print("Sleeping %i seconds..." % t)
    time.sleep(t)
    return t

def work(num_procs):
    print("Creating %i (daemon) workers and jobs in child." % num_procs)
    pool = multiprocessing.Pool(num_procs)

    result = pool.map(sleepawhile,
        [randint(1, 5) for x in range(num_procs)])

    # The following is not really needed, since the (daemon) workers of the
    # child's pool are killed when the child is terminated, but it's good
    # practice to cleanup after ourselves anyway.
    pool.close()
    pool.join()
    return result

def test():
    print("Creating 5 (non-daemon) workers and jobs in main process.")
    pool = MyPool(5)

    result = pool.map(work, [randint(1, 5) for x in range(5)])

    pool.close()
    pool.join()
    print(result)
"""
