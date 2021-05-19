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

from tract_manager import create_tracts, evaluate_tracts, dwi_preprocessing, MyPool
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty

import argparse


parser = argparse.ArgumentParser(epilog="TractEVAL (Compare tracts) version 1.0")
parser.add_argument("--v","--verbose", help="output verbosity", action="store_true")
parser.add_argument("--i", type=str, help="Input file path")
parser.add_argument("--isub","--save",help="Save all transition trk files to specified path")
parser.add_argument("--idwi", type=str,help="Input dwi path")
parser.add_argument("--itrk", type=str, help="Input trk path")
parser.add_argument("--o", type=str,help="Output pickle values path and figures if of unspecified")
parser.add_argument("--of", type=str,help="Output figures path")
parser.add_argument("--vis",help="enable Visualization",action="store_true")
parser.add_argument("--s","--save",help="Save all transition trk files to specified path")

args = parser.parse_args()

if args.v: verbose=True
else: verbose = False
if args.vis: display=True
else: display=False
if not args.idwi:
    if not args.i:
    	#windll.Kernel32.SetConsoleTextAttribute(std_output_hdl, 12)
    	print('ERROR   : Input path not specified')
    	#windll.Kernel32.SetConsoleTextAttribute(std_output_hdl, 7)
    	sys.exit()
    else:
    	dwipath=args.i
else:
    dwipath=args.idwi

if not args.trk:
    if not args.i:
    	#windll.Kernel32.SetConsoleTextAttribute(std_output_hdl, 12)
    	print('ERROR   : Input path not specified')
    	#windll.Kernel32.SetConsoleTextAttribute(std_output_hdl, 7)
    	sys.exit()
    else:
    	trkpath=args.i
else:
    trkpath=args.itrk

if not args.isub:
    if os.path.exists(os.joinpath(trkpath,"subject_list.txt"):
        subject_list=os.joinpath(trkpath,"subject_list.txt")

if not args.o:
    print("No output specified, using trk path for all outputs")
    figspath=trkpath
    picklepath=trkpath
else:
    if not args.of:
        figspath=args.o
    else:
        figspath=args.of

    


l = ['N57433','N57434','N57435','N57436','N57437','N57440']
#l = ['N57433']
#l = ['N54717']

max_processors = 100

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'
#dwipath = '/Users/alex/brain_data/19abb14/C57_RAS/'
#BIGGUS_DISKUS = os.environ.get('BIGGUS_DISKUS')
BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"
BIGGUS_DISKUS = "/mnt/BIAC/munin3.dhe.duke.edu/Badea/Lab/mouse/"
dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_RAS/"
#outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'

#outtrkpath = '/Users/alex/bass/testdata/lifetest/'
outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS/"
figspath = BIGGUS_DISKUS + "/C57_JS/Figures_RAS/"
#figspath = '/Users/alex/bass/testdata/lifetest/'

outpathpickle = BIGGUS_DISKUS + "/C57_JS/PicklesFig_RAS/"

atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"

df = pd.read_excel(atlas_legends, sheet_name='Sheet1')
df['Structure'] = df['Structure'].str.lower()


stepsize = 2
subject_processes = np.size(l)
subject_processes = 1
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

rois = ["hypothalamus","septum"]
ratio = 1
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
labelslist=[]#fimbria

for roi in rois:
    rslt_df = df.loc[df['Structure'] == roi.lower()]
    if roi.lower() == "wholebrain" or roi.lower() == "brain":
        labelslist=None
    else:
        labelslist=np.concatenate((labelslist,np.array(rslt_df.index2)))

if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
    txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
    print(txt)

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

#    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
#                                                            saved_streamlines, savefa, labelslist, bvec_orient, verbose) for subject in
#                                                           l]).get()
    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
                                                          labelslist, outpathpickle, figspath, function_processes,
                                                          allsave, display, strproperty, ratio, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in l:
#        tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
#                                         saved_streamlines, savefa, labelslist, bvec_orient, verbose))
        tract_results.append(evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, labelslist,
                                             outpathpickle, figspath, function_processes, allsave, display, strproperty,
                                             ratio, verbose))

#dwip_results = pool.starmap_async(dwi_preprocessing[(dwipath,outpath,subject,denoise,savefa,function_processes, verbose) for subject in l]).get()

#tract_results = pool.starmap_async(create_tracts,[(dwipath, outpath, subject, stepsize, function_processes,
#                                            saved_streamlines, denoise, savefa, verbose) for subject in l]).get()

subject=l[0]
#dwip_results = dwi_preprocessing(dwipath,dwipath,subject,denoise,savedenoise=savedenoise, savefa=savefa, processes=function_processes, verbose=verbose)
#tract_results = dwi_create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes,
#                                          saved_streamlines, denoise, savefa, verbose)


#tracteval_results = evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, outpathfig=figspath,
#                                    processes=function_processes, doprune=True, display=display, verbose=verbose)

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
example to use parallel in parallel
def sleepawhile(t):
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
