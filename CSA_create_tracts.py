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

l = ['N57433','N57434','N57435','N57436','N57437','N57440']
#l = ['N57433']
#l = ['N54717']

max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'
#dwipath = '/Users/alex/brain_data/19abb14/C57_RAS/'
#BIGGUS_DISKUS = os.environ.get('BIGGUS_DISKUS')
BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"
dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_ARI/"
#outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'

#outtrkpath = '/Users/alex/bass/testdata/lifetest/'
outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_ARI/"
figspath = BIGGUS_DISKUS + "/C57_JS/Figures_ARI/"
#figspath = '/Users/alex/bass/testdata/lifetest/'

outpathpickle = BIGGUS_DISKUS + "/C57_JS/PicklesFig_ARI/"

stepsize = 2
subject_processes = np.size(l)
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

"""
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57433/co_reg_N57433_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57433',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57434/co_reg_N57434_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57434',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57435/co_reg_N57435_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57435',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57436/co_reg_N57436_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57436',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57437/co_reg_N57437_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57437',"all")
"""

saved_streamlines = "all"
savefa="no"
verbose=True
denoise='mpca'
savedenoise=True
display=False
savefig=False
doprune=True
#strproperty = "_pypxmz_wholebrain"
strproperty = "_pxmymz_wholebrain"
#labellist=[163,1163,120,1120]
#labelslist=[121,1121]#corpuscallosum
labelslist=[120,1120]#fimbria
labelslist=None
bvec_orient=[1,-2,-3]
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

    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                                            saved_streamlines, savefa, labelslist, bvec_orient, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          doprune, display, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in l:
        tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                         saved_streamlines, savefa, labelslist, bvec_orient, verbose))
#        tract_results.append(evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, labelslist,
#                                             outpathpickle, figspath, function_processes, doprune, display, verbose))

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
