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


#l = ['N57433']
l = ['N54717']
l = ['H28029']

max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"
BIGGUS_DISKUS = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics/"
#dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_RAS/"
dwipath = BIGGUS_DISKUS
#outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'

#outtrkpath = '/Users/alex/bass/testdata/lifetest/'
#outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS/"
outtrkpath = '/Users/alex/bass/testdata/' + 'btable_sanitycheck/'
outtrkpath = '/Volumes/Data/Badea/Lab/mouse/C57_JS/VBM_whistson_QA'
#figspath = BIGGUS_DISKUS + "/C57_JS/Figures_RAS/"
figspath = outtrkpath

#outpathpickle = BIGGUS_DISKUS + "/C57_JS/PicklesFig_RAS/"

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

saved_streamlines = "large"
savefa="no"
verbose=True
denoise='mpca'
savedenoise=True
display=False
savefig=False
doprune=True
get_params=True
strproperty = "_pzmypx_fimbria"
labelslist = [120,1120]#fimbria
bvec_orient = [1, 2, 3]
ratio = 100
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

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(create_tracts, [(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                                            ratio, savefa, labelslist, bvec_orient, get_params, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                         figspath, function_processes, doprune, display, verbose)
#                                                        for subject in l]).get()
    pool.close()
else:
    for subject in l:
        for bvec_orient in bvec_orient_list:
            strproperty = orient_to_str(bvec_orient)
            tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                              ratio, savefa, labelslist, bvec_orient, get_params, verbose))

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
