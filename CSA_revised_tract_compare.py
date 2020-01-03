#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eleftherios and Serge

Wenlin make some changes to track on the whole brain
Wenlin add for loop to run all the animals 2018-20-25
"""

from time import time
import numpy as np

import multiprocessing as mp
import pickle

from tracking_func import dwi_create_tracts, evaluate_tracts, extractbvec_fromheader, dwi_preprocessing, MyPool

l = ['N57433', 'N57434', 'N57435', 'N57436','N57437']

#l = ['N54859', 'N54860']
print("Running on ", mp.cpu_count(), " processors")
#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'
dwipath = '/Users/alex/brain_data/19abb14/4DNifti'

outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'
figspath = '/Users/alex/bass/braindata_results/figures/'

stepsize = 2
function_processes = 4
saved_streamlines = "small"
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

subject_processes = np.int(mp.cpu_count()/function_processes)

"""
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57433/co_reg_N57433_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57433',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57434/co_reg_N57434_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57434',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57435/co_reg_N57435_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57435',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57436/co_reg_N57436_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57436',"all")
extractbvec_fromheader('/Users/alex/brain_data/19abb14/N57437/co_reg_N57437_m00.headfile','/Users/alex/brain_data/19abb14/4DNifti/N57437',"all")
"""

savefa="no"
verbose=True
denoise='mpca'
savedenoise=True
display=False
savefig=False
doprune=True
# ---------------------------------------------------------
tall = time()
tract_results=[]

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    #tract_results = pool.starmap_async(dwi_create_tracts, [(dwipath, outtrkpath, subject, stepsize, function_processes,
    #                                                        saved_streamlines, denoise, savefa, verbose) for subject in
    #                                                       l]).get()
    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
                                                         figspath, function_processes, doprune, display, verbose)
                                                        for subject in l]).get()
    pool.close()
else:
    for subject in l:
        tract_results.append(dwi_create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes,
                                          saved_streamlines, denoise, savefa, verbose))


#dwip_results = pool.starmap_async(dwi_preprocessing[(dwipath,outpath,subject,denoise,savefa,function_processes, verbose) for subject in l]).get()

#tract_results = pool.starmap_async(create_tracts,[(dwipath, outpath, subject, stepsize, function_processes,
#                                            saved_streamlines, denoise, savefa, verbose) for subject in l]).get()

subject=l[0]
#dwip_results = dwi_preprocessing(dwipath,dwipath,subject,denoise,savedenoise=savedenoise, savefa=savefa, processes=function_processes, verbose=verbose)
#tract_results = dwi_create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes,
#                                          saved_streamlines, denoise, savefa, verbose)


#tracteval_results = evaluate_tracts(dwipath, outtrkpath, subject, stepsize, saved_streamlines, outpathfig=figspath,
#                                    processes=function_processes, doprune=True, display=display, verbose=verbose)

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