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

from tracking_func import create_tracts, MyPool

l = ['N54717','N54718','N54719','N54720','N54722','N54759','N54760','N54761','N54762','N54763','N54764','N54765',
 'N54766','N54770','N54771','N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806','N54807','N54818',
 'N54824','N54825','N54826','N54837','N54838','N54843','N54844','N54856','N54857','N54858','N54859','N54860','N54861',
 'N54873','N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892','N54893','N54897','N54898','N54899',
 'N54900','N54915','N54916','N54917']
# l = ['N54717','N54718']
print("Running on ", mp.cpu_count(), " processors")
pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'  wenlin make this change
mypath = '/Users/alex/code/Wenlin/data/wenlin_data/'

outpath = '/Users/alex/bass/testdata/' + 'results/'  # wenlin make this change

step_size = 2
peak_processes = 1
saved_streamlines = "small"
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

subject_processes = np.int(mp.cpu_count()/peak_processes)
if peak_processes>1:
    pool = MyPool(subject_processes)
else:
    pool = mp.Pool(subject_processes)

if subject_processes<np.size(l):
    subject_processes = np.size(l)

# ---------------------------------------------------------
tall = time()

results = pool.starmap_async(create_tracts,[(mypath, outpath, subject, step_size, peak_processes,
                                           saved_streamlines, "yes", True) for subject in l]).get()

pool.close()
picklepath = '/Users/alex/jacques/allsubjects_test.p'
pickle.dump(results, open(picklepath,"wb"))

# for j in range(np.size(l)):
#    print(j+1)
#    subject = l[j]
    #pool.starmap_async(create_tracts(mypath,outpath,subject,step_size,peakprocesses))



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