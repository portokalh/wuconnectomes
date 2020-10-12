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

from tracking_func import dwi_create_tracts, evaluate_tracts, extractbvec_fromheader, dwi_preprocessing, MyPool, send_mail

l = [ 'N57437', 'N57442', 'N57446', 'N57447', 'N57449', 'N57451', 'N57496', 'N57498', 'N57500', 'N57502', 'N57504', 'N57513', 'N57515', 'N57518', 'N57520', 'N57522', 'N57546', 'N57548', 'N57550', 'N57552', 'N57554', 'N57559', 'N57580', 'N57582', 'N57584', 'N57587', 'N57590', 'N57692', 'N57694', 'N57700', 'N57702', 'N57709']
l =  ['N57437']

max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

#pool = mp.Pool(mp.cpu_count())

# please set the parameter here

# mypath = '/Users/alex/brain_data/E3E4/wenlin/'
#dwipath = '/Users/alex/brain_data/19abb14/4DNifti'
#BIGGUS_DISKUS = os.environ.get('BIGGUS_DISKUS')
BIGGUS_DISKUS = "/mnt/BIAC/munin3.dhe.duke.edu/Badea/Lab/mouse"
dwipath = BIGGUS_DISKUS + "/../19abb14/"

#outtrkpath = '/Users/alex/bass/testdata/' + 'braindata_results/'
outtrkpath = BIGGUS_DISKUS + "/../19abb14/"

figspath = BIGGUS_DISKUS + "/../19abb14/"

stepsize = 2
subject_processes = np.size(l)
if max_processors < subject_processes:
    subject_processes = max_processors
saved_streamlines = "small"
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

savefa="yes"
verbose=True
denoise='none' #mpcas
savedenoise=True
display=False
savefig=False
doprune=True
strproperty = "_zflip"
labellist=None
# ---------------------------------------------------------
tall = time()
tract_results=[]


if verbose:
    txt=("Process running with % d max processes available on % d subjects with % d subjects in parallel each using % d processes"
      % (mp.cpu_count(), np.size(l), subject_processes, function_processes))
    print(txt)
    send_mail(txt,subject="LifE start msg ")

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    dwip_results = pool.starmap_async(dwi_preprocessing, [(dwipath,dwipath,subject,denoise, savefa,function_processes,labelslist, strproperty,
                                                         verbose) for subject in l]).get()

    pool.close()
else:
    for subject in l:
        dwip_results = dwi_preprocessing(dwipath,dwipath,subject,denoise, savefa=savefa, processes=function_processes, labelslist=labellist, strproperty=strproperty, verbose=verbose)

duration_all = time() -  tall
print('All animals tracking finished, running time is {}'.format(duration_all))


