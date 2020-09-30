import numpy as np
from dipy.viz import regtools
from dipy.data import fetch_stanford_hardi, read_stanford_hardi
from dipy.data.fetcher import fetch_syn_data, read_syn_data
from dipy.align.imaffine import (transform_centers_of_mass,
                                 AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)
from dipy.align.transforms import (TranslationTransform3D,
                                   RigidTransform3D,
                                   AffineTransform3D)
from dipy.io.image import load_nifti


subjects = ['N54717','N54718','N54719','N54720']

max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

inputpath = "/Users/alex/bass/testdata/lifetest"
target = basepath + "/DWI/"
outputpath = basepath + "/TRK/"
figspath = basepath + "/Figures/"

registration = ["AffineTransform3D"]
nbins = 32
sampling_prop = None
#metric = MutualInformationMetric(nbins, sampling_prop)
level_iters = [10000, 1000, 100]
sigmas = [3.0, 1.0, 0.0]
factors = [4, 2, 1]
params = registrationparams()

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    tract_results = pool.starmap_async(register, [(inputpath, target, outputpath, figspath, registration, metric, level_iters, labelslist, bvec_orient, verbose) for subject in
                                                           l]).get()
#    tract_results = pool.starmap_async(evaluate_tracts, [(dwipath, outtrkpath, subject, stepsize, saved_streamlines,
#                                                          labelslist, outpathpickle, figspath, function_processes,
#                                                          doprune, display, verbose) for subject in l]).get()
    pool.close()
else:
    for subject in subj_list:
        dwi_results.append(dwi_preprocessing(dwipath, outpath, subject, bvec_orient, denoise, savefa, function_processes,
                          labelslist, str_identifier, verbose=False)
        tract_results.append(create_tracts(dwipath, outtrkpath, subject, stepsize, function_processes, strproperty,
                                         saved_streamlines, savefa, labelslist, bvec_orient, verbose))

