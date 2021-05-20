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

import multiprocessing as mp

from registration_handler import register_save, registrationparams

subjects = ['N57500']

max_processors = 1

if mp.cpu_count() < max_processors:
    max_processors = mp.cpu_count()

print("Running on ", max_processors, " processors")

inputpath = "/Users/alex/jacques/registration_test/"
target = "/Volumes/Data/Badea/Lab/mouse/C57_JS/DWI_RAS_40subj/N57437_nii4D_RAS.nii.gz"
toapply = "TRK"
outputpath = "/Users/alex/jacques/registration_test/registered/"
outputpath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/registration_test/"
figspath = "/Volumes/Data/Badea/Lab/mouse/C57_JS/registration_test/"
figspath = inputpath + "/Figures/"
figspath = None
applydirs = ["/Users/alex/jacques/registration_test/"]

nbins = 32
sampling_prop = None
metric = MutualInformationMetric(nbins, sampling_prop)
level_iters = [10000, 1000, 100]
sigmas = [3.0, 1.0, 0.0]
factors = [4, 2, 1]
params = registrationparams(nbins, sampling_prop, level_iters, sigmas, factors)
verbose = True

subject_processes = np.size(subjects)
subject_processes = 10
if max_processors < subject_processes:
    subject_processes = max_processors
# accepted values are "small" for one in ten streamlines, "all or "large" for all streamlines,
# "none" or None variable for neither and "both" for both of them

function_processes = np.int(max_processors/subject_processes)

register_save_path = []
registration_type = ["center_mass", "AffineRegistration", "RigidTransform3D", "AffineTransform3D"]

if subject_processes>1:
    if function_processes>1:
        pool = MyPool(subject_processes)
    else:
        pool = mp.Pool(subject_processes)

    register_save_path = pool.starmap_async(register_save, [(inputpath, target, subject, outputpath, figspath, params,
                                           registration_type, applydirs, verbose) for subject in
                                                           l]).get()
    """register_apply_path = pool.starmap_async(register_apply, [(inputpath, target, subject, outputpath, figspath,
                                                             registrationparams, registration_type, verbose) for subject
                                                        in l]).get()
    """
    pool.close()
else:
    for subject in subjects:
        register_save_path.append(register_save(inputpath, target, subject, outputpath, figspath, params,
                                           registration_type, applydirs, verbose))
        """
        register_apply_path.append(register_save(inputpath, target, subject, outputpath, figspath, registrationparams,
                                           registration_type, verbose))
        """
