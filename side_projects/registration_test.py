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

"""
fetch_stanford_hardi()
nib_stanford, gtab_stanford = read_stanford_hardi()
nifti_data = np.squeeze(nib_stanford.get_data())[..., 0]
nifti_affine = nib_stanford.affine
"""

from dipy.io.image import load_nifti
mynifti = load_nifti("/Volumes/Data/Badea/Lab/19abb14/N57437_nii4D.nii")
nifti_data = np.squeeze(mynifti[0])[..., 0]
nifti_affine = mynifti[1]

mynifti = load_nifti("/Volumes/Data/Badea/Lab/19abb14/N57500_nii4D.nii")
nifti_data2 = np.squeeze(mynifti[0])[..., 0]
nifti_affine2 = mynifti[1]

identity = np.eye(4)
affine_map = AffineMap(identity,
                       nifti_data.shape, nifti_affine,
                       nifti_data2.shape, nifti_affine2)
resampled = affine_map.transform(nifti_data2)
regtools.overlay_slices(nifti_data, resampled, None, 0,
                        "nifti_data", "nifti_data2", "resampled_0.png")
"""
#regtools.overlay_slices(nifti_data, nifti_data2, None, 0,
#                        "nifti_data", "nifti_data2", "resampled_0.png")
#regtools.overlay_slices(nifti_data, resampled, None, 1,
#                        "nifti_data", "nifti_data2", "resampled_1.png")
#regtools.overlay_slices(nifti_data, resampled, None, 2,
                        "nifti_data", "nifti_data2", "resampled_2.png")
"""

c_of_mass = transform_centers_of_mass(nifti_data, nifti_affine,
                                      nifti_data2, nifti_affine2)

transformed = c_of_mass.transform(nifti_data2)
regtools.overlay_slices(nifti_data, transformed, None, 0,
                        "nifti_data", "Transformed", "transformed_com_0.png")
"""
regtools.overlay_slices(nifti_data, transformed, None, 1,
                        "nifti_data", "Transformed", "transformed_com_1.png")
regtools.overlay_slices(nifti_data, transformed, None, 2,
                        "nifti_data", "Transformed", "transformed_com_2.png")
"""
nbins = 32
sampling_prop = None
metric = MutualInformationMetric(nbins, sampling_prop)

level_iters = [10000, 1000, 100]
sigmas = [3.0, 1.0, 0.0]
factors = [4, 2, 1]
affreg = AffineRegistration(metric=metric,
                            level_iters=level_iters,
                            sigmas=sigmas,
                            factors=factors)

transform = TranslationTransform3D()
params0 = None
starting_affine = c_of_mass.affine
translation = affreg.optimize(nifti_data, nifti_data2, transform, params0,
                              nifti_affine, nifti_affine2,
                              starting_affine=starting_affine)

transformed = translation.transform(nifti_data2)
regtools.overlay_slices(nifti_data, transformed, None, 0,
                        "nifti_data", "Transformed", "transformed_trans_0.png")
"""
regtools.overlay_slices(nifti_data, transformed, None, 1,
                        "nifti_data", "Transformed", "transformed_trans_1.png")
regtools.overlay_slices(nifti_data, transformed, None, 2,
                        "nifti_data", "Transformed", "transformed_trans_2.png")
"""
transform = RigidTransform3D()
params0 = None
starting_affine = translation.affine
rigid = affreg.optimize(nifti_data, nifti_data2, transform, params0,
                        nifti_affine, nifti_affine2,
                        starting_affine=starting_affine)

transformed = rigid.transform(nifti_data2)
regtools.overlay_slices(nifti_data, transformed, None, 0,
                        "nifti_data", "Transformed", "transformed_rigid_0.png")
"""
regtools.overlay_slices(nifti_data, transformed, None, 1,
                        "nifti_data", "Transformed", "transformed_rigid_1.png")
regtools.overlay_slices(nifti_data, transformed, None, 2,
                        "nifti_data", "Transformed", "transformed_rigid_2.png")
"""
transform = AffineTransform3D()
params0 = None
starting_affine = rigid.affine
affine = affreg.optimize(nifti_data, nifti_data2, transform, params0,
                         nifti_affine, nifti_affine2,
                         starting_affine=starting_affine)

transformed = affine.transform(nifti_data2)
regtools.overlay_slices(nifti_data, transformed, None, 0,
                        "nifti_data", "Transformed", "transformed_affine_0.png")
"""
regtools.overlay_slices(nifti_data, transformed, None, 1,
                        "nifti_data", "Transformed", "transformed_affine_1.png")
regtools.overlay_slices(nifti_data, transformed, None, 2,
                        "nifti_data", "Transformed", "transformed_affine_2.png")
"""