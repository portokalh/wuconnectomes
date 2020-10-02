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

class registrationparams(object):

    def __init__(self, nbins=32, sampling_proportion=None, level_iters=[10000, 1000, 100], ):
        r"""Initialize an instance of the Mutual Information metric.

        This class implements the methods required by Optimizer to drive the
        registration process.

        Parameters
        ----------
        nbins : int, optional
            the number of bins to be used for computing the intensity
            histograms. The default is 32.
        sampling_proportion : None or float in interval (0, 1], optional
            There are two types of sampling: dense and sparse. Dense sampling
            uses all voxels for estimating the (joint and marginal) intensity
            histograms, while sparse sampling uses a subset of them. If
            `sampling_proportion` is None, then dense sampling is
            used. If `sampling_proportion` is a floating point value in (0,1]
            then sparse sampling is used, where `sampling_proportion`
            specifies the proportion of voxels to be used. The default is
            None.

        Notes
        -----
        Since we use linear interpolation, images are not, in general,
        differentiable at exact voxel coordinates, but they are differentiable
        between voxel coordinates. When using sparse sampling, selected voxels
        are slightly moved by adding a small random displacement within one
        voxel to prevent sampling points from being located exactly at voxel
        coordinates. When using dense sampling, this random displacement is
        not applied.

        """
        nbins, sampling_prop, level_iters, sigmas, factors
        self.nbins = nbins
        self.sampling_proportion = sampling_proportion
        self.level_iters = level_iters
        self.sigmas = sigmas
        self.factors = factors

def register_save(inputpathdir, target_path, subject, outputpath, figspath, registrationparams, registration_type, verbose)

    mynifti = get_nifti(inputpathdir, subject)
    mynifti = load_nifti(mynifti)
    #mynifti = load_nifti("/Volumes/Data/Badea/Lab/19abb14/N57437_nii4D.nii")
    nifti_data = np.squeeze(mynifti[0])[..., 0]
    nifti_affine = mynifti[1]

    target = load_nifti(target_path)
    target_data = np.squeeze(mynifti[0])[..., 0]
    target_affine = mynifti[1]

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