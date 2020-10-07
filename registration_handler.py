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

from dipy.io.image import load_nifti, save_nifti
from file_catcher import get_anat


class registrationparams(object):

    def __init__(self, nbins=32, sampling_proportion=None, level_iters=[10000, 1000, 100], sigmas = [3.0, 1.0, 0.0],
                 factors = [4, 2, 1]):
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
        self.nbins = nbins
        self.sampling_proportion = sampling_proportion
        self.level_iters = level_iters
        self.sigmas = sigmas
        self.factors = factors

def register_save(inputpathdir, target_path, subject, outputpath, figspath, registrationparams, registration_types, verbose):

    mynifti = get_anat(inputpathdir, subject)
    mynifti = load_nifti(mynifti)
    #mynifti = load_nifti("/Volumes/Data/Badea/Lab/19abb14/N57437_nii4D.nii")
    anat_data = np.squeeze(mynifti[0])[..., 0]
    anat_affine = mynifti[1]

    target = load_nifti(target_path)
    target_data = np.squeeze(target[0])[..., 0]
    target_affine = target[1]

    identity = np.eye(4)

    affine_map = AffineMap(identity,
                           target_data.shape, target_affine,
                           anat_data.shape, anat_affine)
    resampled = affine_map.transform(anat_data)
    regtools.overlay_slices(target_data, resampled, None, 0,
                            "target_data", "anat_data", "resampled_0.png")
    """
    regtools.overlay_slices(target_data, resampled, None, 1,
                            "target_data", "anat_data", "resampled_1.png")
    regtools.overlay_slices(target_data, resampled, None, 2,
                            "target_data", "anat_data", "resampled_2.png")
    """

    c_of_mass = transform_centers_of_mass(target_data, target_affine,
                                          anat_data, anat_affine)

    if "center_mass" in registration_types:
        transformed = c_of_mass.transform(anat_data)
        """
        regtools.overlay_slices(target_data, transformed, None, 0,
                                "target_data", "Transformed", "transformed_com_0.png")
        regtools.overlay_slices(target_data, transformed, None, 1,
                                "target_data", "Transformed", "transformed_com_1.png")
        regtools.overlay_slices(target_data, transformed, None, 2,
                                "target_data", "Transformed", "transformed_com_2.png")
        """
        fname = outputpath + 
        save_nifti(outpath_mpca, denoised_arr, affine, hdr=hdr)

    metric = MutualInformationMetric(registrationparams.nbins, registrationparams.sampling_prop)

    affreg = AffineRegistration(metric=metric,
                                level_iters=registrationparams.level_iters,
                                sigmas=registrationparams.sigmas,
                                factors=registrationparams.factors)

    transform = TranslationTransform3D()
    params0 = None
    starting_affine = c_of_mass.affine
    translation = affreg.optimize(target_data, anat_data, transform, params0,
                                  target_affine, anat_affine,
                                  starting_affine=starting_affine)

    transformed = translation.transform(anat_data)
    regtools.overlay_slices(target_data, transformed, None, 0,
                            "target", "Transformed", "transformed_trans_0.png")
    """
    regtools.overlay_slices(target_data, transformed, None, 1,
                            "target", "Transformed", "transformed_trans_1.png")
    regtools.overlay_slices(target_data, transformed, None, 2,
                            "target", "Transformed", "transformed_trans_2.png")
    """
    transform = RigidTransform3D()
    params0 = None
    starting_affine = translation.affine
    rigid = affreg.optimize(target_data, anat_data, transform, params0,
                            target_affine, anat_affine,
                            starting_affine=starting_affine)

    transformed = rigid.transform(anat_data)
    regtools.overlay_slices(target_data, transformed, None, 0,
                            "target", "Transformed", "transformed_rigid_0.png")
    """
    regtools.overlay_slices(target_data, transformed, None, 1,
                            "nifti", "Transformed", "transformed_rigid_1.png")
    regtools.overlay_slices(target_data, transformed, None, 2,
                            "nifti", "Transformed", "transformed_rigid_2.png")
    """
    transform = AffineTransform3D()
    params0 = None
    starting_affine = rigid.affine
    affine = affreg.optimize(target_data, anat_data, transform, params0,
                             target_affine, anat_affine,
                             starting_affine=starting_affine)

    transformed = affine.transform(anat_data)
    regtools.overlay_slices(target_data, transformed, None, 0,
                            "target", "Transformed", "transformed_affine_0.png")
    """
    regtools.overlay_slices(nifti_data, transformed, None, 1,
                            "nifti_data", "Transformed", "transformed_affine_1.png")
    regtools.overlay_slices(nifti_data, transformed, None, 2,
                            "nifti_data", "Transformed", "transformed_affine_2.png")
    """