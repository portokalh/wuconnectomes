"""
Create by Jacques STOUT
VERY WIP
Used to register TRKs
Mostly taken from dipy advice/tutorials, will probably be completely deprecated in future
"""

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

from dipy.io.image import save_nifti
import nibabel as nib
from file_catcher import get_anat, get_niftis, get_trks
from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.tracking.streamline import deform_streamlines
from dipy.viz import has_fury
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.metrics import CCMetric
import os

from dipy.viz import window, actor


def show_template_bundles(bundles, static, show=True, fname=None):
    scene = window.Scene()
    template_actor = actor.slicer(static)
    scene.add(template_actor)

    lines_actor = actor.streamtube(bundles, window.colors.orange,
                                   linewidth=0.3)
    scene.add(lines_actor)

    if show:
        window.show(scene)
    if fname is not None:
        window.record(scene, n_frames=1, out_path=fname, size=(900, 900))


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
        self.sampling_prop = sampling_proportion
        self.level_iters = level_iters
        self.sigmas = sigmas
        self.factors = factors

def register_save(inputpathdir, target_path, subject, outputpath, figspath, params, registration_types, applydirs,
                  verbose):
    anat_path = get_anat(inputpathdir, subject)
    #myanat = load_nifti(anat_path)
    myanat = nib.load(anat_path)
    anat_data = np.squeeze(myanat.get_data()[..., 0])
    anat_affine = myanat.affine
    anat_hdr = myanat.header
    vox_size = myanat.header.get_zooms()[0]
    #mynifti = load_nifti("/Volumes/Data/Badea/Lab/19abb14/N57437_nii4D.nii")
    #anat_data = np.squeeze(myanat[0])[..., 0]
    #anat_affine = myanat[1]
    #hdr = myanat.header

    mytarget = nib.load(target_path)
    target_data = np.squeeze(mytarget.get_data()[..., 0])
    target_affine = mytarget.affine

    identity = np.eye(4)

    affine_map = AffineMap(identity,
                           target_data.shape, target_affine,
                           anat_data.shape, anat_affine)
    resampled = affine_map.transform(anat_data)

    """
    regtools.overlay_slices(target_data, resampled, None, 0,
                            "target_data", "anat_data", figspath + "resampled_0.png")
    regtools.overlay_slices(target_data, resampled, None, 1,
                            "target_data", "anat_data", figspath + "resampled_1.png")
    regtools.overlay_slices(target_data, resampled, None, 2,
                            "target_data", "anat_data", figspath + "resampled_2.png")
    """
    c_of_mass = transform_centers_of_mass(target_data, target_affine,
                                          anat_data, anat_affine)
    apply_niftis = []
    apply_trks = []
    if inputpathdir in applydirs:
        applyfiles = [anat_path]
    else:
        applyfiles = []
    for applydir in applydirs:
        apply_niftis.extend(get_niftis(applydir, subject))
        apply_trks.extend(get_trks(applydir, subject))

    if "center_mass" in registration_types:

        if apply_trks:
            metric = CCMetric(3)
            level_iters = [10, 10, 5]
            sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
            mapping = sdr.optimize(target_data, anat_data, target_affine, anat_affine,
                                   c_of_mass.affine)

        for apply_nifti in apply_niftis:
            fname = os.path.basename(apply_nifti).split(".")[0]
            fpath = outputpath + fname + "_centermass.nii"
            applynii = nib.load(apply_nifti)
            apply_data = applynii.get_data()
            apply_affine = applynii.affine
            apply_hdr = myanat.header
            
            if len(np.shape(apply_data)) == 4:
                transformed_all = c_of_mass.transform(apply_data,apply4D=True)
                transformed = transformed_all[:, :, :, 0]
            else:
                transformed_all = c_of_mass.transform(apply_data)
                transformed = transformed_all
            save_nifti(fpath, transformed_all, apply_affine, hdr=apply_hdr)
            if figspath is not None:
                regtools.overlay_slices(target_data, transformed, None, 0,
                                        "target_data", "Transformed", figspath + fname + "_centermass_1.png")
                regtools.overlay_slices(target_data, transformed, None, 1,
                                        "target_data", "Transformed", figspath + fname + "_centermass_2.png")
                regtools.overlay_slices(target_data, transformed, None, 2,
                                        "target_data", "Transformed", figspath + fname + "_centermass_3.png")
            if verbose:
                print("Saved the file at " + fpath)
        #mapping = sdr.optimize(target_data, anat_data, target_affine, anat_affine,
        #                       c_of_mass.affine)
        #warped_moving = mapping.transform(anat_data)
        for apply_trk in apply_trks:

            fname = os.path.basename(apply_trk).split(".")[0]
            fpath = outputpath + fname + "_centermass.trk"

            sft = load_tractogram(apply_trk, 'same')
            target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
            origin_affine = affine_map.affine.copy()
            origin_affine[0][3] = -origin_affine[0][3]
            origin_affine[1][3] = -origin_affine[1][3]
            origin_affine[2][3] = origin_affine[2][3] / vox_size

            origin_affine[1][3] = origin_affine[1][3] / vox_size ** 2

            # Apply the deformation and correct for the extents
            mni_streamlines = deform_streamlines(
                sft.streamlines, deform_field=mapping.get_forward_field(),
                stream_to_current_grid=target_isocenter,
                current_grid_to_world=origin_affine, stream_to_ref_grid=target_isocenter,
                ref_grid_to_world=np.eye(4))

            if has_fury:

                show_template_bundles(mni_streamlines, anat_data, show=False,
                                      fname = figspath + fname + '_streamlines_centermass.png')

            sft = StatefulTractogram(mni_streamlines, myanat, Space.RASMM)

            save_tractogram(sft, fpath, bbox_valid_check=False)
            if verbose:
                print("Saved the file at " + fpath)

    metric = MutualInformationMetric(params.nbins, params.sampling_prop)

    if "AffineRegistration" in registration_types:
        affreg = AffineRegistration(metric=metric,
                                    level_iters=params.level_iters,
                                    sigmas=params.sigmas,
                                    factors=params.factors)

        transform = TranslationTransform3D()
        params0 = None
        starting_affine = c_of_mass.affine
        translation = affreg.optimize(target_data, anat_data, transform, params0,
                                      target_affine, anat_affine,
                                      starting_affine=starting_affine)

        if apply_trks:
            metric = CCMetric(3)
            level_iters = [10, 10, 5]
            sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
            mapping = sdr.optimize(target_data, anat_data, target_affine, anat_affine,
                                   translation.affine)

        for apply_nifti in apply_niftis:
            fname = os.path.basename(apply_nifti).split(".")[0]
            fpath = outputpath + fname + "_affinereg.nii"

            applynii = nib.load(apply_nifti)
            apply_data = applynii.get_data()
            apply_affine = applynii.affine
            apply_hdr = myanat.header

            if len(np.shape(apply_data)) == 4:
                transformed_all = translation.transform(apply_data, apply4D=True)
                transformed = transformed_all[:, :, :, 0]
            else:
                transformed_all = translation.transform(apply_data)
                transformed = transformed_all
            save_nifti(fpath, transformed_all, anat_affine, hdr=anat_hdr)
            if figspath is not None:
                regtools.overlay_slices(target_data, transformed, None, 0,
                                        "target_data", "Transformed", figspath + fname + "_affinereg_1.png")
                regtools.overlay_slices(target_data, transformed, None, 1,
                                        "target_data", "Transformed", figspath + fname + "_affinereg_2.png")
                regtools.overlay_slices(target_data, transformed, None, 2,
                                        "target_data", "Transformed", figspath + fname + "_affinereg_3.png")
            if verbose:
                print("Saved the file at " + fpath)

        for apply_trk in apply_trks:

            fname = os.path.basename(apply_trk).split(".")[0]
            fpath = outputpath + fname + "_affinereg.trk"

            sft = load_tractogram(apply_trk, 'same')
            target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
            origin_affine = affine_map.affine.copy()
            origin_affine[0][3] = -origin_affine[0][3]
            origin_affine[1][3] = -origin_affine[1][3]
            origin_affine[2][3] = origin_affine[2][3] / vox_size

            origin_affine[1][3] = origin_affine[1][3] / vox_size ** 2

            # Apply the deformation and correct for the extents
            mni_streamlines = deform_streamlines(
                sft.streamlines, deform_field=mapping.get_forward_field(),
                stream_to_current_grid=target_isocenter,
                current_grid_to_world=origin_affine, stream_to_ref_grid=target_isocenter,
                ref_grid_to_world=np.eye(4))

            if has_fury:

                show_template_bundles(mni_streamlines, anat_data, show=False,
                                      fname= figspath + fname + '_streamlines_affinereg.png')

            sft = StatefulTractogram(mni_streamlines, myanat, Space.RASMM)

            save_tractogram(sft, fpath, bbox_valid_check=False)
            if verbose:
                print("Saved the file at " + fpath)


    if "RigidTransform3D" in registration_types:
        transform = RigidTransform3D()
        params0 = None
        if 'translation' not in locals():
            affreg = AffineRegistration(metric=metric,
                                        level_iters=params.level_iters,
                                        sigmas=params.sigmas,
                                        factors=params.factors)
            translation = affreg.optimize(target_data, anat_data, transform, params0,
                                          target_affine, anat_affine,
                                          starting_affine=c_of_mass.affine)
        starting_affine = translation.affine
        rigid = affreg.optimize(target_data, anat_data, transform, params0,
                                target_affine, anat_affine,
                                starting_affine=starting_affine)

        transformed = rigid.transform(anat_data)

        if apply_trks:
            metric = CCMetric(3)
            level_iters = [10, 10, 5]
            sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
            mapping = sdr.optimize(target_data, anat_data, target_affine, anat_affine,
                                   rigid.affine)

        for apply_nifti in apply_niftis:
            fname = os.path.basename(apply_nifti).split(".")[0]
            fpath = outputpath + fname + "_rigidtransf3d.nii"

            applynii = nib.load(apply_nifti)
            apply_data = applynii.get_data()
            apply_affine = applynii.affine
            apply_hdr = myanat.header

            if len(np.shape(apply_data)) == 4:
                transformed_all = rigid.transform(apply_data, apply4D=True)
                transformed = transformed_all[:, :, :, 0]
            else:
                transformed_all = rigid.transform(apply_data)
                transformed = transformed_all
            save_nifti(fpath, transformed_all, anat_affine, hdr=anat_hdr)
            if figspath is not None:
                regtools.overlay_slices(target_data, transformed, None, 0,
                                        "target_data", "Transformed", figspath + fname + "_rigidtransf3d_1.png")
                regtools.overlay_slices(target_data, transformed, None, 1,
                                        "target_data", "Transformed", figspath + fname + "_rigidtransf3d_2.png")
                regtools.overlay_slices(target_data, transformed, None, 2,
                                        "target_data", "Transformed", figspath + fname + "_rigidtransf3d_3.png")
            if verbose:
                print("Saved the file at " + fpath)

        for apply_trk in apply_trks:

            fname = os.path.basename(apply_trk).split(".")[0]
            fpath = outputpath + fname + "_rigidtransf3d.trk"

            sft = load_tractogram(apply_trk, 'same')
            target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
            origin_affine = affine_map.affine.copy()
            origin_affine[0][3] = -origin_affine[0][3]
            origin_affine[1][3] = -origin_affine[1][3]
            origin_affine[2][3] = origin_affine[2][3] / vox_size

            origin_affine[1][3] = origin_affine[1][3] / vox_size ** 2

            # Apply the deformation and correct for the extents
            mni_streamlines = deform_streamlines(
                sft.streamlines, deform_field=mapping.get_forward_field(),
                stream_to_current_grid=target_isocenter,
                current_grid_to_world=origin_affine, stream_to_ref_grid=target_isocenter,
                ref_grid_to_world=np.eye(4))

            if has_fury:

                show_template_bundles(mni_streamlines, anat_data, show=False,
                                      fname= figspath + fname + '_rigidtransf3d.png')

            sft = StatefulTractogram(mni_streamlines, myanat, Space.RASMM)

            save_tractogram(sft, fpath, bbox_valid_check=False)
            if verbose:
                print("Saved the file at " + fpath)

    if "AffineTransform3D" in registration_types:
        transform = AffineTransform3D()
        params0 = None
        starting_affine = rigid.affine
        affine = affreg.optimize(target_data, anat_data, transform, params0,
                                 target_affine, anat_affine,
                                 starting_affine=starting_affine)

        transformed = affine.transform(anat_data)

        if apply_trks:
            metric = CCMetric(3)
            level_iters = [10, 10, 5]
            sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
            mapping = sdr.optimize(target_data, anat_data, target_affine, anat_affine,
                                   affine.affine)

        for apply_nifti in apply_niftis:
            fname = os.path.basename(apply_nifti).split(".")[0]
            fpath = outputpath + fname + "_affinetransf3d.nii"

            applynii = nib.load(apply_nifti)
            apply_data = applynii.get_data()
            apply_affine = applynii.affine
            apply_hdr = myanat.header

            if len(np.shape(apply_data)) == 4:
                transformed_all = affine.transform(apply_data, apply4D=True)
                transformed = transformed_all[:, :, :, 0]
            else:
                transformed_all = affine.transform(apply_data)
                transformed = transformed_all
            save_nifti(fpath, transformed_all, anat_affine, hdr=anat_hdr)
            if figspath is not None:
                regtools.overlay_slices(target_data, transformed, None, 0,
                                        "target_data", "Transformed", figspath + fname + "_affinetransf3d_1.png")
                regtools.overlay_slices(target_data, transformed, None, 1,
                                        "target_data", "Transformed", figspath + fname + "_affinetransf3d_2.png")
                regtools.overlay_slices(target_data, transformed, None, 2,
                                        "target_data", "Transformed", figspath + fname + "_affinetransf3d_3.png")
            if verbose:
                print("Saved the file at " + fpath)

        for apply_trk in apply_trks:

            fname = os.path.basename(apply_trk).split(".")[0]
            fpath = outputpath + fname + "_affinetransf3d.trk"

            sft = load_tractogram(apply_trk, 'same')
            target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
            origin_affine = affine_map.affine.copy()
            origin_affine[0][3] = -origin_affine[0][3]
            origin_affine[1][3] = -origin_affine[1][3]
            origin_affine[2][3] = origin_affine[2][3] / vox_size

            origin_affine[1][3] = origin_affine[1][3] / vox_size ** 2

            # Apply the deformation and correct for the extents
            mni_streamlines = deform_streamlines(
                sft.streamlines, deform_field=mapping.get_forward_field(),
                stream_to_current_grid=target_isocenter,
                current_grid_to_world=origin_affine, stream_to_ref_grid=target_isocenter,
                ref_grid_to_world=np.eye(4))

            if has_fury:

                show_template_bundles(mni_streamlines, anat_data, show=False,
                                      fname= figspath + fname + '_streamlines_affinetransf3d.png')

            sft = StatefulTractogram(mni_streamlines, myanat, Space.RASMM)

            save_tractogram(sft, fpath, bbox_valid_check=False)
            if verbose:
                print("Saved the file at " + fpath)