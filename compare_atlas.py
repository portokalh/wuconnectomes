

#import SimpleITK as sitk
import nibabel as nib
from dipy.io.image import load_nifti, save_nifti
import numpy as np

#from dipy.align.imaffine import (transform_centers_of_mass,
#                                 AffineMap,
#                                 MutualInformationMetric,
#                                 AffineRegistration)
from dipy.align.imaffine import AffineRegistration, transform_centers_of_mass, MutualInformationMetric, AffineMap
from dipy.align.transforms import (TranslationTransform3D)
from nilearn.image import resample_img

"""
mhdpath = "/Users/alex/jacques/atlases/atlasVolume.mhd"
niipath = "/Users/alex/jacques/atlases/test.nii"
resampledpath = "/Users/alex/jacques/atlases/whiston_resampled.nii"

img = sitk.ReadImage(mhdpath)
sitk.WriteImage(img, niipath)
"""

allen_brain = "/Volumes/Data/Badea/Lab/jacques/average_template_10.nii"
allen_brain = r"C:\Users\Jacques Stout\Documents\Work\average_template_10_reorient.nii.gz"
#resampledpath = "/Volumes/Data/Badea/Lab/jacques/resampled.nii"
#resampledpath = "/Volumes/Data/Badea/Lab/jacques/average_template_10_resampled_dipy.nii"

allen, affine_allen = load_nifti(allen_brain)

chass3_path = ("/Volumes/Data/Badea/Lab/jacques/chass_symmetric3_T1.nii.gz")
chass3_anat = load_nifti(chass3_path)

if np.size(chass3_anat[0].shape) > 3:
    chass3_data = np.squeeze(chass3_anat[0][..., 0])
else:
    chass3_data = chass3_anat[0]
chass3_affine = chass3_anat[1]

#anat_hdr = chass3_anat.header
#vox_size = chass3_anat.header.get_zooms()[0]

allen = load_nifti(allen_brain)
if np.size(allen[0].shape)>3:
    allen_data = np.squeeze(allen[0][..., 0])
else:
    allen_data = allen[0]
allen_affine = allen[1]

allen_resampled = resample_img(allen_data, target_shape=np.shape(chass3_data))

c_of_mass = transform_centers_of_mass(allen_data, allen_affine,
                                      chass3_anat, chass3_affine)

transformed = c_of_mass.transform(chass3_anat)

identity = np.eye(4)

affine_map = AffineMap(identity,
                       allen_data.shape, allen_affine,
                       chass3_data.shape, chass3_affine)
resampled = affine_map.transform(chass3_data)
save_nifti(resampledpath, resampled, chass3_affine, hdr = None)
#(fname, data, affine, hdr=None)
"""
regtools.overlay_slices(allen_data, resampled, None, 0,
                        "allen_data", "chass3_data", figspath + "resampled_0.png")
regtools.overlay_slices(allen_data, resampled, None, 1,
                        "allen_data", "chass3_data", figspath + "resampled_1.png")
regtools.overlay_slices(allen_data, resampled, None, 2,
                        "allen_data", "chass3_data", figspath + "resampled_2.png")
"""
#c_of_mass = transform_centers_of_mass(allen_data, allen_affine,
#                                      chass3_data, chass3_affine)

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







apply_niftis = []
apply_trks = []
if inputpathdir in applydirs:
    applyfiles = [anat_path]
else:
    applyfiles = []
for applydir in applydirs:
    apply_niftis.extend(get_niftis(applydir, subject))
    apply_trks.extend(get_trks(applydir, subject))


"""
if "center_mass" in registration_types:

    if apply_trks:
        metric = CCMetric(3)
        level_iters = [10, 10, 5]
        sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
        mapping = sdr.optimize(allen_data, chass3_data, allen_affine, chass3_affine,
                               c_of_mass.affine)

    for apply_nifti in apply_niftis:
        fname = os.path.basename(apply_nifti).split(".")[0]
        fpath = outputpath + fname + "_centermass.nii"
        applynii = nib.load(apply_nifti)
        apply_data = applynii.get_data()
        apply_affine = applynii.affine
        apply_hdr = myanat.header

        if len(np.shape(apply_data)) == 4:
            transformed_all = c_of_mass.transform(apply_data, apply4D=True)
            transformed = transformed_all[:, :, :, 0]
        else:
            transformed_all = c_of_mass.transform(apply_data)
            transformed = transformed_all
        save_nifti(fpath, transformed_all, apply_affine, hdr=apply_hdr)
        if figspath is not None:
            regtools.overlay_slices(allen_data, transformed, None, 0,
                                    "allen_data", "Transformed", figspath + fname + "_centermass_1.png")
            regtools.overlay_slices(allen_data, transformed, None, 1,
                                    "allen_data", "Transformed", figspath + fname + "_centermass_2.png")
            regtools.overlay_slices(allen_data, transformed, None, 2,
                                    "allen_data", "Transformed", figspath + fname + "_centermass_3.png")
        if verbose:
            print("Saved the file at " + fpath)
    # mapping = sdr.optimize(allen_data, chass3_data, allen_affine, chass3_affine,
    #                       c_of_mass.affine)
    # warped_moving = mapping.transform(chass3_data)
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
            show_template_bundles(mni_streamlines, chass3_data, show=False,
                                  fname=figspath + fname + '_streamlines_centermass.png')

        sft = StatefulTractogram(mni_streamlines, myanat, Space.RASMM)

        save_tractogram(sft, fpath, bbox_valid_check=False)
        if verbose:
            print("Saved the file at " + fpath)

metric = MutualInformationMetric(params.nbins, params.sampling_prop)

    affreg = AffineRegistration(metric=metric,
                                level_iters=params.level_iters,
                                sigmas=params.sigmas,
                                factors=params.factors)
"""