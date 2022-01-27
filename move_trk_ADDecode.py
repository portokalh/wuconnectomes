import os
from transform_handler import get_affine_transform, get_flip_affine, header_superpose, recenter_nii_affine, \
    convert_ants_vals_to_affine, read_affine_txt, recenter_nii_save, add_translation, recenter_nii_save_test, \
    affine_superpose, get_affine_transform_test, recenter_affine_test
import shutil
import nibabel as nib
import numpy as np
from dipy.align.imaffine import (MutualInformationMetric, AffineRegistration,
                                 transform_origins)
from dipy.tracking.streamline import deform_streamlines, transform_streamlines
from streamline_nocheck import load_trk, unload_trk
from dipy.io.utils import create_tractogram_header
from tract_save import save_trk_heavy_duty, make_tractogram_object, save_trk_header
from scipy.io import loadmat
from nifti_handler import extract_nii_info
import socket
from file_tools import mkcdir
from tract_handler import gettrkpath
from tract_manager import get_str_identifier
from file_tools import mkcdir, check_files
from tract_manager import check_dif_ratio
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
from nibabel.streamlines import ArraySequence
import glob, warnings
from dipy.viz import regtools

from dipy.align.imaffine import (transform_centers_of_mass,
                                 AffineMap,
                                 MutualInformationMetric,
                                 AffineRegistration)

subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361',
            'S02363', 'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451',
            'S02469', 'S02473', 'S02485', 'S02491', 'S02490', 'S02506', 'S02523', 'S02524', 'S02535', 'S02654',
            'S02666', 'S02670', 'S02686', 'S02690', 'S02695', 'S02715', 'S02720', 'S02737', 'S02745', 'S02753',
            'S02765', 'S02771', 'S02781', 'S02802', 'S02804', 'S02813', 'S02812', 'S02817', 'S02840', 'S02842',
            'S02871', 'S02877', 'S02898', 'S02926', 'S02938', 'S02939', 'S02954', 'S02967', 'S02987', 'S03010',
            'S03017', 'S03028', 'S03033', 'S03034', 'S03045', 'S03048', 'S03069', 'S03225', 'S03265', 'S03293',
            'S03308', 'S03321']

removed_list = ['S02230', 'S02490', 'S02523', 'S02745']
removed_list = removed_list + ['S03028', 'S03033', 'S03034', 'S03045', 'S03048', 'S03069', 'S03225', 'S03265', 'S03293',
                               'S03308', 'S03321']
for remove in removed_list:
    try:
        subjects.remove(remove)
    except:
        print(f'could not find {remove}')
print(subjects)
subjects = ["S03308"]
#subjects = ["S02967"]

ext = ".nii.gz"
computer_name = socket.gethostname()

if 'samos' in computer_name:
    main_path = '/mnt/paros_MRI/jacques'
elif 'santorini' in computer_name:
    main_path = '/Volumes/Data/Badea/Lab/human/AD_Decode/moving_Testing/'
else:
    raise Exception('No other computer name yet')
project = "AD_Decode"

if project == "AD_Decode":
    path_TRK = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_MPCA_100')
    path_TRK_output = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_MPCA_MDT')
    path_DWI = os.path.join(main_path, 'AD_Decode', 'Analysis', 'DWI')
    path_DWI_MDT = os.path.join(main_path, 'AD_Decode', 'Analysis', 'DWI_MDT')
    path_transforms = os.path.join(main_path, 'AD_Decode', 'Analysis', 'Transforms')
    ref = "md"
    path_trk_tempdir = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_save_alt2')
    DWI_save = os.path.join(main_path, 'AD_Decode', 'Analysis', 'NII_tempsave_alt2')
    mkcdir([path_trk_tempdir, path_TRK_output, DWI_save])
    # Get the values from DTC_launcher_ADDecode. Should probalby create a single parameter file for each project one day
    stepsize = 2
    ratio = 100
    trkroi = ["wholebrain"]
    str_identifier = get_str_identifier(stepsize, ratio, trkroi)

overwrite = False
cleanup = False
verbose = True
save_temp_files = True
recenter = 1
contrast = 'dwi'
prune = True

native_ref = ''

orient_string = os.path.join(path_DWI, 'relative_orientation.txt')
orient_relative = open(orient_string, mode='r').read()
orientation_out = orient_relative.split(',')[0]
orientation_out = orientation_out.split(':')[1]
orientation_in = orient_relative.split(',')[1]
orientation_in = orientation_in.split(':')[1]

nii_test_files = True
overwrite = True

if save_temp_files:
    mkcdir(path_trk_tempdir)

for subj in subjects:
    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    runno_to_MDT = os.path.join(path_transforms, f'{subj}_to_MDT_warp.nii.gz')
    subj_dwi_path = os.path.join(path_DWI, f'{subj}_subjspace_dwi{ext}')

    if nii_test_files:
        mkcdir([DWI_save, path_DWI_MDT])
        SAMBA_preprocess_ref = os.path.join(path_DWI, f'{subj}_labels{ext}')
        SAMBA_coreg_ref = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
        SAMBA_init = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
        SAMBA_init = subj_dwi_path
        SAMBA_preprocess = os.path.join(DWI_save, f'{subj}_{contrast}_preprocess{ext}')

        if recenter and (not os.path.exists(SAMBA_preprocess) or overwrite):
            """
            header_superpose(SAMBA_preprocess_ref, SAMBA_init, outpath=SAMBA_preprocess, verbose=False)
            img_transform_exec(SAMBA_preprocess, 'RAS', 'LPS', output_path=SAMBA_preprocess_2)
            recenter_nii_save(SAMBA_preprocess_2, SAMBA_preprocess_recentered_1, verbose=True)
            recenter_nii_save(SAMBA_preprocess_recentered_1,SAMBA_preprocess_2)
            SAMBA_init = SAMBA_preprocess_2
            """
            recenter_nii_save_test(SAMBA_init, SAMBA_preprocess)
            SAMBA_init = SAMBA_preprocess

        SAMBA_preprocess_test_posttrans = os.path.join(DWI_save, f'{subj}_{contrast}_masked_posttrans{ext}')
        SAMBA_preprocess_test_posttrans_2 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_posttrans_2{ext}')
        SAMBA_preprocess_test_posttrans_3 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_posttrans_3{ext}')

        SAMBA_preprocess_test_rigid = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid{ext}')
        SAMBA_preprocess_test_rigid_affine = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid_affine{ext}')
        SAMBA_preprocess_test_postwarp = os.path.join(path_DWI_MDT, f'{subj}_{contrast}_postwarp{ext}')
        if native_ref == '':
            native_ref = SAMBA_init
        if not os.path.exists(SAMBA_preprocess_test_postwarp) or overwrite:
            cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_init} -r {SAMBA_init}  -n Linear  -o {SAMBA_preprocess_test_posttrans}'
            os.system(cmd)

            affine_superpose(SAMBA_init, SAMBA_preprocess_test_posttrans, outpath=SAMBA_preprocess_test_posttrans_2)

            cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_preprocess_test_posttrans_2} -r {SAMBA_preprocess_test_posttrans_2}  -n Linear  -o {SAMBA_preprocess_test_posttrans_3} -t {trans}'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_posttrans_3} -o {SAMBA_preprocess_test_rigid} ' \
                f'-r {SAMBA_preprocess_test_posttrans_2} -n Linear -t [{rigid},0]'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_rigid} -o {SAMBA_preprocess_test_rigid_affine} ' \
                f'-r {SAMBA_preprocess_test_posttrans_2} -n Linear -t [{affine_orig},0]'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_rigid_affine} -o {SAMBA_preprocess_test_postwarp} ' \
                f'-r {SAMBA_preprocess_test_posttrans_2} -n Linear -t {runno_to_MDT}'
            os.system(cmd)

    trk_preprocess_posttrans = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_posttrans.trk')
    trk_preprocess_postrigid = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_postrigid.trk')
    trk_preprocess_postrigid_affine = os.path.join(path_trk_tempdir,
                                                   f'{subj}{str_identifier}_preprocess_postrigid_affine.trk')
    trk_MDT_space = os.path.join(path_TRK_output, f'{subj}_MDT.trk')

    overwrite = True
    if not os.path.exists(trk_MDT_space) or overwrite:

        subj_dwi_path = os.path.join(path_DWI, f'{subj}_subjspace_dwi{ext}')
        subj_dwi = nib.load(subj_dwi_path)
        subj_dwi_data = subj_dwi.get_data()
        subj_affine = subj_dwi._affine

        SAMBA_input_real_file = os.path.join(path_DWI, f'{subj}_dwi{ext}')

        check_dif_ratio(path_TRK, subj, str_identifier, ratio)
        subj_trk, trkexists = gettrkpath(path_TRK, subj, str_identifier, pruned=prune, verbose=False)

        _, exists = check_files([trans, rigid, runno_to_MDT])
        if np.any(exists == 0):
            raise Exception('missing transform file')
        if not os.path.exists(affine) and not os.path.exists(affine_orig):
            raise Exception('missing transform file')

        subj_affine_new = recenter_affine_test(np.shape(subj_dwi_data), subj_affine)

        subj_centered = np.copy(subj_affine)
        subj_centered[:3, 3] = [0, 0, 0]
        subj_dwi_c = nib.Nifti1Image(subj_dwi_data, subj_centered)
        subj_dwi_c_path = os.path.join(DWI_save, f'{subj}_subjspace_centered{ext}')

        subj_centered_rotated = np.copy(subj_affine_new)
        subj_centered_rotated[:3, 3] = [0, 0, 0]
        subj_dwi_cr = nib.Nifti1Image(subj_dwi_data, subj_centered_rotated)
        subj_dwi_cr_path = os.path.join(DWI_save, f'{subj}_subjspace_centered_rotated{ext}')

        preprocess_target = f'/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/preprocess/{subj}_dwi_masked.nii.gz'
        subj_affine_new = recenter_affine_test(np.shape(subj_dwi_data), subj_affine)
        subj_torecenter_transform_affine = get_affine_transform_test(subj_affine, subj_affine_new)
        streamlines, header = unload_trk(subj_trk)

        c_of_mass = transform_centers_of_mass(subj_dwi_cr._data, subj_dwi_cr.affine,
                                              subj_dwi_c._data, subj_dwi_c.affine)

        subj_affine_eye = np.eye(4)
        subj_affine_eye[:3, :3] = subj_torecenter_transform_affine[:3, :3]
        subj_affine_eye_inv = np.linalg.inv(subj_affine_eye)
        nii_shape = np.array(np.shape(subj_dwi_data))

        subj_affine_eye_inv[:3, 3] = c_of_mass.affine_inv[:3, 3]

        """
        subj_eye_streamlines = transform_streamlines(streamlines, subj_affine_eye_inv,
                                                     in_place=False)
        trk_subjspace_eye = os.path.join(path_trk_tempdir, f'{subj}_subjspace_eye_ctest.trk')

        if (not os.path.exists(trk_subjspace_eye) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_subjspace_eye, streamlines=subj_eye_streamlines, header=header,
                            affine=np.eye(4), verbose=verbose)
        """

        subj_torecenter_inv = np.eye(4)
        subj_torecenter_inv[:3, :3] = np.linalg.inv(subj_torecenter_transform_affine[:3, :3])
        center_array = [element * 0.5 for element in nii_shape + [1, 1, 1]]
        recentering_translation = subj_affine[:3, 3] - np.multiply(center_array, [1, 1, -1])
        subj_torecenter_inv[:3, 3] = c_of_mass.affine_inv[:3, 3] - subj_affine[:3, 3] + np.multiply(center_array,
                                                                                                    [1, 1, -1])

        trk_preprocess_path = os.path.join(path_trk_tempdir, f'{subj}_preprocess.trk')
        streamlines_prepro = transform_streamlines(streamlines, subj_torecenter_inv,
                                                          in_place=False)
        if (not os.path.exists(trk_preprocess_path) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_path, streamlines=streamlines_prepro, header=header,
                            affine=np.eye(4), verbose=verbose)

        #streamlines_prepro, header = unload_trk(trk_preprocess_path)

        mat_struct = loadmat(trans)
        var_name = list(mat_struct.keys())[0]
        later_trans_mat = mat_struct[var_name]
        new_transmat = np.eye(4)
        vox_dim = [1, 1, -1]
        #new_transmat[:3, 3] = np.squeeze(later_trans_mat[3:6]) * vox_dim
        new_transmat[:3, 3] = np.squeeze(np.matmul(subj_affine[:3, :3], later_trans_mat[3:6])) #should be the AFFINE of the current image, to make sure the slight difference in orientation is ACCOUNTED FOR!!!!!!!!!!
        new_transmat[:3, 3] = np.squeeze(np.matmul(subj_torecenter_transform_affine[:3, :3], later_trans_mat[3:6]))
        new_transmat[:3, 3] = np.squeeze(later_trans_mat[3:6])
        #new_transmat[2, 3] = 0
        print(new_transmat)
        streamlines_posttrans = transform_streamlines(streamlines_prepro, (new_transmat))

        if (not os.path.exists(trk_preprocess_posttrans) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_posttrans, streamlines=streamlines_posttrans, header=header,
                            affine=np.eye(4), verbose=verbose)

        rigid_struct = loadmat(rigid)
        var_name = list(rigid_struct.keys())[0]
        rigid_ants = rigid_struct[var_name]
        rigid_mat = convert_ants_vals_to_affine(rigid_ants)

        # streamlines_posttrans, header_posttrans = unload_trk(trk_preprocess_posttrans)
        streamlines_postrigid = transform_streamlines(streamlines_posttrans, np.linalg.inv(rigid_mat))

        if (not os.path.exists(trk_preprocess_postrigid) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_postrigid, streamlines=streamlines_postrigid, header=header,
                            affine=np.eye(4), verbose=verbose)

        if os.path.exists(affine):
            affine_mat_s = read_affine_txt(affine)
        else:
            cmd = f'ConvertTransformFile 3 {affine_orig} {affine} --matrix'
            os.system(cmd)
            affine_mat_s = read_affine_txt(affine)

        affine_struct = loadmat(affine_orig)
        var_name = list(affine_struct.keys())[0]
        affine_ants = affine_struct[var_name]
        affine_mat = convert_ants_vals_to_affine(affine_ants)

        streamlines_postrigidaffine = transform_streamlines(streamlines_postrigid, np.linalg.inv(affine_mat))

        if (not os.path.exists(trk_preprocess_postrigid_affine) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_postrigid_affine, streamlines=streamlines_postrigidaffine,
                            header=header,
                            affine=np.eye(4), verbose=verbose)

        # streamlines_postrigidaffine, header_postrigidaffine = unload_trk(trk_preprocess_postrigid_affine)

        warp, affine, vox_size, header_warp, ref_info = extract_nii_info(runno_to_MDT)
        warp = warp[:, :, :, 0, :]

        vox_size = tuple(np.linalg.eigvals(affine[:3, :3]))
        vox_size = vox_size[0]
        target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
        streamlines_post_warp = deform_streamlines(
            streamlines_postrigidaffine, deform_field=warp,
            stream_to_current_grid=target_isocenter,
            current_grid_to_world=np.eye(4), stream_to_ref_grid=target_isocenter,
            ref_grid_to_world=np.eye(4))

        if (not os.path.exists(trk_MDT_space) or overwrite):
            save_trk_header(filepath=trk_MDT_space, streamlines=streamlines_post_warp, header=header,
                            affine=np.eye(4), fix_streamlines=False, verbose=verbose)
    else:
        print(f'{trk_MDT_space} already exists')