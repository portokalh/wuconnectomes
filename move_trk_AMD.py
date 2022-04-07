import os
from transform_handler import get_affine_transform, get_flip_affine, header_superpose, recenter_nii_affine, \
    convert_ants_vals_to_affine, read_affine_txt, recenter_nii_save, add_translation, recenter_nii_save_test, \
    affine_superpose, get_affine_transform_test
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
from file_tools import mkcdir, check_files, getfromfile
from tract_manager import check_dif_ratio
from computer_nav import get_mainpaths, load_nifti_remote, load_trk_remote, loadmat_remote, checkfile_exists_remote
import glob, warnings
import random

project= 'AMD'
subjects = ["H26578", "H29060", "H26637", "H29264", "H26765", "H29225", "H26660", "H29304", "H26890",
            "H29556", "H26862", "H29410", "H26966", "H29403", "H26841", "H21593", "H27126", "H29618", "H27111", "H29627",
            "H27164", "H29502", "H27100", "H27381", "H21836", "H27391", "H21850", "H27495", "H21729", "H27488", "H21915",
            "H27682", "H21956", "H27686", "H22331", "H28208", "H21990", "H28955", "H29878", "H27719", "H22102", "H27841",
            "H22101", "H27842", "H22228", "H28029", "H22140", "H27852", "H22276", "H27999", "H22369", "H28115", "H22644",
            "H28308", "H22574", "H28377", "H22368", "H28325", "H22320", "H28182", "H22898", "H28748", "H22683", "H28373",
            "H22536", "H28433", "H22825", "H28662", "H22864", "H28698", "H23143", "H28861", "H23157", "H28820", "H23028",
            "H29002", "H23210", "H29020", "H23309", "H29161", "H26841", "H26862", "H26949", "H27100", "H27126",
            "H27163", "H27246", "H27488", "H27682", "H27686", "H27719", "H27841", "H27842", "H27852", "H27869", "H27999",
            "H28029", "H28068", "H28208", "H28262", "H28325", "H28820", "H28856", "H28869", "H28955", "H29002", "H29044",
            "H29089", "H29127", "H29161", "H29242", "H29254", "H26578", "H26637", "H26660", "H26745", "H26765", "H26850",
            "H26880", "H26890", "H26958", "H26974", "H27017", "H27111", "H27164", "H27381", "H27391", "H27495", "H27610",
            "H27640", "H27680", "H27778", "H27982", "H28115", "H28308", "H28338", "H28373", "H28377", "H28433", "H28437",
            "H28463", "H28532", "H28662", "H28698", "H28748", "H28809", "H28857", "H28861", "H29013", "H29020", "H29025"
            ,"H29056"]
#subjects = ["H26660"]
#subjects = ["H29410"]
subjects = sorted(subjects)
random.shuffle(subjects)
#temporarily removing "H29056" to recalculate it
ext = ".nii.gz"
computer_name = socket.gethostname()

username=None
passwd = None
if 'samos' not in computer_name:
    remote=True
else:
    remote=False
if remote:
    username, passwd = getfromfile('/Users/jas/samos_connect.rtf')
inpath, outpath, _, sftp = get_mainpaths(remote,project = 'AMD', username=username,password=passwd)


testing = False
if testing:
    inpath = '/Users/alex/jacques/AMD_testing_zone/AMD_TRK_testing'

if project == "AMD":
    path_TRK = os.path.join(inpath, 'TRK')
    path_DWI = os.path.join(inpath, 'DWI')
    path_transforms = os.path.join(inpath, 'Transforms')
    ref = "md"
    path_trk_tempdir = os.path.join(outpath, 'TRK_transition')
    path_TRK_output = os.path.join(outpath, 'TRK_rigidaff')
    DWI_save = os.path.join(outpath, 'NII_toMDT')
    #Get the values from DTC_launcher_ADDecode. Should probably create a single parameter file for each project one day

mkcdir([path_trk_tempdir,path_TRK_output, DWI_save],sftp)

stepsize = 2
ratio = 1
trkroi = ["wholebrain"]
str_identifier = get_str_identifier(stepsize, ratio, trkroi)
prune= True
overwrite = False
cleanup = False
verbose=True
recenter=True

save_temp_trk_files = False
nii_test_files = False

contrast='dwi'
native_ref=''


if save_temp_trk_files:
    mkcdir(path_trk_tempdir, sftp)

for subj in subjects:

    print(f'running move for subject {subj}')
    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    runno_to_MDT = os.path.join(path_transforms, f'{subj}_to_MDT_warp.nii.gz')
    subj_dwi = os.path.join(path_DWI, f'{subj}_subjspace_dwi{ext}')

    if nii_test_files:
        mkcdir(DWI_save)
        SAMBA_preprocess_ref = os.path.join(path_DWI, f'{subj}_labels{ext}')
        SAMBA_coreg_ref = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
        SAMBA_init = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
        SAMBA_init = subj_dwi
        SAMBA_preprocess = os.path.join(DWI_save, f'{subj}_{contrast}_preprocess{ext}')
        #SAMBA_preprocess_2 = os.path.join(DWI_save, f'{subj}_{contrast}_preprocess_2{ext}')
        #SAMBA_preprocess_recentered_1 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_1{ext}')
        #SAMBA_preprocess_recentered_2 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_2{ext}')

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
        SAMBA_preprocess_test_postwarp = os.path.join(DWI_save, f'{subj}_{contrast}_postwarp{ext}')
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
    #trk_preprocess_postrigid_affine = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_postrigid_affine.trk')
    trk_preprocess_postrigid_affine = os.path.join(path_TRK_output, f'{subj}{str_identifier}_preprocess_postrigid_affine.trk')
    trk_MDT_space = os.path.join(path_TRK_output, f'{subj}_MDT.trk')

    #final_img_exists = checkfile_exists_remote(trk_MDT_space)
    final_img_exists = checkfile_exists_remote(trk_preprocess_postrigid_affine, sftp)

    if not final_img_exists or overwrite:

        _, exists = check_files([trans, rigid, runno_to_MDT], sftp)
        print('reaching point 1')

        subj_dwi = os.path.join(path_DWI, f'{subj}_coreg_diff{ext}')
        nii_data, subj_affine, _, _, _= load_nifti_remote(subj_dwi,sftp)
        subj_affine_new = subj_affine

        #subj_torecenter_transform_affine = get_affine_transform_test(subj_affine, subj_affine_new)
        #added_trans = subj_affine[:3, 3] + np.multiply(preprocess_affine[:3, 3], [1,1,-1]) + [-1,-1,0]
        #added_trans = subj_affine[:3, 3] + np.multiply(subjtorecenter_affine[:3, 3], [-1,-1,-1])
        #subj_torecenter_transform_affine[:3, 3] = reorient_trans + added_trans

        #SAMBA_input_real_file =  os.path.join(path_DWI, f'{subj}_dwi{ext}')
        #new_affine, translation, translate_affine = recenter_nii_affine(SAMBA_input_real_file, return_translation=True)

        #subj_trk, _ = gettrkpath(path_TRK, subj, str_identifier, pruned=True, verbose=verbose)
        print('reaching point 2')
        check_dif_ratio(path_TRK, subj, str_identifier, ratio, sftp)
        subj_trk, trkexists = gettrkpath(path_TRK, subj, str_identifier, pruned=prune, verbose=False, sftp=sftp)
        if not trkexists:
            txt = f'Could not find TRK file for subject {subj}'
            raise Exception(txt)
        _, exists = check_files([trans, rigid, runno_to_MDT], sftp)
        if np.any(exists==0):
            raise Exception('missing transform file')
        _, exists = check_files([affine, affine_orig], sftp)
        if np.any(exists==0):
            raise Exception('missing transform file')
        #streamlines_prepro, header = unload_trk(subj_trk)
        streamlines_data = load_trk_remote(subj_trk, 'same', sftp)
        streamlines_prepro = streamlines_data.streamlines
        header = streamlines_data.space_attributes
        #streamlines_prepro, header_prepro = unload_trk(trk_preprocess)
        mat_struct = loadmat_remote(trans, sftp)
        var_name = list(mat_struct.keys())[0]
        later_trans_mat = mat_struct[var_name]
        new_transmat = np.eye(4)
        vox_dim = [1, 1, -1]
        #new_transmat[:3, 3] = np.squeeze(later_trans_mat[3:6]) * vox_dim
        new_transmat[:3, 3] = np.squeeze(np.matmul(subj_affine[:3, :3], later_trans_mat[3:6])) #should be the AFFINE of the current image, to make sure the slight difference in orientation is ACCOUNTED FOR!!!!!!!!!!
        new_transmat[2, 3] = 0
        print(new_transmat)
        streamlines_posttrans = transform_streamlines(streamlines_prepro, (new_transmat))
        print('reaching point 3')
        if (not checkfile_exists_remote(trk_preprocess_posttrans, sftp) or overwrite) and save_temp_trk_files:
            save_trk_header(filepath= trk_preprocess_posttrans, streamlines = streamlines_posttrans, header = header,
                    affine=np.eye(4), verbose=verbose, sftp=sftp)

        rigid_struct = loadmat_remote(rigid,sftp)
        var_name = list(rigid_struct.keys())[0]
        rigid_ants = rigid_struct[var_name]
        rigid_mat = convert_ants_vals_to_affine(rigid_ants)

        #streamlines_posttrans, header_posttrans = unload_trk(trk_preprocess_posttrans)
        streamlines_postrigid = transform_streamlines(streamlines_posttrans, np.linalg.inv(rigid_mat))
        print('reaching point 4')
        if (not checkfile_exists_remote(trk_preprocess_postrigid, sftp) or overwrite) and save_temp_trk_files:
            save_trk_header(filepath=trk_preprocess_postrigid, streamlines=streamlines_postrigid, header=header,
                    affine=np.eye(4), verbose=verbose, sftp=sftp)

        """
        if os.path.exists(affine):
            affine_mat_s = read_affine_txt(affine)
        else:
            cmd = f'ConvertTransformFile 3 {affine_orig} {affine} --matrix'
            os.system(cmd)
            affine_mat_s = read_affine_txt(affine)
        """

        affine_struct = loadmat_remote(affine_orig,sftp)
        var_name = list(affine_struct.keys())[0]
        affine_ants = affine_struct[var_name]
        affine_mat = convert_ants_vals_to_affine(affine_ants)

        streamlines_postrigidaffine = transform_streamlines(streamlines_postrigid, np.linalg.inv(affine_mat))
        print('reaching point 5')
        if (not checkfile_exists_remote(trk_preprocess_postrigid_affine, sftp) or overwrite):
            save_trk_header(filepath=trk_preprocess_postrigid_affine, streamlines=streamlines_postrigidaffine, header=header,
                    affine=np.eye(4), verbose=verbose, sftp=sftp)
            print('did this')
        print('finished')
        #streamlines_postrigidaffine, header_postrigidaffine = unload_trk(trk_preprocess_postrigid_affine)
        """
        warp, affine, vox_size, header_warp, ref_info = extract_nii_info(runno_to_MDT)
        warp = warp[:,:,:,0,:]

        vox_size = tuple(np.linalg.eigvals(affine[:3,:3]))
        vox_size = vox_size[0]
        target_isocenter = np.diag(np.array([-vox_size, vox_size, vox_size, 1]))
        streamlines_post_warp = deform_streamlines(
            streamlines_postrigidaffine, deform_field=warp,
            stream_to_current_grid=target_isocenter,
            current_grid_to_world=np.eye(4), stream_to_ref_grid=target_isocenter,
            ref_grid_to_world=np.eye(4))

        if (not os.path.exists(trk_MDT_space) or overwrite):
            save_trk_header(filepath=trk_MDT_space, streamlines=streamlines_post_warp, header=header,
                    affine=np.eye(4), fix_streamlines=False, verbose=verbose, sftp=sftp)
        """
    else:
        print(f'{trk_MDT_space} already exists')