import os
from transform_handler import get_affine_transform, get_flip_affine, header_superpose, recenter_nii_affine, \
    convert_ants_vals_to_affine, read_affine_txt, recenter_nii_save, add_translation
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
import glob, warnings

subjects = ["H26578", "H29060", "H26637", "H29264", "H26765", "H29225", "H26660", "H29304", "H26890",
            "H29556", "H26862", "H29410", "H26966", "H29403", "H26841", "H21593", "H27126", "H29618", "H27111", "H29627",
            "H27164", "H29502", "H27100", "H27381", "H21836", "H27391", "H21850", "H27495", "H21729", "H27488", "H21915",
            "H27682", "H21956", "H27686", "H22331", "H28208", "H21990", "H28955", "H29878", "H27719", "H22102", "H27841",
            "H22101", "H27842", "H22228", "H28029", "H22140", "H27852", "H22276", "H27999", "H22369", "H28115", "H22644",
            "H28308", "H22574", "H28377", "H22368", "H28325", "H22320", "H28182", "H22898", "H28748", "H22683", "H28373",
            "H22536", "H28433", "H22825", "H28662", "H22864", "H28698", "H23143", "H28861", "H23157", "H28820", "H23028",
            "H29002", "H23210", "H29020", "H23309", "H29161", "H26841", "H26862", "H26949", "H26966", "H27100", "H27126",
            "H27163", "H27246", "H27488", "H27682", "H27686", "H27719", "H27841", "H27842", "H27852", "H27869", "H27999",
            "H28029", "H28068", "H28208", "H28262", "H28325", "H28820", "H28856", "H28869", "H28955", "H29002", "H29044",
            "H29089", "H29127", "H29161", "H29242", "H29254", "H26578", "H26637", "H26660", "H26745", "H26765", "H26850",
            "H26880", "H26890", "H26958", "H26974", "H27017", "H27111", "H27164", "H27381", "H27391", "H27495", "H27610",
            "H27640", "H27680", "H27778", "H27982", "H28115", "H28308", "H28338", "H28373", "H28377", "H28433", "H28437",
            "H28463", "H28532", "H28662", "H28698", "H28748", "H28809", "H28857", "H28861", "H29013", "H29020", "H29025"]

#temporarily removing "H29056" to recalculate it
ext = ".nii.gz"
computer_name = socket.gethostname()

if computer_name == 'samos':
    main_path = '/mnt/paros_MRI/jacques/'
elif 'santorini' in computer_name:
    main_path = '/Users/alex/jacques/'
else:
    raise Exception('No other computer name yet')

project = "AMD"

testing = False
if testing:
    main_path = '/Users/alex/jacques/AMD_testing_zone/AMD_TRK_testing'

if project == "AMD":
    path_TRK = os.path.join(main_path, 'AMD', 'TRK')
    path_DWI = os.path.join(main_path, 'AMD', 'DWI')
    path_transforms = os.path.join(main_path, 'AMD','Transforms')
    ref = "md"
    path_trk_tempdir = os.path.join(main_path, 'AMD', 'TRK_tempsave')
    path_trk_MDT = os.path.join(main_path, 'AMD', 'TRK_MDT')
    DWI_save = os.path.join(main_path, 'AMD', 'NII_tempsave')
    #Get the values from DTC_launcher_ADDecode. Should probably create a single parameter file for each project one day

stepsize = 2
ratio = 1
trkroi = ["wholebrain"]
str_identifier = get_str_identifier(stepsize, ratio, trkroi)
prune= True
overwrite = False
cleanup = False
verbose=True
recenter=1

save_temp_files = False
nii_test_files = False

contrast='dwi'
native_ref=''


if save_temp_files:
    mkcdir(path_trk_tempdir)

for subj in subjects:
    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    runno_to_MDT = os.path.join(path_transforms, f'{subj}_to_MDT_warp.nii.gz')

    if nii_test_files:
        mkcdir(DWI_save)
        SAMBA_preprocess = os.path.join(DWI_save, f'{subj}_{contrast}_masked{ext}')
        if recenter:
            SAMBA_preprocess_recentered_1 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_1{ext}')
            recenter_nii_save(SAMBA_preprocess, SAMBA_preprocess_recentered_1,verbose=True)
            SAMBA_preprocess_recentered_2 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_2{ext}')
            add_translation(SAMBA_preprocess_recentered_1, SAMBA_preprocess_recentered_2, translation = [0,0,-33], verbose=True)
            SAMBA_preprocess = SAMBA_preprocess_recentered_2
        SAMBA_preprocess_test_posttrans = os.path.join(DWI_save, f'{subj}_{contrast}_masked_posttrans{ext}')
        SAMBA_preprocess_test_rigid = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid{ext}')
        SAMBA_preprocess_test_rigid_affine = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid_affine{ext}')
        SAMBA_preprocess_test_postwarp = os.path.join(DWI_save, f'{subj}_{contrast}_postwarp{ext}')
        if native_ref == '':
            native_ref = SAMBA_preprocess
        if not os.path.exists(SAMBA_preprocess_test_postwarp) or overwrite:
            cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_preprocess} -r {native_ref}  -n Linear  -o {SAMBA_preprocess_test_posttrans} -t {trans}'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_posttrans} -o {SAMBA_preprocess_test_rigid} ' \
                f'-r {native_ref} -n Linear -t [{rigid},0]'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_rigid} -o {SAMBA_preprocess_test_rigid_affine} ' \
                f'-r {native_ref} -n Linear -t [{affine_orig},0]'
            os.system(cmd)

            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_rigid_affine} -o {SAMBA_preprocess_test_postwarp} ' \
                f'-r {native_ref} -n Linear -t {runno_to_MDT}'
            os.system(cmd)

    trk_preprocess_posttrans = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_posttrans.trk')
    trk_preprocess_postrigid = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_postrigid.trk')
    trk_preprocess_postrigid_affine = os.path.join(path_trk_tempdir, f'{subj}{str_identifier}_preprocess_postrigid_affine.trk')
    trk_MDT_space = os.path.join(path_trk_tempdir, f'{subj}_MDT.trk')

    if not os.path.exists(trk_MDT_space) or overwrite:
        orient_string = os.path.join(path_DWI, f'{subj}_relative_orientation.txt')
        if os.path.exists(orient_string):
            orient_relative = open(orient_string, mode='r').read()
        else:      
            warnings.warn(f'could not find {orient_string}')
            orient_string_varpath = glob.glob(os.path.join(path_DWI, '*relative_orientation.txt'))[0]  
            orient_relative = open(orient_string_varpath, mode='r').read()
        orientation_out = orient_relative.split(',')[0]
        orientation_out = orientation_out.split(':')[1]
        orientation_in = orient_relative.split(',')[1]
        orientation_in = orientation_in.split(':')[1]
        subj_dwi = os.path.join(path_DWI, f'{subj}_subjspace_dwi{ext}')

        SAMBA_input_real_file =  os.path.join(path_DWI, f'{subj}_dwi{ext}')

        #new_affine, translation, translate_affine = recenter_nii_affine(SAMBA_input_real_file, return_translation=True)

        #subj_trk, _ = gettrkpath(path_TRK, subj, str_identifier, pruned=True, verbose=verbose)

        check_dif_ratio(path_TRK, subj, str_identifier, ratio)
        subj_trk, trkexists = gettrkpath(path_TRK, subj, str_identifier, pruned=prune, verbose=False)

        _, exists = check_files([trans, rigid, runno_to_MDT])
        if np.any(exists==0):
            raise Exception('missing transform file')
        if not os.path.exists(affine) and not os.path.exists(affine_orig):
            raise Exception('missing transform file')
        streamlines_prepro, header = unload_trk(subj_trk)

        #if recenter:
            #new_affine, translation, translate_affine = recenter_nii_affine(SAMBA_preprocess,
        #                                                                    return_translation=True)
            #streamlines_prepro_recenter = transform_streamlines(streamlines_prepro, translate_affine)
            #save_trk_header(filepath= trk_preprocess_recentered_1, streamlines = streamlines_prepro_recenter, header = header,
            #        affine=np.eye(4), verbose=verbose)
            #translate_affine = np.eye(4)
            #translate_affine[:3,3] = [0,0,-33]
            #streamlines_prepro_recenter = transform_streamlines(streamlines_prepro_recenter, translate_affine)
            #save_trk_header(filepath= trk_preprocess_recentered_2, streamlines = streamlines_prepro_recenter, header = header,
            #        affine=np.eye(4), verbose=verbose)

        #streamlines_prepro, header_prepro = unload_trk(trk_preprocess)
        mat_struct = loadmat(trans)
        var_name = list(mat_struct.keys())[0]
        later_trans_mat = mat_struct[var_name]
        new_transmat = np.eye(4)
        vox_dim = [1,1,-1]
        new_transmat[:3,3] = np.squeeze(later_trans_mat[3:6]) * vox_dim

        #pix_dim = np.eye(4)
        #vox_dim = [1,1,2]
        #for i in np.arange(3):
        #    pix_dim[i,i]=vox_dim[i]
        #new_transmat = np.matmul(pix_dim, new_transmat)
        streamlines_posttrans = transform_streamlines(streamlines_prepro, (new_transmat))
        if (not os.path.exists(trk_preprocess_posttrans) or overwrite) and save_temp_files:
            save_trk_header(filepath= trk_preprocess_posttrans, streamlines = streamlines_posttrans, header = header,
                    affine=np.eye(4), verbose=verbose)

        rigid_struct = loadmat(rigid)
        var_name = list(rigid_struct.keys())[0]
        rigid_ants = rigid_struct[var_name]
        rigid_mat = convert_ants_vals_to_affine(rigid_ants)

        #streamlines_posttrans, header_posttrans = unload_trk(trk_preprocess_posttrans)
        streamlines_postrigid = transform_streamlines(streamlines_posttrans, np.linalg.inv(rigid_mat))

        if (not os.path.exists(trk_preprocess_postrigid) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_postrigid, streamlines=streamlines_postrigid, header=header,
                    affine=np.eye(4), verbose=verbose)

        overwrite=True
        if os.path.exists(affine):
            affine_mat_s = read_affine_txt(affine)
        else:
            cmd = f'ConvertTransformFile 3 {affine_orig} {affine} --matrix'
            os.system(cmd)
            affine_mat_s = read_affine_txt(affine)

        affine_mat = np.eye(4)
        affine_mat[:3, :3] = affine_mat_s
        #streamlines_postrigid, header_postrigid = unload_trk(trk_preprocess_postrigid)
        vox_dim = [1,1,-1]
        vox_mat = np.eye(4)
        for i in np.arange(3):
            vox_mat[i,i] = vox_dim[i]
        #affine_mat = np.matmul(vox_mat, affine_mat)
        streamlines_postrigidaffine = transform_streamlines(streamlines_postrigid, (affine_mat))

        if (not os.path.exists(trk_preprocess_postrigid_affine) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_preprocess_postrigid_affine, streamlines=streamlines_postrigidaffine, header=header,
                    affine=np.eye(4), verbose=verbose)

        #streamlines_postrigidaffine, header_postrigidaffine = unload_trk(trk_preprocess_postrigid_affine)

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
                    affine=np.eye(4), verbose=verbose)
    else:
        print(f'{trk_MDT_space} already exists')
