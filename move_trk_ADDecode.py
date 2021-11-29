import os
from transform_handler import get_affine_transform, get_affine_transform_test, get_flip_affine, header_superpose, recenter_nii_affine, \
    convert_ants_vals_to_affine, read_affine_txt, recenter_nii_save, add_translation, recenter_affine, \
    img_transform_exec, recenter_nii_save_test,recenter_nii_save_pure, recenter_affine_test
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

subjects = ["S02654", "S02666",  "S02670",  "S02686", "S02695", "S02720", "S02737", "S02753", "S02765", "S02781", "S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02926", "S02938", "S02939", "S02954", "S02967",
                "S02987", "S03010", "S03017", "S03033", "S03034", "S03045", "S03048","S01912", "S02110", "S02224", "S02227", "S02231", "S02266", "S02289", "S02320", "S02361", "S02363", "S02373", "S02386", "S02390", "S024S02", "S02410", "S02421", "S02424", "S02446", "S02451", "S02469", "S02473", "S02485", "S02491", "S02506"]
subjects = ["S03010", "S03017", "S03033", "S03034", "S03045", "S03048","S01912", "S02110", "S02224", "S02227", "S02231", "S02266", "S02289", "S02320", "S02361", "S02363", "S02373", "S02386", "S02390", "S024S02", "S02410", "S02421", "S02424", "S02446", "S02451", "S02469", "S02473", "S02485", "S02491", "S02506"]
subjects = ["S02320", "S02361", "S02363", "S02373", "S02386", "S02390", "S024S02", "S02410", "S02421", "S02424", "S02446", "S02451", "S02469", "S02473", "S02485", "S02491", "S02506"]
subjects= ["S02802", "S02813", "S02817", "S02840", "S02877", "S02898", "S02926", "S02938", "S02939", "S02954", "S02967",
                "S02987"]
subjects = ["S02802"]
#subjects = ["S02938", "S02939", "S02954", "S02967",
#                "S02987", "S03010", "S03017", "S03033", "S03034", "S03045", "S03048","S01912", "S02110", "S02224", "S02227", "S02231", "S02266", "S02289", "S02320", "S02361", "S02363", "S02373", "S02386", "S02390", "S024S02", "S02410", "S02421", "S02424", "S02446", "S02451", "S02469", "S02473", "S02485", "S02491", "S02506"]
# "S02715", S02690, S027701. S02926, S02804
ext = ".nii.gz"
computer_name = socket.gethostname()

if computer_name == 'samos':
	main_path = '/Volumes/Data/Badea/'
elif 'santorini' in computer_name:
	main_path = '/Volumes/Data/Badea/Lab/human/'
else:
	raise Exception('No other computer name yet')

project = "AD_Decode"

if project == "AD_Decode":
	path_TRK = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_MPCA_100')
	path_TRK_output = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_MPCA_MDT_100')
	path_DWI = os.path.join(main_path, 'AD_Decode', 'Analysis', 'DWI')
	path_transforms = os.path.join(main_path, 'AD_Decode', 'Analysis','Transforms')
	ref = "md"
	path_trk_tempdir = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_save_alt2')
	DWI_save = os.path.join(main_path, 'AD_Decode', 'Analysis', 'NII_tempsave_alt2')

	mkcdir([path_trk_tempdir,path_TRK_output, DWI_save])
	#Get the values from DTC_launcher_ADDecode. Should probalby create a single parameter file for each project one day
	stepsize = 2
	ratio = 100
	trkroi = ["wholebrain"]
	str_identifier = get_str_identifier(stepsize, ratio, trkroi)

overwrite = True
cleanup = False
verbose=True
save_temp_files = True
recenter=1
contrast='dwi'
native_ref=''

orient_string = os.path.join(path_DWI, 'relative_orientation.txt')
orient_relative = open(orient_string, mode='r').read()
orientation_out = orient_relative.split(',')[0]
orientation_out = orientation_out.split(':')[1]
orientation_in = orient_relative.split(',')[1]
orientation_in = orientation_in.split(':')[1]

nii_test_files = 0

for subj in subjects:
    subj_trk, _ = gettrkpath(path_TRK, subj, str_identifier, pruned=True, verbose=verbose)
    subj_dwi = os.path.join(path_DWI, f'{subj}_subjspace_dwi{ext}')

    #if not os.path.exists(subj_trk):
    #	print(f'could not find {subj_trk}, skipping')
    #continue
    trkname = os.path.basename(subj_trk)
    trk_MDT_space = os.path.join(path_TRK_output, trkname)

    #preprocess_target = os.path.join(path_DWI,{subj}+'_dwi_masked.nii.gz')
    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    runno_to_MDT = os.path.join(path_transforms, f'{subj}_to_MDT_warp.nii.gz')

    print(f'Beginning the process to transfer trk file {subj_trk} to {trk_MDT_space}')
    if nii_test_files:
        mkcdir(DWI_save)
        SAMBA_preprocess_ref = os.path.join(path_DWI, f'{subj}_reference{ext}')
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
        SAMBA_preprocess_test_rigid = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid{ext}')
        SAMBA_preprocess_test_rigid_affine = os.path.join(DWI_save, f'{subj}_{contrast}_postrigid_affine{ext}')
        SAMBA_preprocess_test_postwarp = os.path.join(DWI_save, f'{subj}_{contrast}_postwarp{ext}')
        if native_ref == '':
            native_ref = SAMBA_init
        if not os.path.exists(SAMBA_preprocess_test_postwarp) or overwrite:
            cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_init} -r {SAMBA_init}  -n Linear  -o {SAMBA_preprocess_test_posttrans} -t {trans}'
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

    overwrite = True
    if not os.path.exists(trk_MDT_space) or overwrite:
        reference = os.path.join(path_DWI, f'{subj}_reference{ext}')
        SAMBA_input_real_file =  os.path.join(path_DWI, f'{subj}_dwi{ext}')

        subj_trk, _ = gettrkpath(path_TRK, subj, str_identifier, pruned=True, verbose=verbose)
        if not os.path.exists(subj_trk):
            print(f'could not find {subj_trk}, skipping')
            continue

        nii = nib.load(subj_dwi)
        nii_data = nii.get_data()
        subj_affine = nii._affine
        subjtorecenter_affine = recenter_affine(np.shape(nii_data), subj_affine)
        subjtorecenter_affine = recenter_affine_test(np.shape(nii_data), subj_affine)
        subj_torecenter_transform_affine = get_affine_transform_test(subj_affine, subjtorecenter_affine)
        #subj_torecenter_transform_affine_inverse = np.eye(4)
        #subj_torecenter_transform_affine_inverse[:3,:3] = np.linalg.inv(subj_torecenter_transform_affine[:3,:3])
        #subj_torecenter_transform_affine_inverse[:3,3] = -subj_torecenter_transform_affine[:3,3]
        #subj_torecenter_transform_affine_inverse[:3,3] = [-4.5, -28.98, -23.92]

        subj_affine_new = np.eye(4)
        subj_affine_new[0,0]=-1
        subj_affine_new[1,1]=-1
        subj_affine_new[:,3]=subj_affine[:,3]
        subj_dwi_new = nib.Nifti1Image(nii_data,subj_affine_new)
        subj_dwi_new_path = os.path.join(DWI_save, f'{subj}_subjspace_dwi_eye{ext}')
        nib.save(subj_dwi_new, subj_dwi_new_path)
        subj_affine_new_trk = get_affine_transform_test(subj_affine, subj_affine_new)
        reorient_trans = np.matmul(subj_torecenter_transform_affine[:3,:3],np.shape(nii_data))-np.shape(nii_data)
        reorient_trans = np.array([x / 2 for x in reorient_trans])
        reorient_trans = np.multiply(reorient_trans, [1, 1, -1])
        subj_affine_new_trk[:3,3] = reorient_trans

        preprocess_target = f'/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/preprocess/{subj}_dwi_masked.nii.gz'
        nii_preprocess = nib.load(preprocess_target)
        preprocess_affine = nii_preprocess._affine

        subj_torecenter_transform_affine = get_affine_transform_test(subj_affine, subj_affine_new)
        added_trans = subj_affine[:3, 3] + np.multiply(preprocess_affine[:3, 3], [1,1,-1]) + [-1,-1,0]
        added_trans = subj_affine[:3, 3] + np.multiply(subjtorecenter_affine[:3, 3], [-1,-1,-1])
        subj_torecenter_transform_affine[:3, 3] = reorient_trans + added_trans

        #subj_torecenter_transform_affine_test = get_affine_transform_test(subj_affine, preprocess_affine)

        trk_recentered = os.path.join(path_trk_tempdir, f'{subj}_recentered.trk')
        trk_posttrans = os.path.join(path_trk_tempdir, f'{subj}_preprocess_posttrans.trk')
        trk_postrigid = os.path.join(path_trk_tempdir, f'{subj}_preprocess_postrigid.trk')
        trk_postrigid_affine = os.path.join(path_trk_tempdir, f'{subj}_preprocess_postrigid_affine.trk')

        _, exists = check_files([trans, rigid, affine, runno_to_MDT])
        if np.any(exists==0):
            raise Exception('missing transform file')

        streamlines, header = unload_trk(subj_trk)
        #subj_tocenter_streamlines = transform_streamlines(streamlines, np.linalg.inv(subj_torecenter_transform_affine), in_place=False)

        print(subj_affine_new_trk)
        subj_eye_streamlines = transform_streamlines(streamlines, np.linalg.inv(subj_affine_new_trk),
                                                     in_place=False)

        trk_subjspace_eye = os.path.join(path_trk_tempdir, f'{subj}_subjspace_eye.trk')
        if (not os.path.exists(trk_subjspace_eye) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_subjspace_eye, streamlines=subj_eye_streamlines, header=header,
                            affine=np.eye(4), verbose=verbose)

        print(subj_torecenter_transform_affine)
        subj_tocenter_streamlines = transform_streamlines(streamlines, np.linalg.inv(subj_torecenter_transform_affine),
                                                          in_place=False)
        if (not os.path.exists(trk_recentered) or overwrite) and save_temp_files:
            save_trk_header(filepath= trk_recentered, streamlines = subj_tocenter_streamlines, header = header, affine=np.eye(4), verbose=verbose)

        mat_struct = loadmat(trans)
        var_name = list(mat_struct.keys())[0]
        later_trans_mat = mat_struct[var_name]
        new_transmat = np.eye(4)
        vox_dim = [1, 1, -1]
        new_transmat[:3, 3] = np.squeeze(later_trans_mat[3:6]) * vox_dim
        #new_transmat[:3, 3] = np.squeeze(np.matmul(subj_torecenter_transform_affine[:3,:3],later_trans_mat[3:6]))
        #new_transmat[:3, 3] = np.squeeze(np.matmul(subj_torecenter_transform_affine[:3, :3], later_trans_mat[3:6])) + [-1,-1,0]
        new_transmat[:3, 3] = np.squeeze(np.matmul(subj_torecenter_transform_affine[:3, :3], later_trans_mat[3:6]))
        new_transmat[2, 3] = 0
        # pix_dim = np.eye(4)
        # vox_dim = [1,1,2]
        # for i in np.arange(3):
        #    pix_dim[i,i]=vox_dim[i]
        # new_transmat = np.matmul(pix_dim, new_transmat)
        print(new_transmat)
        streamlines_posttrans = transform_streamlines(subj_tocenter_streamlines, (new_transmat))
        if (not os.path.exists(trk_posttrans) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_posttrans, streamlines=streamlines_posttrans, header=header,
                            affine=np.eye(4), verbose=verbose)

        rigid_struct = loadmat(rigid)
        var_name = list(rigid_struct.keys())[0]
        rigid_ants = rigid_struct[var_name]
        rigid_mat = convert_ants_vals_to_affine(rigid_ants)

        # streamlines_posttrans, header_posttrans = unload_trk(trk_preprocess_posttrans)
        print(rigid_mat)
        streamlines_postrigid = transform_streamlines(streamlines_posttrans, np.linalg.inv(rigid_mat))

        if (not os.path.exists(trk_postrigid) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_postrigid, streamlines=streamlines_postrigid, header=header,
                            affine=np.eye(4), verbose=verbose)

        overwrite = True

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

        # streamlines_postrigid, header_postrigid = unload_trk(trk_preprocess_postrigid)
        print(affine_mat)
        streamlines_postrigidaffine = transform_streamlines(streamlines_postrigid, np.linalg.inv(affine_mat))

        if (not os.path.exists(trk_postrigid_affine) or overwrite) and save_temp_files:
            save_trk_header(filepath=trk_postrigid_affine, streamlines=streamlines_postrigidaffine,
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
                            affine=np.eye(4), verbose=verbose)
    else:
        print(f'{trk_MDT_space} already exists')
