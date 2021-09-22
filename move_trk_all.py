import os
from transform_handler import get_affine_transform, get_flip_affine, header_superpose, recenter_nii_affine, \
    convert_ants_vals_to_affine, read_affine_txt
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

subj='S02654'
#subj = 'S01912'
ext = ".nii.gz"
contrast='dwi'

computer_name = socket.gethostname()

if computer_name == 'samos':
    main_path = '/mnt/paros_MRI/jacques/'
else:
    raise Exception('No other computer name yet')

project = "AD_Decode"

if project == "AD_Decode":
    path_TRK = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK')
    path_DWI = os.path.join(main_path, 'AD_Decode', 'Analysis', 'DWI')
    ref = "md"
    save_trk_tempdir = os.path.join(main_path, 'AD_Decode', 'Analysis', 'TRK_save')
    mkcdir(save_trk_tempdir)

if ref=="md":
    reference= os.path.join(path_DWI},f'{subj}_tmp_md{ext}')
elif ref=="coreg":
    reference= os.path.join(path_DWI,f'Reg_{D_subj}_nii4D{ext}')
else:
    raise Exception('No other ref yet')

gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
overwrite = False
cleanup = False
verbose=True
recenter=0
skip = False

labpath = '/Volumes/Data/Badea/Lab/'
work_dir = os.path.join(labpath, f'mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work')

orient_string = os.path.join(path_DWI, 'relative_orientation.txt')
orient_relative = open(orient_string, mode='r').read()
orientation_out = orient_relative.split(',')[0]
orientation_out = orientation_out.split(':')[1]
orientation_in = orient_relative.split(',')[1]
orientation_in = orientation_in.split(':')[1]

subj_dwi = os.path.join(test_dir, f'{subj}_tmp_{contrast}{ext}')
SAMBA_input_real_file = os.path.join(labpath,f'mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-inputs/{subj}_{contrast}.nii.gz')

subj_trk = os.path.join(test_dir,f'{subj}_stepsize_2_ratio_100_wholebrain_pruned.trk')
trk_filepath_tmp2 = os.path.join(test_dir, f'{subj}_ratio_100_tmp2_notdirect.trk')
inputS_trk_new = os.path.join(test_dir, f'{subj}_ratio_100_final_direct.trk')
trk_preprocess = os.path.join(test_dir, f'{subj}_preprocess_direct.trk')
trk_preprocess_posttrans = os.path.join(test_dir, f'{subj}_{contrast}_preprocess_posttrans_direct.trk')
trk_preprocess_postrigid = os.path.join(test_dir, f'{subj}_{contrast}_preprocess_postrigid_direct.trk')
trk_preprocess_postrigid_affine = os.path.join(test_dir, f'{subj}_{contrast}_preprocess_postrigid_affine_direct.trk')
trk_preprocess_postwarp = os.path.join(test_dir, f'{subj}_{contrast}_postwarp_direct_2.trk')

trans = os.path.join(work_dir, "preprocess", "base_images", "translation_xforms", f"{subj}_0DerivedInitialMovingTranslation.mat")
rigid = os.path.join(work_dir, "dwi", f"{subj}_rigid.mat")
affine = os.path.join(work_dir, "dwi", f"{subj}_affine.mat")
runno_to_MDT = os.path.join(work_dir, f'dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/reg_diffeo/{subj}_to_MDT_warp.nii.gz')

save_temp_files = False
overwrite = False

affine_map_test = get_affine_transform(reference, subj_dwi)
streamlines, header = unload_trk(subj_trk)

tmp2_streamlines = transform_streamlines(streamlines, np.linalg.inv(affine_map_test), in_place=False)

if (not os.path.exists(trk_filepath_tmp2) or overwrite) and save_temp_files:
    save_trk_header(filepath= trk_filepath_tmp2, streamlines = tmp2_streamlines, header = header, affine=np.eye(4), verbose=verbose)

#tmp2_streamlines, header_new = unload_trk(trk_filepath_tmp2)
affine_transform, newaffine = get_flip_affine(orientation_in, orientation_out)
affine_transform_new = affine_transform

#affine_transform[:3,3] = affine_transform.diagonal()[:3]* header[0][:3,3] - header[0][:3,3]

center1 = header[0][:3, 3]
center2 = affine_transform.diagonal()[:3] * center1[:3]
affine_transform_new[:3, 3] = center1 - center2
affine_transform_new[1, 3] = affine_transform[0, 3]

#affine_transform_new = np.eye(4)
#affine_transform_new[:3,:3] = np.dot(affine_transform[:3,:3], affine_map_test[:3,:3])
project_input_streamlines = transform_streamlines(tmp2_streamlines, np.linalg.inv(affine_transform_new))

if (not os.path.exists(inputS_trk_new) or overwrite) and save_temp_files:
    save_trk_header(filepath=inputS_trk_new, streamlines=project_input_streamlines, header=header, affine=np.eye(4),
                                   verbose=verbose)

new_affine, translation, translate_affine = recenter_nii_affine(SAMBA_input_real_file, return_translation=True)
#project_input_streamlines, _ = unload_trk(inputS_trk_new)
streamlines_prepro = transform_streamlines(project_input_streamlines, translate_affine)

if (not os.path.exists(trk_preprocess) or overwrite) and save_temp_files:
    save_trk_header(filepath= trk_preprocess, streamlines = streamlines_prepro, header = header,
                    affine=np.eye(4), verbose=verbose)

#streamlines_prepro, header_prepro = unload_trk(trk_preprocess)
mat_struct = loadmat(trans)
var_name = list(mat_struct.keys())[0]
later_trans_mat = mat_struct[var_name]
new_transmat = np.eye(4)
new_transmat[:3,3] = np.squeeze(later_trans_mat[3:6])
streamlines_posttrans = transform_streamlines(streamlines_prepro, new_transmat)

if (not os.path.exists(trk_preprocess_posttrans) or overwrite) and save_temp_files:
    save_trk_header(filepath= trk_preprocess_posttrans, streamlines = streamlines_posttrans, header = header,
                    affine=np.eye(4), verbose=verbose)

rigid_struct = loadmat(rigid)
var_name = list(rigid_struct.keys())[0]
rigid_ants = rigid_struct[var_name]
rigid_mat = convert_ants_vals_to_affine(rigid_ants)

#streamlines_posttrans, header_posttrans = unload_trk(trk_preprocess_posttrans)
streamlines_postrigid = transform_streamlines(streamlines_posttrans, np.linalg.inv(rigid_mat))

if not os.path.exists(trk_preprocess_postrigid) or overwrite:
    save_trk_header(filepath=trk_preprocess_postrigid, streamlines=streamlines_postrigid, header=header,
                    affine=np.eye(4), verbose=verbose)

affine_mat_path = os.path.join(test_dir, 'new_affine.txt')
cmd = f'ConvertTransformFile 3 {affine} {affine_mat_path} --matrix'
os.system(cmd)
affine_mat_s = read_affine_txt(affine_mat_path)
affine_mat = np.eye(4)
affine_mat[:3, :3] = affine_mat_s
#streamlines_postrigid, header_postrigid = unload_trk(trk_preprocess_postrigid)
streamlines_postrigidaffine = transform_streamlines(streamlines_postrigid, np.linalg.inv(affine_mat))

if not os.path.exists(trk_preprocess_postrigid_affine) or overwrite and save_temp_files:
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

if (not os.path.exists(trk_preprocess_postwarp) or overwrite):
    save_trk_header(filepath=trk_preprocess_postwarp, streamlines=streamlines_post_warp, header=header,
                    affine=np.eye(4), verbose=verbose)