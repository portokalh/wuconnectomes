import os, glob

MDT_folder = '/Volumes/Data/Badea/Lab/human/AMD/ref/'
DWI_folder = '/Volumes/Data/Badea/Lab/human/AMD/DWI'

files = glob.glob(os.path.join(DWI_folder, '*'))

for ref_file in files:

    MDT_to_runno = ''
    SAMBA_preprocess_prewarp = ''

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {ref_file} -o {SAMBA_preprocess_prewarp} ' \
          f'-r {ref_file} -n Linear -t {MDT_to_runno}'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_prewarp} -o {SAMBA_preprocess_test_rigid_affine} ' \
          f'-r {SAMBA_preprocess_test_posttrans_2} -n Linear -t [{affine_orig},0]'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_test_posttrans_3} -o {SAMBA_preprocess_test_rigid} ' \
          f'-r {SAMBA_preprocess_test_posttrans_2} -n Linear -t [{rigid},0]'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_preprocess_test_posttrans_2} -r {SAMBA_preprocess_test_posttrans_2}  -n Linear  -o {SAMBA_preprocess_test_posttrans_3} -t {trans}'
    os.system(cmd)





mkcdir(DWI_save)
SAMBA_preprocess_ref = os.path.join(path_DWI, f'{subj}_labels{ext}')
SAMBA_coreg_ref = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
SAMBA_init = os.path.join(path_DWI, f'{subj}_{contrast}{ext}')
SAMBA_init = subj_dwi
SAMBA_preprocess = os.path.join(DWI_save, f'{subj}_{contrast}_preprocess{ext}')
# SAMBA_preprocess_2 = os.path.join(DWI_save, f'{subj}_{contrast}_preprocess_2{ext}')
# SAMBA_preprocess_recentered_1 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_1{ext}')
# SAMBA_preprocess_recentered_2 = os.path.join(DWI_save, f'{subj}_{contrast}_masked_recenter_2{ext}')

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