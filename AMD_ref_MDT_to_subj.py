import os, glob
from file_tools import mkcdir
from transform_handler import recenter_nii_save_test, affine_superpose, get_affine_transform_test

mainpath = '/Volumes/Data/Badea/Lab/human/AMD/'
ref_subj_folder = os.path.join(mainpath,'ref_subj')
ref_MDT_folder = os.path.join(mainpath,'ref_MDT')
#DWI_folder = os.path.join(mainpath,'MDT_to_subj_testzone')
path_transforms = os.path.join(mainpath,'Transforms')
MDT_to_subj = os.path.join(mainpath,'MDT_to_subj_testzone')
subj_to_MDT = os.path.join(mainpath,'subj_to_MDT_testzone')
mkcdir([MDT_to_subj,subj_to_MDT])
ext='.nii.gz'


files = glob.glob(os.path.join(ref_subj_folder, '*'))
files = ['/Volumes/Data/Badea/Lab/human/AMD/ref_subj/H26578_subjspace_fa.nii.gz']
recenter = True
overwrite = True

for ref_file in files:
    subj = os.path.basename(ref_file)[:6]
    ref = os.path.basename(ref_file)[17:19]
    SAMBA_init = ref_file
    print(ref_file)
    preprocess = os.path.join(subj_to_MDT, f'{subj}_prepro.nii.gz')
    if recenter and (not os.path.exists(preprocess) or overwrite):
        recenter_nii_save_test(SAMBA_init, preprocess)
        SAMBA_init = preprocess

    ref_posttrans = os.path.join(subj_to_MDT, f'{subj}_{ref}_masked_posttrans{ext}')
    ref_posttrans_2 = os.path.join(subj_to_MDT, f'{subj}_{ref}_masked_posttrans_2{ext}')
    ref_posttrans_3 = os.path.join(subj_to_MDT, f'{subj}_{ref}_masked_posttrans_3{ext}')

    ref_rigid = os.path.join(subj_to_MDT, f'{subj}_{ref}_postrigid{ext}')
    ref_rigid_affine = os.path.join(subj_to_MDT, f'{subj}_{ref}_postrigid_affine{ext}')
    ref_postwarp = os.path.join(subj_to_MDT, f'{subj}_{ref}_postwarp{ext}')
    runno_to_MDT_warp = os.path.join(path_transforms, f'{subj}_to_MDT_warp.nii.gz')

    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    MDT_to_runno_warp = os.path.join(path_transforms,f'{subj}_to_MDT_warp.nii.gz')


    if not os.path.exists(ref_postwarp) or overwrite:
        #if not os.path.exists(ref_posttrans):
        #    cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_init} -r {SAMBA_init}  -n Linear  -o {ref_posttrans}'
        #    os.system(cmd)
        #if not os.path.exists(ref_posttrans_2) or overwrite:
        #    affine_superpose(SAMBA_init, ref_posttrans, outpath=ref_posttrans_2)

        if not os.path.exists(ref_posttrans) or overwrite:
            cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_init} -r {SAMBA_init}  -n Linear  -o {ref_posttrans} -t {trans}'
            os.system(cmd)

        if not os.path.exists(ref_rigid) or overwrite:
            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {ref_posttrans} -o {ref_rigid} ' \
                  f'-r {SAMBA_init} -n Linear -t [{rigid},0]'
            os.system(cmd)

        if not os.path.exists(ref_rigid_affine) or overwrite:
            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {ref_rigid} -o {ref_rigid_affine} ' \
                  f'-r {SAMBA_init} -n Linear -t [{affine_orig},0]'
            os.system(cmd)

        if not os.path.exists(ref_postwarp) or overwrite:
            cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {ref_rigid_affine} -o {ref_postwarp} ' \
                  f'-r {SAMBA_init} -n Linear -t {runno_to_MDT_warp}'
        print(cmd)
        os.system(cmd)


files = glob.glob(os.path.join(ref_MDT_folder, '*'))
files = ['/Volumes/Data/Badea/Lab/human/AMD/ref_MDT/H26578_fa_to_MDT.nii.gz']
ext = '.nii.gz'

for ref_file in files:
    subj = os.path.basename(ref_file)[:6]
    ref = os.path.basename(ref_file)[7:9]
    trans = os.path.join(path_transforms, f"{subj}_0DerivedInitialMovingTranslation.mat")
    rigid = os.path.join(path_transforms, f"{subj}_rigid.mat")
    affine_orig = os.path.join(path_transforms, f"{subj}_affine.mat")
    affine = os.path.join(path_transforms, f"{subj}_affine.txt")
    MDT_to_runno_warp = os.path.join(path_transforms,f'{subj}_to_MDT_warp.nii.gz')

    SAMBA_preprocess_prewarp = os.path.join(MDT_to_subj,f'{subj}_{ref}_prewarp.nii.gz')
    SAMBA_preprocess_preaffine = os.path.join(MDT_to_subj,f'{subj}_{ref}_preaffine.nii.gz')
    SAMBA_preprocess_prerigid = os.path.join(MDT_to_subj,f'{subj}_{ref}_prerigid.nii.gz')
    SAMBA_preprocess_pretrans = os.path.join(MDT_to_subj,f'{subj}_{ref}_pretrans.nii.gz')

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {ref_file} -o {SAMBA_preprocess_prewarp} ' \
          f'-r {ref_file} -n Linear -t {MDT_to_runno_warp}'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_prewarp} -o {SAMBA_preprocess_preaffine} ' \
          f'-r {ref_file} -n Linear -t [{affine_orig},1]'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 --float -d 3 -i {SAMBA_preprocess_preaffine} -o {SAMBA_preprocess_prerigid} ' \
          f'-r {ref_file} -n Linear -t [{rigid},1]'
    os.system(cmd)

    cmd = f'antsApplyTransforms -v 1 -d 3  -i {SAMBA_preprocess_prerigid} -r {ref_file}  -n Linear  -o {SAMBA_preprocess_pretrans} -t [{trans},1]'
    os.system(cmd)

