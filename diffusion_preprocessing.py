
import os
import shutil
import subprocess
from file_tools import largerfile, mkcdir, getext, buildlink
from transform_handler import img_transform_exec, space_transpose, affine_superpose, header_superpose
import glob
from basic_LPCA_denoise import basic_LPCA_denoise_func
from mask_handler import applymask_samespace, median_mask_make
import numpy as np
import warnings
import time

def launch_preprocessing(subj, raw_nii, outpath, cleanup=False, nominal_bval=4000, SAMBA_inputs_folder=None,
                         shortcuts_all_folder = None, gunniespath="~/gunnies/", processes=1, masking="bet", ref=None,
                         transpose=None, overwrite=False, denoise='None', recenter=0, verbose=False):

    proc_name ="diffusion_prep_" # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline
    work_dir=os.path.join(outpath,proc_name+subj)
    """
    for filePath in glob.glob(os.path.join(work_dir,'*')):
        modTimesinceEpoc = os.path.getmtime(filePath)
        modificationTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(modTimesinceEpoc))
        if modificationTime[5:7]=='09' and int(modificationTime[8:10])>9:
            os.remove(filePath)
    """
    if verbose:
        print(f"Processing diffusion data with runno/subj: {subj}")
        print(f"Work directory is {work_dir}")
    mkcdir(work_dir)

    sbatch_folder =os.path.join(work_dir,"sbatch")
    mkcdir(sbatch_folder)
    #nii_path = os.path.join(work_dir,'nii4D_'+subj + '.nii.gz')
    #if not os.path.exists(nii_path):
    #    shutil.copy(raw_nii, nii_path)
    nii_name=os.path.split(raw_nii)[1]
    niifolder = os.path.dirname(raw_nii)
    ext = ".nii.gz"
    nii_ext=getext(nii_name)
    bxheader = nii_name.replace(nii_ext,".bxh")
    bxheader = os.path.join(niifolder, bxheader)
    bvecs = os.path.join(work_dir, subj+"_bvecs.txt")
    bvals =bvecs.replace("bvecs","bvals")

    if verbose:
        print(f"Original nifti is at {nii_name}\nbvecs are at {bvecs}\nbvals are at {bvals}\n")
    if not os.path.exists(bvecs):
        if verbose:
            print("Extracting diff directions")
        #print("Bvals and bvecs not found, using extractdiffdirs, however it it NOT RELIABLE, beware!")
        bvec_cmd = (f"extractdiffdirs --colvectors --writebvals --fieldsep='\t' --space=RAI {bxheader} {bvecs} {bvals}")
        os.system(bvec_cmd)

    # Make dwi for mask generation purposes.
    tmp_mask = os.path.join(work_dir,f"{subj}_tmp_mask{ext}")
    raw_dwi = os.path.join(work_dir,f"{subj}_raw_dwi.nii.gz")
    b0_dwi = os.path.join(work_dir,f"{subj}_b0_dwi.nii.gz")  #test average of the b0 images to make a better mask
    orient_string = os.path.join(work_dir,"relative_orientation.txt")

    if shortcuts_all_folder is not None:
        #nii_path_link = os.path.join(shortcuts_all_folder, f"{subj}_rawnii{ext}")
        #if not os.path.exists(nii_path_link) or overwrite:
        #    buildlink(nii_path, nii_path_link)
        bvecs_new = os.path.join(shortcuts_all_folder, subj + "_bvecs.txt")
        bvals_new = os.path.join(shortcuts_all_folder, subj + "_bvals.txt")
        if not os.path.exists(bvecs_new) or not os.path.exists(bvals_new) or overwrite:
            shutil.copyfile(bvecs, bvecs_new)
            shutil.copyfile(bvals, bvals_new)

    final_mask = os.path.join(work_dir, f'{subj}_mask{ext}')

    #if (not os.path.exists(final_mask) and not os.path.exists(tmp_mask)) or overwrite:
    if not os.path.exists(tmp_mask) or overwrite:
        if not os.path.exists(raw_dwi) or overwrite:
            select_cmd = f"select_dwi_vols {raw_nii} {bvals} {raw_dwi} {nominal_bval} -m"
            os.system(select_cmd)
        if not os.path.exists(b0_dwi) or overwrite:
            select_cmd = f"select_dwi_vols {raw_nii} {bvals} {b0_dwi} 0 -m"
            os.system(select_cmd)
        if not os.path.exists(tmp_mask) or overwrite:
            if 'median' in masking:
                tmp = tmp_mask.replace("_mask", "")
                if np.size(masking.split('_'))>1:
                    median_radius = int(masking.split('_')[1])
                else:
                    median_radius = 4
                median_mask_make(b0_dwi, tmp, outpathmask=tmp_mask, median_radius = median_radius, numpass=median_radius)
                #median_mask_make(b0_dwi, tmp, outpathmask='/Users/jas/jacques/Chavez_test_temp/b0_test.nii.gz', median_radius = median_radius, numpass=median_radius)
                #median_mask_make(raw_dwi, tmp, outpathmask='/Users/jas/jacques/Chavez_test_temp/007_mask_rad7.nii.gz',
                #                 median_radius=7, numpass=7)
            elif masking=="bet":
                tmp=tmp_mask.replace("_mask", "")
                bet_cmd = f"bet {raw_dwi} {tmp} -m -n -R"
                os.system(bet_cmd)
            else:
                raise Exception("Unrecognized masking type")

    # I think this part is done later more properly:     if create_subj_space_files: for contrast in ['dwi', 'b0', 'mask']:
    #if SAMBA_inputs_folder is not None:
    #    mask_subj_link = os.path.join(SAMBA_inputs_folder,f'{subj}_subjspace_mask{ext}')
    #    if not os.path.exists(mask_subj_link) or overwrite:
    #        shutil.copy(tmp_mask, mask_subj_link)

    #if cleanup and (os.path.exists(tmp_mask) and os.path.exists(raw_dwi)):
    #    os.remove(raw_dwi)
    #overwrite=False
    # Run Local PCA Denoising algorithm on 4D nifti:
    masked_nii = os.path.join(work_dir, nii_name)
    if not "nii.gz" in masked_nii:
        masked_nii = masked_nii.replace(".nii", ".nii.gz")
    masked_nii = masked_nii.replace(ext, "_masked" + ext)

    if denoise.lower()=='lpca':
        D_subj=f'LPCA_{subj}';
    elif denoise.lower()=='mpca':
        D_subj=f'MPCA_{subj}';
    elif denoise=="None" or denoise is None:
        D_subj = f'{subj}'

    if denoise=="None" or denoise is None:
        denoised_nii = masked_nii
        if not os.path.exists(masked_nii) or overwrite:
            fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
            os.system(fsl_cmd)
    else:
        denoised_nii = os.path.join(work_dir,f"{D_subj}_nii4D.nii.gz")
        if not os.path.exists(denoised_nii) or overwrite:
            if not os.path.exists(masked_nii) or overwrite:
                fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
                os.system(fsl_cmd)
            basic_LPCA_denoise_func(subj,masked_nii,bvecs,denoised_nii, processes=processes,
                                    denoise=denoise, verbose=False) #to improve and make multiprocessing

    #if cleanup and os.path.exists(denoised_nii) and os.path.exists(masked_nii) and denoised_nii!=masked_nii:
    #    os.remove(masked_nii)

    # Run coregistration/eddy current correction:

    coreg_nii_old = f'{outpath}/co_reg_{D_subj}_m00-results/Reg_{D_subj}_nii4D{ext}';
    coreg_nii = os.path.join(work_dir,f'Reg_{D_subj}_nii4D{ext}')
    if not cleanup:
        coreg_nii=coreg_nii_old
    if not os.path.exists(coreg_nii) or overwrite:
        if not os.path.exists(coreg_nii_old) or overwrite:
            temp_cmd = os.path.join(gunniespath,'co_reg_4d_stack_tmpnew.bash')+f' {denoised_nii} {D_subj} 0 {outpath} 0';
            os.system(temp_cmd)
        if cleanup:
            shutil.move(coreg_nii_old,coreg_nii)

    if shortcuts_all_folder is not None:
        coreg_link = os.path.join(shortcuts_all_folder,f'{subj}_subjspace_coreg{ext}')
        if not os.path.exists(coreg_link) or overwrite:
            buildlink(coreg_nii, coreg_link)

    toeddy=False
    if toeddy:
        #fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
        #os.system(fsl_cmd)
        eddy_cmd = f"eddy --imain={coreg_nii} --mask={tmp_mask} --acqp=acq_params.txt --index={os.path.join(work_dir,'index.txt')} --bvecs={bvecs} --bvals={bvals} --topup=topup_results --repol --out = {os.path.join(work_dir,f'Reg_{D_subj}_nii4D_eddy{ext}')}"
        os.system(eddy_cmd)

    coreg_inputs=os.path.join(outpath,f'co_reg_{D_subj}_m00-inputs')
    coreg_work=coreg_inputs.replace('-inputs','-work')
    coreg_results=coreg_inputs.replace('-inputs','-results')
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_inputs):
        shutil.rmtree(coreg_inputs)
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_work):
        shutil.rmtree(coreg_work)
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_results):
        shutil.rmtree(coreg_results)

    # Generate tmp DWI:

    tmp_dwi_out=os.path.join(work_dir, f'{subj}_tmp_dwi{ext}')
    dwi_out=os.path.join(work_dir,f'{subj}_dwi{ext}')

    if not os.path.exists(tmp_dwi_out) or overwrite:
        cmd=f'select_dwi_vols {coreg_nii} {bvals} {tmp_dwi_out} {nominal_bval}  -m'
        os.system(cmd)

    # Generate tmp B0:
    tmp_b0_out=os.path.join(work_dir,f'{subj}_tmp_b0{ext}')
    b0_out = os.path.join(work_dir, f'{subj}_b0{ext}')
    if (not os.path.exists(b0_out) and not os.path.exists(tmp_b0_out)) or overwrite:
        cmd=f'select_dwi_vols {coreg_nii} {bvals} {tmp_b0_out} 0  -m;'
        os.system(cmd)
    #overwrite=False
    #elif cleanup and os.path.exists(tmp_b0_out):
    #    os.remove(tmp_b0_out)

    # Generate DTI contrasts and perform some tracking QA:
    if cleanup:
        c_string=' --cleanup '
    else:
        c_string=''

    #Important note: this is what first creates the fa, md, etc
    if len(glob.glob(os.path.join(work_dir,f'*.fib.gz.md{ext}'))) == 0 or overwrite:
        if overwrite:
            oldfiles = glob.glob(os.path.join(work_dir, f'*.fib.gz*'))
            for oldfile in oldfiles:
                os.remove(oldfile)
        cmd = 'bash ' + os.path.join(gunniespath,'dti_qa_with_dsi_studio_weirdcall.bash')+f' {coreg_nii} {bvecs} {tmp_mask} {work_dir} {c_string}';
        os.system(cmd)

    #Save the subject space dti results

    #Generate tmp MD:
    for contrast in ['md']:
        real_file=largerfile(os.path.join(work_dir,f'*.fib.gz.{contrast}{ext}'))  #Catch the 'real file' for each contrast
        tmp_file = f'{work_dir}/{subj}_tmp_{contrast}{ext}';
        if not os.path.exists(tmp_file):
            shutil.copy(real_file,tmp_file)

    tmp_md = f'{work_dir}/{subj}_tmp_md{ext}';

    if ref=="md" or ref is None:
        reference=tmp_md
    elif ref=="coreg":
        reference=coreg_nii
    elif os.path.exists(ref):
        reference=ref

    reference_file = os.path.join(work_dir, f'{subj}_reference{ext}')
    if not os.path.exists(reference_file):
        shutil.copy(reference, reference_file)

    if shortcuts_all_folder is not None:
        bonus_ref_link = os.path.join(shortcuts_all_folder, f'{subj}_reference{ext}')
        if not os.path.exists(bonus_ref_link) or overwrite:
            buildlink(reference_file,bonus_ref_link)

    #give new header to the non-dti files using md as reference


    for contrast in ['dwi', 'b0', 'mask']:
        tmp_file=os.path.join(work_dir,f'{subj}_tmp_{contrast}{ext}')
        tmp2_file=os.path.join(work_dir,f'{subj}_tmp2_{contrast}{ext}')
        final_file=os.path.join(work_dir,f'{subj}_{contrast}{ext}')
        if ((not os.path.exists(tmp2_file) and not os.path.exists(final_file)) or overwrite):
            if not os.path.exists(tmp_file):
                raise Exception("Tmp file was not created, need to rerun previous processes")
            else:
                header_superpose(reference, tmp_file, outpath=tmp2_file)

    create_subj_space_files = True
    if create_subj_space_files:
        for contrast in ['dwi', 'b0', 'mask']:
            tmp_file = os.path.join(work_dir, f'{subj}_tmp_{contrast}{ext}')
            subj_file = os.path.join(work_dir, f'{subj}_subjspace_{contrast}{ext}')
            if not os.path.exists(subj_file) or overwrite:
                if not os.path.exists(tmp_file):
                    raise Exception("Tmp file was not created, need to rerun previous processes")
                else:
                    header_superpose(raw_dwi, tmp_file, outpath=subj_file)
            if shortcuts_all_folder is not None:
                subj_link = os.path.join(shortcuts_all_folder, f'{subj}_subjspace_{contrast}{ext}')
                if not os.path.exists(subj_link) or overwrite:
                    buildlink(subj_file, subj_link)

    #write the relative orientation file here
    if not os.path.isfile(orient_string) or overwrite:
        if os.path.isfile(orient_string):
            os.remove(orient_string)
        file = os.path.join(work_dir,subj+'_tmp_mask'+ext);
        cmd = 'bash ' + os.path.join(gunniespath,'find_relative_orientation_by_CoM.bash') + f' {reference_file} {file}'
        orient_relative = subprocess.getoutput(cmd)

        with open(orient_string, 'w') as f:
            f.write(orient_relative)
    else:
        orient_relative = open(orient_string, mode='r').read()

    if SAMBA_inputs_folder is not None:
        subj_orient_string = os.path.join(SAMBA_inputs_folder, f'{subj}_relative_orientation.txt')
        shutil.copy(orient_string, subj_orient_string)

    if shortcuts_all_folder is not None:
        subj_orient_string = os.path.join(shortcuts_all_folder, f'{subj}_relative_orientation.txt')
        shutil.copy(orient_string, subj_orient_string)

    #check extracted values from relative orientation vals
    orientation_out = orient_relative.split(',')[0]
    orientation_out = orientation_out.split(':')[1]
    orientation_in = orient_relative.split(',')[1]
    orientation_in = orientation_in.split(':')[1]
    if verbose:
        print(f'flexible orientation: {orientation_in}');
        print(f'reference orientation: {orientation_out}');

    #apply the orientation modification to specified contrasts
    for contrast in ['dwi', 'b0', 'mask']:
        img_in=os.path.join(work_dir,f'{subj}_tmp2_{contrast}{ext}')
        img_out=os.path.join(work_dir,f'{subj}_{contrast}{ext}')
        if not os.path.isfile(img_out) or overwrite:
            if orientation_out != orientation_in:
                print('TRYING TO REORIENT...b0 and dwi and mask')
                if os.path.exists(img_in) and (not os.path.exists(img_out) or overwrite):
                    img_transform_exec(img_in, orientation_in, orientation_out, img_out)
                    if os.path.exists(img_out):
                        os.remove(img_in)
                elif os.path.exists(img_out) and cleanup:
                    os.remove(img_in)
            else:
                shutil.move(img_in,img_out)

        if SAMBA_inputs_folder is not None:
            inputs_space_link = os.path.join(SAMBA_inputs_folder, f'{subj}_{contrast}{ext}')
            if not os.path.exists(inputs_space_link) or overwrite:
                buildlink(img_out, inputs_space_link)

        if shortcuts_all_folder is not None:
            inputs_space_link = os.path.join(shortcuts_all_folder, f'{subj}_{contrast}{ext}')
            if not os.path.exists(inputs_space_link) or overwrite:
                buildlink(img_out, inputs_space_link)


    mask = os.path.join(work_dir,f'{subj}_mask{ext}')
    b0 = os.path.join(work_dir,f'{subj}_b0{ext}')

    #if cleanup and os.path.exists(dwi_out) and os.path.exists(tmp_dwi_out):
    #    os.remove(tmp_dwi_out)

    for contrast in ['fa0', 'rd', 'ad', 'md']:
        real_file=largerfile(os.path.join(work_dir,f'*.fib.gz.{contrast}{ext}'))  # It will be fun times if we ever have more than one match to this pattern...
        #inputspace = real_file
        inputspace = os.path.join(work_dir, f'{subj}_inputspace_{contrast}{ext}')

        contrast=contrast.replace('0','')
        #linked_file=os.path.join(shortcut_dir,f'{subj}_{contrast}{ext}')
        linked_file_w=os.path.join(work_dir,f'{subj}_{contrast}{ext}')

        made_newfile = affine_superpose(dwi_out, real_file, outpath = inputspace, transpose=transpose)
        if not made_newfile:
            inputspace = real_file
        if not os.path.isfile(linked_file_w) or overwrite:
            buildlink(inputspace, linked_file_w)
        if SAMBA_inputs_folder is not None:
            #warnings.warn('should reach this!')
            blinked_file = os.path.join(SAMBA_inputs_folder, f'{subj}_{contrast}{ext}')
            if not os.path.exists(blinked_file) or overwrite:
                buildlink(inputspace, blinked_file)
                print(f'build link from {inputspace} to {blinked_file}')

        if shortcuts_all_folder is not None:
            #warnings.warn('should reach this!')
            blinked_file = os.path.join(shortcuts_all_folder, f'{subj}_{contrast}{ext}')
            if not os.path.exists(blinked_file) or overwrite:
                buildlink(inputspace, blinked_file)
                print(f'build link from {inputspace} to {blinked_file}')




    if create_subj_space_files:
        for contrast in ['fa0', 'rd', 'ad', 'md']:

            real_file = largerfile(os.path.join(work_dir,
                                                f'*.fib.gz.{contrast}{ext}'))  # It will be fun times if we ever have more than one match to this pattern...
            contrast = contrast.replace('0', '')
            subj_file_tmp = os.path.join(work_dir, f'{subj}_subjspace_tmp_{contrast}{ext}')
            subj_file = os.path.join(work_dir, f'{subj}_subjspace_{contrast}{ext}')
            if not os.path.exists(subj_file) or overwrite:
                if not os.path.exists(subj_file_tmp):
                    if orientation_out != orientation_in:
                        print('TRYING TO REORIENT...b0 and dwi and mask')
                        if os.path.exists(real_file) and (not os.path.exists(subj_file) or overwrite):
                            img_transform_exec(real_file, orientation_out, orientation_in, subj_file_tmp)
                    else:
                        shutil.copy(real_file,subj_file_tmp)
                header_superpose(raw_dwi, subj_file_tmp, outpath=subj_file)

            if shortcuts_all_folder is not None:
                subj_link = os.path.join(shortcuts_all_folder, f'{subj}_subjspace_{contrast}{ext}')
                buildlink(subj_file, subj_link)

    #if cleanup:
    #    tmp_files = glob.glob(os.path.join(work_dir, '*tmp*'))
    #    for file in tmp_files:
    #        os.remove(file)

    nii_path = os.path.join(work_dir,'nii4D_'+subj + ext)
    if os.path.exists(nii_path):
        os.remove(nii_path)

    if cleanup:
        reg_src = os.path.join(work_dir,'Reg_' + subj + f'_nii4D{ext}.src.gz')
        if os.path.exists(reg_src):
            os.remove(reg_src)
        reg_src_fib = os.path.join(work_dir,'Reg_' + subj + f'_nii4D{ext}.src.gz.dti.fib.gz')
        if os.path.exists(reg_src):
            os.remove(reg_src)

"""
    bonusmask=False

    if bonusmask:
        mask_new = mask.replace('_mask', '_new')
        bet_cmd = f'bet {b0} {mask_new} -m -f 0.3'
        os.system(bet_cmd)
        mask_new = mask_new.replace('_new','_new_mask')

        outpath=dwi_out.replace(".nii","_newmask.nii")
        applymask_samespace(dwi_out, mask_new, outpath) #(add outpath if you want to save to different location rather than overwrite)
        dwi_out=outpath

        outpath=b0_out.replace(".nii","_newmask.nii")
        applymask_samespace(b0_out, mask_new,outpath=outpath)


    for contrast in ['fa0', 'rd', 'ad', 'md']:
        real_file=largerfile(os.path.join(work_dir,f'*.fib.gz.{contrast}{ext}'))  # It will be fun times if we ever have more than one match to this pattern...
        #inputspace = real_file
        inputspace = os.path.join(work_dir, f'{subj}_inputspace_{contrast}{ext}')

        if bonusmask:
            real_file_new=real_file.replace("D.nii.gz","D_newmask.nii.gz")
            applymask_samespace(real_file, mask_new,outpath=real_file_new)
            real_file=real_file_new
        contrast=contrast.replace('0','')
        #linked_file=os.path.join(shortcut_dir,f'{subj}_{contrast}{ext}')
        linked_file_w=os.path.join(work_dir,f'{subj}_{contrast}{ext}')
        header_fix=True
        if header_fix:
            affine_superpose(dwi_out, real_file, outpath = inputspace, transpose=transpose)

        if not os.path.isfile(linked_file_w) or overwrite:
            buildlink(inputspace, linked_file_w)
        if SAMBA_inputs_folder is not None:
            blinked_file = os.path.join(SAMBA_inputs_folder, f'{subj}_{contrast}{ext}')
            if not os.path.exists(blinked_file) or overwrite:
                buildlink(inputspace, blinked_file)

"""