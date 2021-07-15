
import os
import numpy as np
import shutil
import subprocess
from file_tools import largerfile, mkcdir, getext, getrelativepath
from img_transform_exec import img_transform_exec, space_transpose, header_superpose
import nibabel as nib
import glob
from dipy.segment.mask import median_otsu
from dipy.io.image import load_nifti, save_nifti
from basic_LPCA_denoise import basic_LPCA_denoise_func


def buildlink(linked_file, real_file):
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        relpath = getrelativepath(real_file, linked_file)
        link_cmd=f"ln -s ./{relpath} {linked_file}"
        os.system(link_cmd)

def applymask_samespace(file, mask, outpath=None):
    #note: there should definitely be a bet option which would probably end up faster, but for some insane reason there is no
    #obvious documentation on it, so this will do for now -_-
    img_nii= nib.load(file)
    mask_nii = nib.load(mask)
    img_data = img_nii.get_fdata()
    img_data_new = img_data
    mask_data = mask_nii.get_fdata()
    img_shape = img_nii.shape
    dims = np.size(img_shape)
    if (dims == 3 or dims == 4) and mask_nii.shape[0:3] == img_nii.shape[0:3]:
        for i in np.arange(img_shape[0]):
            for j in np.arange(img_shape[1]):
                for k in np.arange(img_shape[2]):
                    if mask_data[i, j, k] == 1:
                        if dims == 3:
                            img_data_new[i, j, k] = img_data[i, j, k]
                        else:
                            for l in range(img_shape[3]):
                                img_data_new[i, j, k, l] = img_data[i, j, k, l]
                    else:
                        if dims == 3:
                            img_data_new[i, j, k] = 0
                        else:
                            for l in range(img_shape[3]):
                                img_data_new[i, j, k, l] = 0
    elif dims not in [3,4]:
        raise TypeError("Could not interpret the dimensions of the entering image")
    elif mask_nii.shape[0:3] != img_nii.shape[0:3]:
        raise TypeError("The shape of the mask and the image are not equal, therefore readjustments are needed")
    img_nii_new = nib.Nifti1Image(img_data_new, img_nii.affine, img_nii.header)
    if outpath is None:
        outpath = file
    nib.save(img_nii_new, outpath)


def median_mask_make(inpath, outpath, outpathmask=None, median_radius=4, numpass=4):

    if outpathmask is None:
        outpathmask=outpath.replace(".nii","_mask.nii")
    data, affine = load_nifti(inpath)
    data = np.squeeze(data)
    data_masked, mask = median_otsu(data, median_radius=median_radius, numpass=numpass)
    save_nifti(outpath, data_masked.astype(np.float32), affine)
    save_nifti(outpathmask, mask.astype(np.float32), affine)

def launch_preprocessing(id, raw_nii, outpath, cleanup=False, nominal_bval=4000, bonusshortcutfolder=None, gunniespath="~/gunnies/", processes=1, atlas=None, transpose=None, overwrite=False, verbose=False):

    #img_transform_exec('/Volumes/Data/Badea/Lab/human/Sinha_epilepsy/diffusion_prep_locale/diffusion_prep_00393/00393_tmp_b0.nii.gz','LPS', 'RAS',
    #'/Volumes/Data/Badea/Lab/human/Sinha_epilepsy/diffusion_prep_locale/diffusion_prep_00393/00393_b0.nii.gz')
    proc_name ="diffusion_prep_" # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline

    ## 15 June 2020, BJA: I still need to figure out the best way to pull out the non-zero bval(s) from the bval file.
    ## For now, hardcoding it to 1000 to run Whitson data. # 8 September 2020, BJA: changing to 800 for Sinha data.
    work_dir=os.path.join(outpath,proc_name+id)

    if verbose:
        print(f"Processing diffusion data with runno/id: {id}")
        print(f"Work directory is {work_dir}")
    mkcdir(work_dir)

    sbatch_folder =os.path.join(work_dir,"sbatch")
    mkcdir(sbatch_folder)

    #nii_name =${raw_nii  ##*/};
    nii_name=os.path.split(raw_nii)[1]
    niifolder = os.path.dirname(raw_nii)
    ext = ".nii.gz"
    nii_ext=getext(nii_name)
    bxheader = nii_name.replace(nii_ext,".bxh")
    bxheader = os.path.join(niifolder, bxheader)
    bvecs = os.path.join(work_dir, id+"_bvecs.txt")
    bvals =bvecs.replace("bvecs","bvals");

    if verbose:
        print(f"Original nifti is at {nii_name}\nbvecs are at {bvecs}\nbvals are at {bvals}\n")
    if not os.path.exists(bvecs):
        if verbose:
            print("Extracting diff directions")
        bvec_cmd = (f"extractdiffdirs --colvectors --writebvals --fieldsep=\t --space=RAI {bxheader} {bvecs} {bvals}")
        os.system(bvec_cmd)

    # Make dwi for mask generation purposes.
    orient_string = os.path.join(work_dir,'relative_orientation.txt')
    tmp_mask = os.path.join(work_dir,f"{id}_tmp_mask{ext}")
    raw_dwi = os.path.join(work_dir,f"{id}_raw_dwi.nii.gz")
    orient_string = os.path.join(work_dir,"relative_orientation.txt")

    if bonusshortcutfolder is not None:
        raw_nii_link = os.path.join(bonusshortcutfolder, f"{id}_rawnii{ext}")
        if not os.path.exists(raw_nii_link) or overwrite:
            buildlink(raw_nii_link, raw_nii)
        bvecs_new = os.path.join(bonusshortcutfolder,id+"_bvecs.txt")
        bvals_new = os.path.join(bonusshortcutfolder,id+"_bvals.txt")
        if not os.path.exists(bvecs_new) or not os.path.exists(bvals_new) or overwrite:
            shutil.copyfile(bvecs,bvecs_new)
            shutil.copyfile(bvals,bvals_new)
    final_mask = os.path.join(work_dir, f'{id}_mask{ext}')

    if (not os.path.exists(final_mask) and not os.path.exists(tmp_mask)) or overwrite:
        if not os.path.exists(raw_dwi) or overwrite:
            select_cmd = f"select_dwi_vols {raw_nii} {bvals} {raw_dwi} {nominal_bval} -m"
            os.system(select_cmd)

        if not os.path.exists(tmp_mask) or overwrite:
            tmp = tmp_mask.replace("_mask", "")
            median_mask_make(raw_dwi, tmp, outpathmask=tmp_mask)
            #tmp=tmp_mask.replace("_mask", "")
            #bet_cmd = f"bet {raw_dwi} {tmp} -m -n"
            #os.system(bet_cmd)
    #if cleanup and (os.path.exists(tmp_mask) and os.path.exists(raw_dwi)):
    #    os.remove(raw_dwi)

    # Run Local PCA Denoising algorithm on 4D nifti:
    denoised_nii = os.path.join(work_dir,f"LPCA_{id}_nii4D.nii.gz")
    masked_nii = os.path.join(work_dir,nii_name)
    if not "nii.gz" in masked_nii:
        masked_nii = masked_nii.replace(".nii",".nii.gz")
    masked_nii = masked_nii.replace(ext,"_masked"+ext)

    if not os.path.exists(denoised_nii) or overwrite:
        if not os.path.exists(masked_nii) or overwrite:
            fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
            os.system(fsl_cmd)
        basic_LPCA_denoise_func(id,masked_nii,bvecs,work_dir, processes=processes, verbose=False) #to improve and make multiprocessing

    #if cleanup and os.path.exists(denoised_nii) and os.path.exists(masked_nii):
    #    os.remove(masked_nii)

    # Run coregistration/eddy current correction:

    L_id=f'LPCA_{id}';
    coreg_nii_old = f'{outpath}/co_reg_{L_id}_m00-results/Reg_{L_id}_nii4D{ext}';
    coreg_nii = os.path.join(work_dir,f'Reg_{L_id}_nii4D{ext}')
    if not cleanup:
        coreg_nii=coreg_nii_old
    if not os.path.exists(coreg_nii) or overwrite:
        if not os.path.exists(coreg_nii_old) or overwrite:
            temp_cmd = os.path.join(gunniespath,'co_reg_4d_stack_tmpnew.bash')+f' {denoised_nii} {L_id} 0 {outpath}';
            os.system(temp_cmd)
        if cleanup:
            shutil.move(coreg_nii_old,coreg_nii)
    if bonusshortcutfolder is not None:
        coreg_link = os.path.join(bonusshortcutfolder,f'{id}_coreg{ext}')
        if not os.path.exists(coreg_link) or overwrite:
            buildlink(coreg_link, coreg_nii)

    coreg_inputs=os.path.join(outpath,f'co_reg_{L_id}_m00-inputs')
    coreg_work=coreg_inputs.replace('-inputs','-work')
    coreg_results=coreg_inputs.replace('-inputs','-results')
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_inputs):
        shutil.rmtree(coreg_inputs)
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_work):
        shutil.rmtree(coreg_work)
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_results):
        shutil.rmtree(coreg_results)

    # Generate tmp DWI:

    tmp_dwi_out=os.path.join(work_dir, f'{id}_tmp_dwi{ext}')
    dwi_out=os.path.join(work_dir,f'{id}_dwi{ext}')
    if not os.path.exists(dwi_out) or overwrite:
        if not os.path.exists(tmp_dwi_out) or overwrite:
            cmd=f'select_dwi_vols {coreg_nii} {bvals} {tmp_dwi_out} {nominal_bval}  -m'
            os.system(cmd)

    # Generate tmp B0:
    tmp_b0_out=os.path.join(work_dir,f'{id}_tmp_b0{ext}')
    #b0_out_link=os.path.join(shortcut_dir,f'{id}_b0{ext}')
    #if not os.path.exists(b0_out_link):
    b0_out = os.path.join(work_dir, f'{id}_b0{ext}')
    if (not os.path.exists(b0_out) and not os.path.exists(tmp_b0_out)) or overwrite:
        cmd=f'select_dwi_vols {coreg_nii} {bvals} {tmp_b0_out} 0  -m;'
        os.system(cmd)
    #elif cleanup and os.path.exists(tmp_b0_out):
    #    os.remove(tmp_b0_out)

    # Generate DTI contrasts and perform some tracking QA:
    if cleanup:
        c_string=' --cleanup '
    else:
        c_string=''

    if len(glob.glob(os.path.join(work_dir,f'*.fib.gz.md{ext}'))) == 0:
        cmd = os.path.join(gunniespath,'dti_qa_with_dsi_studio.bash')+f' {coreg_nii} {bvecs} {tmp_mask} {work_dir} {c_string}';
        os.system(cmd)

    #md is used as reference for header splicer
    md=f'{work_dir}/{id}_md{ext}';

    #create link for md
    for contrast in ['md']:
        real_file=largerfile(os.path.join(work_dir,f'*.fib.gz.{contrast}{ext}'))  #Catch the 'real file' for each contrast
        #contrast=contrast.replace('0','')   #some little nonsense mostly to deal with the 'fa0' name to become 'fa', probably sohould change it
        linked_file=os.path.join(work_dir,f'{id}_{contrast}{ext}')
        if not os.path.exists(linked_file) or overwrite:
            buildlink(linked_file, real_file)
        if bonusshortcutfolder is not None:
            bonus_link = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(bonus_link) or overwrite:
                buildlink(bonus_link, real_file)

    #give real header to the temp files using md as reference
    for contrast in ['tmp_dwi', 'tmp_b0', 'tmp_mask']:
        file = os.path.join(work_dir, f'{id}_{contrast}{ext}')
        if os.path.exists(file) or overwrite:
            #cmd = os.path.join(gunniespath,'nifti_header_splicer.bash')+f' {md} {file} {file}'
            cmd = os.path.join(gunniespath, 'nifti_header_splicer.bash') + f' {md} {file} {file}'
            os.system(cmd)

    #write the relative orientation file here

    if not os.path.isfile(orient_string) or overwrite:
        if os.path.isfile(orient_string):
            os.remove(orient_string)
        file = os.path.join(work_dir,id+'_tmp_mask'+ext);
        if atlas is None:
            atlas = md
        cmd = 'bash ' + os.path.join(gunniespath,'find_relative_orientation_by_CoM.bash') + f' {atlas} {file}'
        #curious, need to look at that find relative orientation code a bit later
        #orient_test = subprocess.run([cmd], stdout=subprocess.PIPE).stdout.decode('utf-8'
        orient_relative = subprocess.getoutput(cmd)
        with open(orient_string, 'w') as f:
            f.write(orient_relative)
    else:
        orient_relative = open(orient_string, mode='r').read()

    #check extracted values from relative orientation vals
    orientation_out = orient_relative.split(',')[0]
    orientation_out = orientation_out.split(':')[1]
    orientation_in = orient_relative.split(',')[1]
    orientation_in = orientation_in.split(':')[1]
    if verbose:
        print(f'flexible orientation: {orientation_in}');
        print(f'reference orientation: {orientation_out}');

    overwrite=True
    #apply the orientation modification to specified contrasts
    for contrast in ['dwi', 'b0', 'mask', 'md']:
        img_in=os.path.join(work_dir,f'{id}_tmp_{contrast}{ext}')
        img_out=os.path.join(work_dir,f'{id}_{contrast}{ext}')
        if not os.path.isfile(img_out) or overwrite:
            if orientation_out != orientation_in:
                print('TRYING TO REORIENT...b0 and dwi and mask')
                if os.path.exists(img_in) and not os.path.exists(img_out):
                    img_transform_exec(img_in, orientation_in, orientation_out, img_out)
                    if os.path.exists(img_out):
                        os.remove(img_in)
                elif os.path.exists(img_out):
                    os.remove(img_in)
            else:
                shutil.move(img_in,img_out)

        #linked_file=os.path.join(shortcut_dir,f'{id}_{contrast}{ext}')
        #shutil.copyfile(img_out, linked_file)
        #buildlink(linked_file, img_out)
        if bonusshortcutfolder is not None:
            bonus_link = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(bonus_link) or overwrite:
                buildlink(bonus_link, img_out)

    mask = os.path.join(work_dir,f'{id}_mask{ext}')
    b0 = os.path.join(work_dir,f'{id}_b0{ext}')

    if cleanup and os.path.exists(dwi_out) and os.path.exists(tmp_dwi_out):
        os.remove(tmp_dwi_out)

    #refine mask (to check!!!!)

    bonusmask=False
    if bonusmask:
        mask_new = mask.replace('_mask', '_new')
        bet_cmd = f'bet {b0} {mask_new} -m -f 0.3'
        os.system(bet_cmd)
        mask_new = mask_new.replace('_new','_new_mask')

    if bonusshortcutfolder is not None:
        dwi_out_blink=os.path.join(bonusshortcutfolder,f'{id}_dwi{ext}')
        if not os.path.exists(dwi_out_blink) or overwrite:
            buildlink(dwi_out_blink,dwi_out)

    if bonusmask:
        outpath=dwi_out.replace(".nii","_newmask.nii")
        applymask_samespace(dwi_out, mask_new, outpath) #(add outpath if you want to save to different location rather than overwrite)
        dwi_out=outpath

        outpath=b0_out.replace(".nii","_newmask.nii")
        applymask_samespace(b0_out, mask_new,outpath=outpath)

    for contrast in ['fa0', 'rd', 'ad', 'md']:
        real_file=largerfile(os.path.join(work_dir,f'*.fib.gz.{contrast}{ext}'))  # It will be fun times if we ever have more than one match to this pattern...
        if bonusmask:
            real_file_new=real_file.replace("D.nii.gz","D_newmask.nii.gz")
            applymask_samespace(real_file, mask_new,outpath=real_file_new)
            real_file=real_file_new
        contrast=contrast.replace('0','')
        #linked_file=os.path.join(shortcut_dir,f'{id}_{contrast}{ext}')
        linked_file_w=os.path.join(work_dir,f'{id}_{contrast}{ext}')
        header_fix=True
        if header_fix:
            outpath = linked_file_w.replace('.nii','_test.nii')
            header_superpose(dwi_out, real_file, transpose=None,outpath=outpath)

        #shutil.copyfile(real_file,linked_file)
        buildlink(linked_file_w,real_file)
        if bonusshortcutfolder is not None:
            blinked_file = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(blinked_file) or overwrite:
                buildlink(blinked_file, real_file)
