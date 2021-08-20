
import os
import shutil
import subprocess
from file_tools import largerfile, mkcdir, getext, buildlink
from img_transform_exec import img_transform_exec, space_transpose, affine_superpose
import glob
from basic_LPCA_denoise import basic_LPCA_denoise_func
from mask_handler import applymask_samespace, median_mask_make


def launch_preprocessing(id, raw_nii, outpath, cleanup=False, nominal_bval=4000, bonusshortcutfolder=None,
                         gunniespath="~/gunnies/", processes=1, masking="bet", ref=None, transpose=None,
                         overwrite=False, denoise='None', recenter=0, verbose=False):

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
        #print("Bvals and bvecs not found, using extractdiffdirs, however it it NOT RELIABLE, beware!")
        bvec_cmd = (f"extractdiffdirs --colvectors --writebvals --fieldsep='\t' --space=RAI {bxheader} {bvecs} {bvals}")
        os.system(bvec_cmd)

    # Make dwi for mask generation purposes.
    orient_string = os.path.join(work_dir,'relative_orientation.txt')
    tmp_mask = os.path.join(work_dir,f"{id}_tmp_dwi_mask{ext}")
    raw_dwi = os.path.join(work_dir,f"{id}_raw_dwi.nii.gz")
    orient_string = os.path.join(work_dir,"relative_orientation.txt")

    if bonusshortcutfolder is not None:
        raw_nii_link = os.path.join(bonusshortcutfolder, f"{id}_rawnii{ext}")
        if not os.path.exists(raw_nii_link) or overwrite:
            buildlink(raw_nii, raw_nii_link)
        bvecs_new = os.path.join(bonusshortcutfolder,id+"_bvecs.txt")
        bvals_new = os.path.join(bonusshortcutfolder,id+"_bvals.txt")
        if not os.path.exists(bvecs_new) or not os.path.exists(bvals_new) or overwrite:
            shutil.copyfile(bvecs,bvecs_new)
            shutil.copyfile(bvals,bvals_new)
    final_mask = os.path.join(work_dir, f'{id}_dwi_mask{ext}')
    #overwrite=True
    if (not os.path.exists(final_mask) and not os.path.exists(tmp_mask)) or overwrite:
        if not os.path.exists(raw_dwi) or overwrite:
            select_cmd = f"select_dwi_vols {raw_nii} {bvals} {raw_dwi} {nominal_bval} -m"
            os.system(select_cmd)
        if not os.path.exists(tmp_mask) or overwrite:
            if masking=="median":
                tmp = tmp_mask.replace("_mask", "")
                median_mask_make(raw_dwi, tmp, outpathmask=tmp_mask)
            elif masking=="bet":
                tmp=tmp_mask.replace("_mask", "")
                bet_cmd = f"bet {raw_dwi} {tmp} -m -n -R"
                os.system(bet_cmd)
            else:
                raise Exception("Unrecognized masking type")
    #if cleanup and (os.path.exists(tmp_mask) and os.path.exists(raw_dwi)):
    #    os.remove(raw_dwi)
    #overwrite=False
    # Run Local PCA Denoising algorithm on 4D nifti:
    masked_nii = os.path.join(work_dir, nii_name)
    if not "nii.gz" in masked_nii:
        masked_nii = masked_nii.replace(".nii", ".nii.gz")
    masked_nii = masked_nii.replace(ext, "_masked" + ext)

    if denoise=="None" or denoise is None:
        denoised_nii = masked_nii
        if not os.path.exists(masked_nii) or overwrite:
            fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
            os.system(fsl_cmd)
    else:
        denoised_nii = os.path.join(work_dir,f"LPCA_{id}_nii4D.nii.gz")
        if not os.path.exists(denoised_nii) or overwrite:
            if not os.path.exists(masked_nii) or overwrite:
                fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
                os.system(fsl_cmd)
            basic_LPCA_denoise_func(id,masked_nii,bvecs,work_dir, processes=processes,
                                    denoise=denoise, verbose=False) #to improve and make multiprocessing

    #if cleanup and os.path.exists(denoised_nii) and os.path.exists(masked_nii) and denoised_nii!=masked_nii:
    #    os.remove(masked_nii)

    # Run coregistration/eddy current correction:

    if denoise.lower()=='lpca':
        D_id=f'LPCA_{id}';
    elif denoise=="None" or denoise is None:
        D_id = f'{id}'
    coreg_nii_old = f'{outpath}/co_reg_{D_id}_m00-results/Reg_{D_id}_nii4D{ext}';
    coreg_nii = os.path.join(work_dir,f'Reg_{D_id}_nii4D{ext}')
    if not cleanup:
        coreg_nii=coreg_nii_old
    if not os.path.exists(coreg_nii) or overwrite:
        if not os.path.exists(coreg_nii_old) or overwrite:
            temp_cmd = os.path.join(gunniespath,'co_reg_4d_stack_tmpnew.bash')+f' {denoised_nii} {D_id} 0 {outpath} 0';
            os.system(temp_cmd)
        if cleanup:
            shutil.move(coreg_nii_old,coreg_nii)
    if bonusshortcutfolder is not None:
        coreg_link = os.path.join(bonusshortcutfolder,f'{id}_coreg{ext}')
        if not os.path.exists(coreg_link) or overwrite:
            buildlink(coreg_nii, coreg_link)

    toeddy=False
    if toeddy:
        #fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
        #os.system(fsl_cmd)
        eddy_cmd = f"eddy --imain={coreg_nii} --mask={tmp_mask} --acqp=acq_params.txt --index={os.path.join(work_dir,'index.txt')} --bvecs={bvecs} --bvals={bvals} --topup=topup_results --repol --out = {os.path.join(work_dir,f'Reg_{D_id}_nii4D_eddy{ext}')}"
        os.system(eddy_cmd)

    coreg_inputs=os.path.join(outpath,f'co_reg_{D_id}_m00-inputs')
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
    #overwrite=True
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
            buildlink(real_file, linked_file)
        if bonusshortcutfolder is not None:
            bonus_link = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(bonus_link) or overwrite:
                buildlink(real_file, bonus_link)

    #give real header to the temp files using md as reference
    #overwrite=True

    if ref=="md" or ref is None:
        reference=md
    elif ref=="coreg":
        reference=coreg_nii
    elif os.path.exists(ref):
        reference=ref
    #overwrite=True
    for contrast in ['dwi', 'b0', 'dwi_mask']:
        tmp_file=os.path.join(work_dir,f'{id}_tmp_{contrast}{ext}')
        tmp2_file=os.path.join(work_dir,f'{id}_tmp2_{contrast}{ext}')
        final_file=os.path.join(work_dir,f'{id}_{contrast}{ext}')
        if ((not os.path.exists(tmp2_file) and not os.path.exists(final_file)) or overwrite):
            if not os.path.exists(tmp_file):
                raise Exception("Tmp file was not created, need to rerun previous processes")
            else:
                #cmd = os.path.join(gunniespath,'nifti_header_splicer.bash')+f' {md} {file} {file}'
                cmd = os.path.join(gunniespath, 'nifti_header_splicer.bash') + f' {reference} {tmp_file} {tmp2_file}'
                os.system(cmd)
                #header_superpose(reference, tmp_file, outpath=tmp2_file)

    #write the relative orientation file here
    if not os.path.isfile(orient_string) or overwrite:
        if os.path.isfile(orient_string):
            os.remove(orient_string)
        file = os.path.join(work_dir,id+'_tmp_dwi_mask'+ext);
        cmd = 'bash ' + os.path.join(gunniespath,'find_relative_orientation_by_CoM.bash') + f' {reference} {file}'
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

    #apply the orientation modification to specified contrasts
    for contrast in ['dwi', 'b0', 'dwi_mask']:
        img_in=os.path.join(work_dir,f'{id}_tmp2_{contrast}{ext}')
        img_out=os.path.join(work_dir,f'{id}_{contrast}{ext}')
        if not os.path.isfile(img_out) or overwrite:
            if orientation_out != orientation_in:
                print('TRYING TO REORIENT...b0 and dwi and mask')
                if os.path.exists(img_in) and (not os.path.exists(img_out) or overwrite):
                    img_transform_exec(img_in, orientation_in, orientation_out, img_out, recenter=recenter)
                    if os.path.exists(img_out):
                        os.remove(img_in)
                elif os.path.exists(img_out) and cleanup:
                    os.remove(img_in)
            else:
                shutil.move(img_in,img_out)

        #linked_file=os.path.join(shortcut_dir,f'{id}_{contrast}{ext}')
        #shutil.copyfile(img_out, linked_file)
        #buildlink(img_out, linked_file)
        if bonusshortcutfolder is not None:
            bonus_link = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(bonus_link) or overwrite:
                buildlink(img_out, bonus_link)


    mask = os.path.join(work_dir,f'{id}_mask{ext}')
    b0 = os.path.join(work_dir,f'{id}_b0{ext}')

    #if cleanup and os.path.exists(dwi_out) and os.path.exists(tmp_dwi_out):
    #    os.remove(tmp_dwi_out)

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
            buildlink(dwi_out, dwi_out_blink)

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
            affine_superpose(dwi_out, real_file, transpose=transpose)

        if not os.path.isfile(linked_file) or overwrite:
            buildlink(real_file, linked_file_w)
        if bonusshortcutfolder is not None:
            blinked_file = os.path.join(bonusshortcutfolder, f'{id}_{contrast}{ext}')
            if not os.path.exists(blinked_file) or overwrite:
                buildlink(real_file, blinked_file)
