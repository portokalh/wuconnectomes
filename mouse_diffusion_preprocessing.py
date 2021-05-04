
import os
import numpy as np
import shutil
import subprocess
from file_tools import largerfile, mkcdir, getext
from dipy.io.image import load_nifti, save_nifti


def launch_preprocessing(id, raw_nii, outpath, cleanup=False):


    gunniespath = "/Users/alex/bass/gitfolder/gunnies/"
    proc_name ="diffusion_prep_" # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline

    ## 15 June 2020, BJA: I still need to figure out the best way to pull out the non-zero bval(s) from the bval file.
    ## For now, hardcoding it to 1000 to run Whitson data. # 8 September 2020, BJA: changing to 800 for Sinha data.

    nominal_bval=4000;

    print(f"Processing diffusion data with runno/id: {id}")


    work_dir=os.path.join(outpath,proc_name+id)

    print(work_dir)
    mkcdir(work_dir)

    sbatch_folder =os.path.join(work_dir,"sbatch")
    mkcdir(sbatch_folder)

    #nii_name =${raw_nii  ##*/};
    nii_name=os.path.split(raw_nii)[1]
    ext = ".nii.gz"
    nii_ext=getext(nii_name)
    bxheader = nii_name.replace(nii_ext,".bxh")

    bvecs = os.path.join(work_dir, id+"_bvecs.txt")
    bvals =bvecs.replace("bvecs","bvals");

    print(nii_name)
    print(bvecs)
    print(bvals)

    if not os.path.exists(bvecs):
        bvec_cmd = (f"extractdiffdirs --colvectors --writebvals --fieldsep=\t --space=RAI {bxheader} {bvecs} {bvals}")
        os.system(bvec_cmd)

    # Make dwi for mask generation purposes.
    tmp_mask = os.path.join(work_dir,f"{id}_tmp_mask{ext}")
    mask=os.path.join(work_dir,f"{id}_mask{ext}");
    raw_dwi = os.path.join(work_dir,f"{id}_raw_dwi.nii.gz")
    if not os.path.exists(mask):
        if not os.path.exists(raw_dwi):
            select_cmd = f"select_dwi_vols {raw_nii} {bvals} {raw_dwi} {nominal_bval} -m"
            os.system(select_cmd)

        if not os.path.exists(tmp_mask):
            tmp=tmp_mask.replace("_mask", "")
            bet_cmd = f"bet {raw_dwi} {tmp} -m -n"
            os.system(bet_cmd)

    if cleanup and (os.path.exists(tmp_mask) and os.path.exists(raw_dwi)):
        os.remove(raw_dwi)

    # Run Local PCA Denoising algorithm on 4D nifti:
    # Note: the LPCA python has a logical switch that will produce another, confounding mask (though arguably better). This is currently switched off, but can cause confusion and delay if switched on and not properly documented.
    denoised_nii = os.path.join(work_dir,f"LPCA_{id}_nii4D.nii.gz")
    masked_nii = os.path.join(work_dir,nii_name)
    masked_nii = masked_nii.replace(".nii",".nii.gz")
    masked_nii = masked_nii.replace(ext,"_masked"+ext)

    tempcheck=False
    if tempcheck:
        from basic_LPCA_denoise_new import basic_LPCA_denoise_func
        basic_LPCA_denoise_func(id,masked_nii,bvecs,work_dir)
    if not os.path.exists(denoised_nii):
        if not os.path.exists(masked_nii):
            fsl_cmd = f"fslmaths {raw_nii} -mas {tmp_mask} {masked_nii} -odt 'input'";
            os.system(fsl_cmd)
        from basic_LPCA_denoise_new import basic_LPCA_denoise_func
        basic_LPCA_denoise_func(id,masked_nii,bvecs,work_dir)

    if cleanup and os.path.exists(denoised_nii) and os.path.exists(masked_nii):
        os.remove(masked_nii)

    # Run coregistration/eddy current correction:

    L_id=f"LPCA_{id}";
    coreg_nii=f"{outpath}/co_reg_{L_id}_m00-results/Reg_{L_id}_nii4D{ext}";
    if not os.path.exists(coreg_nii):
        temp_cmd = os.path.join(gunniespath,"co_reg_4d_stack_tmpnew.bash")+f" {denoised_nii} {L_id} 0 {outpath}";
        os.system(temp_cmd)

    coreg_inputs=os.path.join(outpath,f"co_reg_{L_id}_m00-inputs")
    coreg_work=coreg_inputs.replace("-inputs","-work")
    if cleanup and os.path.exists(coreg_nii) and os.path.isdir(coreg_inputs):
        os.remove(coreg_inputs)

    if cleanup and os.path.exists(coreg_nii) and os.path.isdir (coreg_work):
        os.remove(coreg_work)

    # Generate tmp DWI:
    tmp_dwi_out=os.path.join(work_dir, f"{id}_tmp_dwi{ext}")
    dwi_out=os.path.join(work_dir,f"{id}_dwi{ext}")
    if not os.path.exists(dwi_out):
        if not os.path.exists(tmp_dwi_out):
            cmd=f"select_dwi_vols {coreg_nii} {bvals} {tmp_dwi_out} {nominal_bval}  -m"
            os.system(cmd)
    #elif cleanup and os.path.exists(tmp_dwi_out):
    #    os.remove(tmp_dwi_out)

    # Generate tmp B0:
    tmp_b0_out=os.path.join(work_dir,f"{id}_tmp_b0{ext}")
    b0_out=os.path.join(work_dir,f"{id}_b0{ext}")
    if not os.path.exists(b0_out):
        if not os.path.exists(tmp_b0_out):
            cmd=f"select_dwi_vols {coreg_nii} {bvals} {tmp_b0_out} 0  -m;"
            os.system(cmd)
    elif cleanup and os.path.exists(tmp_b0_out):
        os.remove(tmp_b0_out)

    # Generate DTI contrasts and perform some tracking QA:
    c_string='';
    if cleanup:
        c_string=' --cleanup '

    cmd = os.path.join(gunniespath,"dti_qa_with_dsi_studio.bash")+f" {coreg_nii} {bvecs} {tmp_mask} {work_dir} {c_string}";
    os.system(cmd)
    # Create RELATIVE links to DSI_studio outputs:
    for contrast in ["fa0", "rd", "ad", "md"]:
        real_file=largerfile(os.path.join(work_dir,f"*.fib.gz.{contrast}{ext}"))  # It will be fun times if we ever have more than one match to this pattern...
        contrast=contrast.replace("0","")
        linked_file=os.path.join(work_dir,f"{id}_{contrast}{ext}")
        # We need to make sure that we didn't accidentally link to a non-existent file.
        if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
            os.unlink(linked_file)

        if not os.path.islink(linked_file) and os.path.isfile(real_file):
            real_file_name=os.path.basename(real_file)
            link_cmd=f"ln -s ./{real_file_name} {linked_file}"
            os.system(link_cmd)

    if os.path.exists(dwi_out):
        os.remove(dwi_out)
    linked_file=dwi_out
    real_file=tmp_dwi_out
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        real_file_name=os.path.basename(real_file)
        link_cmd=f"ln -s ./{real_file_name} {linked_file}"
        os.system(link_cmd)

    linked_file=b0_out
    real_file=tmp_b0_out
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        real_file_name=os.path.basename(real_file)
        link_cmd=f"ln -s ./{real_file_name} {linked_file}"
        os.system(link_cmd)

    linked_file=mask
    real_file=tmp_mask
    if os.path.islink(linked_file) and not os.path.exists(os.readlink(linked_file)):
        os.unlink(linked_file)
    if not os.path.islink(linked_file) and os.path.isfile(real_file):
        real_file_name=os.path.basename(real_file)
        link_cmd=f"ln -s ./{real_file_name} {linked_file}"
        os.system(link_cmd)
    # This should prevent linking to a non-existent file in the future.
        #if os.path.islink(linked_file) and os.path.exists(real_file):
        #    ln -s "./{real_file##*/}" {linked_file};

    # For Whitson data (at least) we empirically know that we need to move the data from LAS to RAS:
    # Enforce header consistency, based on fa from DSI Studio before reorienting data (will break if there's an affine xform in header)

    # Hopefully the new auto-detect-orientation code will work...

    img_xform_exec='/Users/alex/MATLAB/MATLAB_BJ/matlab_execs_for_SAMBA//img_transform_executable/run_img_transform_exec.sh';
    mat_library='/Users/alex/MATLAB/MATLAB_BJ/MATLAB2015b_runtime/v90';

    # Enforce header consistency, based on mask--no, scratch that--based on fa from DSI Studio
    # No scratch that--fa doesn't necessarily have a consistent center of mass!
    md=f"{work_dir}/{id}_md{ext}";

    # I'm not going to say specifically why, but, we only need to compare orientations between the mask and the md.
    # This only works on clean runs, and won't catch b0s or dwis that have been previously incorrectly turned.

    if not os.path.exists(tmp_mask):
        for contrast in ["tmp_dwi", "tmp_b0", "tmp_mask"]:
            file=os.path.join(work_dir,f"{id}_{contrast}{ext}")
            if os.path.exists(file):
                cmd=f"/Users/alex/bass/gitfolder/gunnies/nifti_header_splicer.bash {md} {file} {file}"
                os.system(cmd)

        file=os.path.join(work_dir,f"{id}_tmp_mask{ext}")
        orient_test=os.path.join(gunniespath, "find_relative_orientation_by_CoM.bash")+f" {md} {file}";
        orientation_output = subprocess.check_output(orient_test, shell=True)
        orientation_out = orientation_output.split(",")[0]
        orientation_out = orientation_out.split(":")[1]
        orientation_in = orientation_output.split(",")[1]
        orientation_in = orientation_in.split(":")[1]
        #orientation_in=$(echo ${orient_test} | cut -d ',' -f2 | cut -d ':' -f2);
        #orientation_out=$(echo ${orient_test} | cut -d ',' -f1 | cut -d ':' -f2);
        print(f"flexible orientation: {orientation_in}");
        print(f"reference orientation: {orientation_out}");

        for contrast in ["dwi", "b0", "mask"]:
            img_in=os.path.join(work_dir,f"{id}_tmp_{contrast}{ext}")
            img_out=os.path.join(work_dir,f"{id}_{contrast}{ext}")

            if orientation_out != orientation_in:
                print("TRYING TO REORIENT...b0 and dwi and mask")
                if os.path.exists(img_in) and not os.path.exists(img_out):
                    reorient_cmd=f"{img_xform_exec} {mat_library} {img_in} {orientation_in} {orientation_out} {img_out}"
                    os.system(reorient_cmd)
                    if os.path.exists(img_out):
                        os.remove(img_in)
                elif os.path.exists(img_out):
                    os.remove(img_in)
            else:
                shutil.move(img_in,img_out)

    # Turning off the following code block for now...
    """         
    if (( 0 ));then
    file=${work_dir} / ${id}_mask.${ext};
    orient_test=f"( / Users / alex / bass / gitfolder / gunnies / find_relative_orientation_by_CoM.bash {md} {file});

    orientation_in=$(echo ${orient_test} | cut -d ',' -f2 | cut -d ':' -f2);
    orientation_out=$(echo ${orient_test} | cut -d ',' -f1 | cut -d ':' -f2);
    echo "flexible orientation: ${orientation_in}";
    echo "reference orientation: ${orientation_out}";
    if[["${orientation_out}" != "${orientation_in}"]];then
    echo "TRYING TO REORIENT...MASK!";
    img_in=${work_dir} / ${id}_tmp2_mask.${ext};
    img_out=${work_dir} / ${id}_mask.${ext};

    mv $img_out $img_in;
    reorient_cmd="${img_xform_exec} ${mat_library} ${img_in} ${orientation_in} ${orientation_out} ${img_out}";
    ${reorient_cmd};
    if os.path.existsimg_out}]];then
    rm $img_in;
    fi
    fi
    fi
    """
    mddata, mdaffine = load_nifti(md)
    for contrast in ["dwi", "b0", "mask"]:
        file=os.path.join(work_dir,f"{id}_{contrast}{ext}")
        olddata, oldaffine=load_nifti(file)
        save_nifti(file, olddata, mdaffine)
        #cmd = os.path.join(gunniespath,"nifti_header_splicer.bash")+f" {md} {file} {file}";
        #os.system(cmd)


    # Apply updated mask to dwi and b0.
    mask=os.path.join(work_dir,f"{id}_mask{ext}")

    # for contrast in dwi b0;do
    #    file=${work_dir}/${id}_${contrast}.${ext};
    #    fslmaths ${file} -mas ${mask} ${file} -odt "input";
    # done

    #if cleanup and os.path.exists(tmp_mask) and os.path.exists(mask):
    #    os.remove(tmp_mask)

    #if cleanup and os.path.exists(tmp_dwi_out) and os.path.exists(dwi_out):
    #    os.remove(tmp_dwi_out)

    #if cleanup and os.path.exists(tmp_b0_out) and os.path.exists(b0_out):
    #    os.remove(tmp_b0_out);