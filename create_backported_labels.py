
from transform_handler import img_transform_exec, space_transpose, header_superpose
import os, re, sys, io, struct, socket, datetime
from pathlib import Path
import os
from file_tools import mkcdir, check_files
import glob
import warnings
import shutil


def get_info_SAMBA_headfile(SAMBA_headfile, verbose=False):

    if os.path.exists(SAMBA_headfile):
        with open(SAMBA_headfile, 'rb') as source:
            if verbose: print('INFO    : Extracting acquisition parameters')
            header_size = source.read(4)
            header_size = struct.unpack('I', header_size)
            if verbose: print('INFO    : Header size = ', int(header_size[0]))
            i = 0
            stopsign = 200
            for line in source:
                pattern1 = 'original_study_orientation'
                rx1 = re.compile(pattern1, re.IGNORECASE | re.MULTILINE | re.DOTALL)
                pattern2 = 'working_image_orientation'
                rx2 = re.compile(pattern2, re.IGNORECASE | re.MULTILINE | re.DOTALL)
                discount = re.compile('#', re.IGNORECASE | re.MULTILINE | re.DOTALL)
                i += 1
                if i == stopsign:
                    print("hi")
                for a in rx1.findall(str(line)):
                    if len(discount.findall(str(line))) == 0:
                        orig_orientation = str(line).split('=')[1]
                        orig_orientation = orig_orientation.split('\\')[0]
                for a in rx2.findall(str(line)):
                    if len(discount.findall(str(line))) == 0:
                        working_orientation = str(line).split('=')[1]
                        working_orientation = working_orientation.split('\\')[0]
    else:
        raise Exception('Cannot find SAMBA_headfile')

    return orig_orientation, working_orientation


def create_backport_labels(subject, mainpath, project_name, atlas_labels, orient_string, preppath = None, headfile=None, overwrite=False, verbose=True):

    #mainpath = "/Volumes/Data/Badea/Lab/mouse/"
    #gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
    #project_name = "VBM_21ADDecode03_IITmean_RPI_fullrun"
    #atlas_labels = "/Volumes/Data/Badea/Lab/atlas/IITmean_RPI/IITmean_RPI_labels.nii.gz"

    out_dir = os.path.join(mainpath, f"{project_name}-results","connectomics",subject)
    work_dir = os.path.join(mainpath, f"{project_name}-work")
    dirty_dir = os.path.join(mainpath,"burn_after_reading")
    mkcdir(dirty_dir)

    template_type_prefix = os.path.basename(os.path.dirname(glob.glob(os.path.join(work_dir,"dwi","SyN*/"))[0]))
    template_runs = glob.glob((os.path.join(work_dir,"dwi",template_type_prefix,"*/")))
    mymax=-1
    for template_run in template_runs:
        if "NoNameYet" in template_run and template_run[-4:-2]=="_i":
            if int(template_run[-2])>mymax:
                mymax=int(template_run[-2])
                final_template_run=template_run
    if mymax==-1:
        raise Exception(f"Could not find template runs in {os.path.join(mainpath, f'{project_name}-work','dwi',template_type_prefix)}")
    #template_type_prefix = "faMDT_NoNameYet_n37_i6"
    #final_template_run = "SyN_0p5_3_0p5_fa"

    #final_ref = f"/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale/diffusion_prep_{subject}/Reg_LPCA_{subject}_nii4D.nii.gz";
    if os.path.exists(orient_string):
        orient_relative = open(orient_string, mode='r').read()
        orientation_out = orient_relative.split(',')[0]
        orientation_out = orientation_out.split(':')[1]
        orientation_in = orient_relative.split(',')[1]
        orientation_in = orientation_in.split(':')[1]
    else:
        orientation_out = "*"
        warnings.warn("Could not find orientation file, may cause errors later")
        orientation_in = "LPS"
        orientation_out = "RAS"

    if headfile is not None:
        SAMBA_orientation_in, SAMBA_orientation_out = get_info_SAMBA_headfile(headfile)
    else:
        SAMBA_orientation_in, SAMBA_orientation_out = 'RAS', 'RAS'


    trans = os.path.join(work_dir,"preprocess","base_images","translation_xforms",f"{subject}_0DerivedInitialMovingTranslation.mat")
    rigid =os.path.join(work_dir,"dwi",f"{subject}_rigid.mat")
    affine = os.path.join(work_dir,"dwi",f"{subject}_affine.mat")
    MDT_to_subject = os.path.join(final_template_run,"reg_diffeo",f"MDT_to_{subject}_warp.nii.gz")
    label_name = os.path.basename(atlas_labels)
    label_name = label_name.split("_labels")[0]
    MDT_to_atlas_affine = os.path.join(final_template_run,"stats_by_region","labels","transforms",f"MDT_*_to_{label_name}_affine.mat")
    atlas_to_MDT = os.path.join(final_template_run,"stats_by_region","labels","transforms",f"{label_name}_to_MDT_warp.nii.gz")

    [trans, rigid, affine, MDT_to_subject, MDT_to_atlas_affine, atlas_to_MDT], _ = check_files([trans,rigid,affine,MDT_to_subject,MDT_to_atlas_affine,atlas_to_MDT])


    preprocess_ref = os.path.join(work_dir,"preprocess",f"{subject}_fa_masked.nii.gz")
    preprocess_labels = os.path.join(dirty_dir,f"{subject}_preprocess_labels.nii.gz")
    fixed_preprocess_labels = os.path.join(dirty_dir,f"{subject}_fixed_preprocess_labels.nii.gz")
    coreg_labels = os.path.join(dirty_dir,f"{subject}_{orientation_in}_labels.nii.gz")
    coreg_reorient_labels = os.path.join(dirty_dir,f"{subject}_{SAMBA_orientation_in}_labels.nii.gz")

    # final_ref="/mnt/munin6/Badea/Lab/mouse/co_reg_LPCA_${subject:1:5}_m00-results/Reg_LPCA_${subject:1:5}_nii4D.nii.gz";
    if preppath is None:
        inputsfolder = os.path.join(mainpath, f"{project_name}-inputs")
        final_ref = os.path.join(inputsfolder, f"{subject}_coreg.nii.gz")
        if os.path.exists(final_ref):
            if Path(final_ref).is_symlink():
                final_ref=str(Path(final_ref).resolve())
        else:
            raise Exception("Could not find final registered subject file")


    symbolic_ref = os.path.join(out_dir,f"{subject}_Reg_LPCA_nii4D.nii.gz")
    final_labels = os.path.join(out_dir,f"{subject}_IITmean_preprocess_labels.nii.gz")
    final_labels_backup = os.path.join(dirty_dir,f"{subject}_IITmean_preprocess_labels.nii.gz")
    superpose=True
    if not os.path.exists(final_labels) or overwrite:
        if verbose:
            print(f"Backporting labels to raw space for subject: {subject} to {final_labels}")
        if not os.path.exists(preprocess_labels) or overwrite:
            #Note, it used to be MultiLabel[1.0x1.0x1.0,2] but it seems like the default parameters tend to be based on the image themselves and work fine,
            #so be careful but I removed it for now so it would work on different image sets
            cmd = f"antsApplyTransforms -v 1 -d 3 -i {atlas_labels} -o {preprocess_labels} -r {preprocess_ref} -n MultiLabel -t [{trans},1] [{rigid},1] [{affine},1] {MDT_to_subject} [{MDT_to_atlas_affine},1] {atlas_to_MDT}"
            if verbose:
                print(f"Runnings the Ants apply transforms to {atlas_labels}")
                print(cmd)
            os.system(cmd)
        if os.path.exists(preprocess_labels) and ((not os.path.exists(fixed_preprocess_labels)) or overwrite):
            header_superpose(final_ref, preprocess_labels, outpath=fixed_preprocess_labels)

        if os.path.exists(fixed_preprocess_labels) and (not os.path.exists(coreg_labels) or overwrite):

            if orientation_out != orientation_in:
                if verbose:
                    print(f"Reorientation from {orientation_out} back to {orientation_in}")
                img_transform_exec(fixed_preprocess_labels, orientation_out, orientation_in, coreg_labels)
            else:
                if verbose:
                    print(f"Orientations are the same, skipping reorientation")
                shutil.copy(fixed_preprocess_labels, coreg_labels)

        if os.path.exists(coreg_labels) and (not os.path.exists(coreg_reorient_labels) or overwrite):

            if SAMBA_orientation_out != SAMBA_orientation_in:
                if verbose:
                    print(f"Reorientation from {SAMBA_orientation_out} back to {SAMBA_orientation_in}")
                img_transform_exec(coreg_labels, SAMBA_orientation_out, SAMBA_orientation_in,
                                   output_path=coreg_reorient_labels)
            else:
                if verbose:
                    print(f"Orientations are the same, skipping reorientation")
                if coreg_labels != coreg_reorient_labels:
                    warnings.warn('This should never happen, investigate if it does')
                    shutil.copy(fixed_preprocess_labels, coreg_labels)

        if os.path.exists(coreg_reorient_labels):
            header_superpose(final_ref, coreg_reorient_labels, outpath=final_labels)
            #cmd = f"{os.path.join(gunniespath, 'nifti_header_splicer.bash')} {final_ref} {coreg_labels} {final_labels}"
            #os.system(cmd)

        if os.path.exists(final_labels):
            if verbose:
                print(f"Applying fsl maths to {final_labels}")
            cmd = f"fslmaths {final_labels} -add 0 {final_labels} -odt short"
            os.system(cmd)
            shutil.copy(final_labels, final_labels_backup)
    else:
        print(f"Already calculated the label file for subject {subject}")

    skip_making_data_package_for_tractography = True

    if not skip_making_data_package_for_tractography:

        print('if you see this, make the bval fixes')
        if not os.path.exists9(symbolic_ref):
            cmd=f"ln - s {final_ref} {symbolic_ref}"

        bvals = os.path.join(mainpath, f"diffusion_prep_{subject}", f"{subject}_bvals.txt")
        bvecs = os.path.join(mainpath, f"diffusion_prep_{subject}", f"{subject}_bvecs.txt")

        bval_copy = os.path.join(out_dir, f"{subject}_bvals.txt")
        bvec_copy = os.path.join(out_dir, f"{subject}_bvecs.txt")

        if not os.path.exists(bval_copy):
            shutil.copy(bvals, bval_copy)

        if not os.path.exists(bvec_copy):
            shutil.copy(bvecs, bvec_copy)


def convert_images_templatespace(subject, mainpath, project_name, subject_img_path, orient_string, preppath = None, gunniespath="/Users/alex/bass/gitfolder/wuconnectomes/gunnies/", recenter=0, verbose=True):

    #mainpath = "/Volumes/Data/Badea/Lab/mouse/"
    #gunniespath = "/Users/alex/bass/gitfolder/wuconnectomes/gunnies/"
    #project_name = "VBM_21ADDecode03_IITmean_RPI_fullrun"
    #atlas_labels = "/Volumes/Data/Badea/Lab/atlas/IITmean_RPI/IITmean_RPI_labels.nii.gz"

    out_dir = os.path.join(mainpath, f"{project_name}-results","connectomics",subject)
    work_dir = os.path.join(mainpath, f"{project_name}-work")
    dirty_dir = os.path.join(mainpath,"burn_after_reading")
    mkcdir(dirty_dir)

    template_type_prefix = os.path.basename(os.path.dirname(glob.glob(os.path.join(work_dir,"dwi","SyN*/"))[0]))
    template_runs = glob.glob((os.path.join(work_dir,"dwi",template_type_prefix,"*/")))
    mymax=-1
    for template_run in template_runs:
        if "NoNameYet" in template_run and template_run[-4:-2]=="_i":
            if int(template_run[-2])>mymax:
                mymax=int(template_run[-2])
                final_template_run=template_run
    if mymax==-1:
        raise Exception(f"Could not find template runs in {os.path.join(mainpath, f'{project_name}-work','dwi',template_type_prefix)}")
    #template_type_prefix = "faMDT_NoNameYet_n37_i6"
    #final_template_run = "SyN_0p5_3_0p5_fa"

    #final_ref = f"/mnt/munin6/Badea/Lab/human/AD_Decode/diffusion_prep_locale/diffusion_prep_{subject}/Reg_LPCA_{subject}_nii4D.nii.gz";
    if os.path.exists(orient_string):
        orient_relative = open(orient_string, mode='r').read()
        orientation_out = orient_relative.split(',')[0]
        orientation_out = orientation_out.split(':')[1]
        orientation_in = orient_relative.split(',')[1]
        orientation_in = orientation_in.split(':')[1]
    else:
        orientation_out = "*"
        warnings.warn("Could not find orientation file, may cause errors later")
        orientation_in = "LPS"
        orientation_out = "RAS"

    trans = os.path.join(work_dir,"preprocess","base_images","translation_xforms",f"{subject}_0DerivedInitialMovingTranslation.mat")
    rigid =os.path.join(work_dir,"dwi",f"{subject}_rigid.mat")
    affine = os.path.join(work_dir,"dwi",f"{subject}_affine.mat")
    #MDT_to_subject = os.path.join(final_template_run,"reg_diffeo",f"MDT_to_{subject}_warp.nii.gz")
    #MDT_to_atlas_affine = os.path.join(final_template_run,"stats_by_region","labels","transforms",f"MDT_*_to_{label_name}_affine.mat")
    #atlas_to_MDT = os.path.join(final_template_run,"stats_by_region","labels","transforms",f"{label_name}_to_MDT_warp.nii.gz")

    [trans, rigid, affine], _ = check_files([trans,rigid,affine])


    preprocess_ref = os.path.join(work_dir,"preprocess",f"{subject}_fa_masked.nii.gz")
    preprocess_labels = os.path.join(dirty_dir,f"{subject}_preprocess_labels.nii.gz")
    fixed_preprocess_labels = os.path.join(dirty_dir,f"{subject}_fixed_preprocess_labels.nii.gz")
    coreg_labels = os.path.join(dirty_dir,f"{subject}_{orientation_in}_labels.nii.gz")

    # final_ref="/mnt/munin6/Badea/Lab/mouse/co_reg_LPCA_${subject:1:5}_m00-results/Reg_LPCA_${subject:1:5}_nii4D.nii.gz";
    if preppath is None:
        inputsfolder = os.path.join(mainpath, f"{project_name}-inputs")
        final_ref = os.path.join(inputsfolder, f"{subject}_coreg.nii.gz")
        if os.path.exists(final_ref):
            if Path(final_ref).is_symlink():
                final_ref=str(Path(final_ref).resolve())
        else:
            raise Exception("Could not find final registered subject file")

    mytype="fa"
    symbolic_ref = os.path.join(out_dir,f"{subject}_Reg_LPCA_nii4D.nii.gz")
    final_image_test = os.path.join(out_dir,f"{subject}_{mytype}_to_MDT_test.nii.gz")

    overwrite = False
    if not os.path.exists(final_image_test) or overwrite:
        if verbose:
            print(f"Porting subject to MDT for subject: {subject} to {final_image_test}")

        img_transposed_test = os.path.join(dirty_dir,f"{subject}_{type}_transposed_test.nii.gz")
        img_transposed_test_2 = os.path.join(dirty_dir,f"{subject}_{type}_transposed_test_2.nii.gz")

        if os.path.exists(subject_img_path) and preprocess_ref:
            header_superpose(preprocess_ref, subject_img_path, outpath=img_transposed_test)
            cmd = f"antsApplyTransforms -v 1 -d 3 -i {subject_img_path} -o {img_transposed_test_2} -r {preprocess_ref} -n MultiLabel -t [{trans}] [{rigid},1] [{affine},1]"
            if verbose:
                print(f"Runnings the Ants apply transforms to {subject_img_path}")
                print(cmd)
            os.system(cmd)

        if not os.path.exists(preprocess_labels) or overwrite:
            #Note, it used to be MultiLabel[1.0x1.0x1.0,2] but it seems like the default parameters tend to be based on the image themselves and work fine,
            #so be careful but I removed it for now so it would work on different image sets
            cmd = f"antsApplyTransforms -v 1 -d 3 -i {atlas_labels} -o {preprocess_labels} -r {preprocess_ref} -n MultiLabel -t [{trans},1] [{rigid},1] [{affine},1] {MDT_to_subject} [{MDT_to_atlas_affine},1] {atlas_to_MDT}"
            if verbose:
                print(f"Runnings the Ants apply transforms to {atlas_labels}")
                print(cmd)
            os.system(cmd)

        if os.path.exists(preprocess_labels) and not ((os.path.exists(fixed_preprocess_labels)) or overwrite):
            header_superpose(final_ref, preprocess_labels, outpath=fixed_preprocess_labels)

        if os.path.exists(fixed_preprocess_labels) and not (os.path.exists(coreg_labels) or overwrite):

            if orientation_out != orientation_in:
                if verbose:
                    print(f"Reorientation from {orientation_out} back to {orientation_in}")
                img_transform_exec(fixed_preprocess_labels, orientation_out, orientation_in, coreg_labels, recenter=recenter)
            else:
                if verbose:
                    print(f"Orientations are the same, skipping reorientation")
                shutil.copy(fixed_preprocess_labels, coreg_labels)

        if os.path.exists(coreg_labels):
            header_superpose(final_ref, coreg_labels, outpath=final_labels)
            #cmd = f"{os.path.join(gunniespath, 'nifti_header_splicer.bash')} {final_ref} {coreg_labels} {final_labels}"
            #os.system(cmd)

        if os.path.exists(final_labels):
            if verbose:
                print(f"Applying fsl maths to {final_labels}")
            cmd = f"fslmaths {final_labels} -add 0 {final_labels} -odt short"
            os.system(cmd)
    else:
        print(f"Already calculated the label file for subject {subject}")

    skip_making_data_package_for_tractography = True

    if not skip_making_data_package_for_tractography:

        print('if you see this, make the bval fixes')
        if not os.path.exists9(symbolic_ref):
            cmd=f"ln - s {final_ref} {symbolic_ref}"

        bvals = os.path.join(mainpath, f"diffusion_prep_{subject}", f"{subject}_bvals.txt")
        bvecs = os.path.join(mainpath, f"diffusion_prep_{subject}", f"{subject}_bvecs.txt")

        bval_copy = os.path.join(out_dir, f"{subject}_bvals.txt")
        bvec_copy = os.path.join(out_dir, f"{subject}_bvecs.txt")

        if not os.path.exists(bval_copy):
            shutil.copy(bvals, bval_copy)

        if not os.path.exists(bvec_copy):
            shutil.copy(bvecs, bvec_copy)