
proc_name="diffusion_prep"; # Not gonna call it diffusion_calc so we don't assume it does the same thing as the civm pipeline

## 15 June 2020, BJA: I still need to figure out the best way to pull out the non-zero bval(s) from the bval file.
## For now, hardcoding it to 1000 to run Whitson data. # 8 September 2020, BJA: changing to 800 for Sinha data.

nominal_bval='4000';

id=$1;
raw_nii=$2;
BIGGUS_DISKUS=$3;
no_cleanup=$4;

if [[ "x1x" == "x${no_cleanup}x" ]];then
    cleanup=0;
else
    cleanup=1;
fi

echo "Processing diffusion data with runno/id: ${id}.";


work_dir=${BIGGUS_DISKUS}/${proc_name}_${id};

echo "$work_dir"
if [[ ! -d ${work_dir} ]];then
    mkdir -pm 775 ${work_dir};
fi

sbatch_folder=${work_dir}/sbatch;
if [[ ! -d ${sbatch_folder} ]];then
    mkdir -pm 775 ${sbatch_folder};
fi

nii_name=${raw_nii##*/};
ext="nii.gz";

bxheader=${raw_nii/nii.gz/bxh}
#    echo "Input nifti = ${raw_nii}.";
#    echo "Input header = ${bxheader}.";

bvecs=${work_dir}/${id}_bvecs.txt;
bvals=${bvecs/bvecs/bvals};

echo "Nii name"
echo $nii_name
echo $bvecs
echo $bvals


if [[ ! -f ${bvecs} ]];then
    echo "sweet home alabama"
    bvec_cmd="extractdiffdirs --colvectors --writebvals --fieldsep=\t --space=RAI ${bxheader} ${bvecs} ${bvals}";
    echo $bvec_cmd
    $bvec_cmd;
fi


# Make dwi for mask generation purposes.
tmp_mask="${work_dir}/${id}_tmp_mask.${ext}";
#mask="${work_dir}/${id}_mask.${ext}";
raw_dwi="${work_dir}/${id}_raw_dwi.nii.gz";
if [[ ! -f ${mask} ]] && [[ ! -f ${tmp_mask} ]];then
   
    if [[ ! -f ${raw_dwi} ]];then
	   echo  select_dwi_vols ${raw_nii} $bvals ${raw_dwi} ${nominal_bval}  -m;
	   select_dwi_vols ${raw_nii} $bvals ${raw_dwi} ${nominal_bval}  -m;
    fi
    
    if [[ ! -f ${tmp_mask} ]];then
	bet ${raw_dwi} ${tmp_mask/_mask/} -m -n;    
    fi
fi


if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_mask} ]] && [[ -f ${raw_dwi} ]];then
    rm ${raw_dwi};
fi


# Run Local PCA Denoising algorithm on 4D nifti:
# Note: the LPCA python has a logical switch that will produce another, confounding mask (though arguably better). This is currently switched off, but can cause confusion and delay if switched on and not properly documented.
denoised_nii="${work_dir}/LPCA_${id}_nii4D.nii.gz";
masked_nii="${work_dir}/${nii_name}";
masked_nii=${masked_nii/.$ext/_masked.$ext};

if [[ ! -f ${denoised_nii} ]];then

    if [[ ! -f ${masked_nii} ]];then
	fslmaths ${raw_nii} -mas ${tmp_mask} ${masked_nii} -odt "input";
    fi

    /Users/alex/bass/gitfolder/gunnies/basic_LPCA_denoise.py ${id} ${masked_nii} ${bvecs} ${work_dir};
fi

if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${denoised_nii} ]] && [[ -f ${masked_nii} ]];then
    rm ${masked_nii};
fi

# Run coregistration/eddy current correction:

L_id="LPCA_${id}";
coreg_nii="${BIGGUS_DISKUS}/co_reg_${L_id}_m00-results/Reg_${L_id}_nii4D.${ext}";
if [[ ! -f ${coreg_nii} ]];then
    /Users/alex/bass/gitfolder/gunnies/co_reg_4d_stack.bash ${denoised_nii} ${L_id} 0;
fi

coreg_inputs="${BIGGUS_DISKUS}/co_reg_${L_id}_m00-inputs";
coreg_work=${coreg_inputs/-inputs/-work};
if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${coreg_nii} ]] && [[ -d ${coreg_inputs} ]];then
    rm -r ${coreg_inputs};
fi

if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${coreg_nii} ]] && [[ -d ${coreg_work} ]];then
    rm -r ${coreg_work};
fi

# Generate tmp DWI:
tmp_dwi_out=${work_dir}/${id}_tmp_dwi.${ext};
dwi_out=${work_dir}/${id}_dwi.${ext};
if [[ ! -f $dwi_out ]];then
    if [[ ! -f $tmp_dwi_out ]];then
	select_dwi_vols ${coreg_nii} $bvals ${tmp_dwi_out} ${nominal_bval}  -m;
    fi
elif [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_dwi_out} ]];then
    rm ${tmp_dwi_out};
fi

# Generate tmp B0:
tmp_b0_out=${work_dir}/${id}_tmp_b0.${ext};
b0_out=${work_dir}/${id}_b0.${ext};
if [[ ! -f ${b0_out} ]];then
    if [[ ! -f ${tmp_b0_out} ]];then
	select_dwi_vols ${coreg_nii} $bvals ${tmp_b0_out} 0  -m;
    fi
elif [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_b0_out} ]];then
    rm ${tmp_b0_out};
fi

# Generate DTI contrasts and perform some tracking QA:
c_string='';
if [[ "x${cleanup}x" == "x1x" ]];then
    c_string=' --cleanup ';
fi

/Users/alex/bass/gitfolder/gunnies/dti_qa_with_dsi_studio.bash ${coreg_nii} $bvecs ${tmp_mask} $work_dir $c_string;

# Create RELATIVE links to DSI_studio outputs:
for contrast in fa0 rd ad md; do
    real_file=$(ls ${work_dir}/*.fib.gz.${contrast}.${ext} | head -1 ); # It will be fun times if we ever have more than one match to this pattern...
    contrast=${contrast/0/};
    linked_file="${work_dir}/${id}_${contrast}.${ext}";

    # We need to make sure that we didn't accidentally link to a non-existent file.
    if [[ -L ${linked_file} ]] && [[ ! -f $(readlink ${linked_file}) ]];then
        unlink ${linked_file};
    fi

    # This should prevent linking to a non-existent file in the future.
    if [[ ! -L ${linked_file} ]] && [[ -f ${real_file} ]];then
	ln -s "./${real_file##*/}" ${linked_file};
    fi
done


# For Whitson data (at least) we empirically know that we need to move the data from LAS to RAS:
# Enforce header consistency, based on fa from DSI Studio before reorienting data (will break if there's an affine xform in header)

# Hopefully the new auto-detect-orientation code will work...

img_xform_exec='/Users/alex/MATLAB/MATLAB_BJ/matlab_execs_for_SAMBA//img_transform_executable/run_img_transform_exec.sh';
mat_library='/Users/alex/MATLAB/MATLAB_BJ/MATLAB2015b_runtime/v90';

# Enforce header consistency, based on mask--no, scratch that--based on fa from DSI Studio
# No scratch that--fa doesn't necessarily have a consistent center of mass!
md="${work_dir}/${id}_md.${ext}";

# I'm not going to say specifically why, but, we only need to compare orientations between the mask and the md.
# This only works on clean runs, and won't catch b0s or dwis that have been previously incorrectly turned.

if [[ -f $tmp_mask ]];then

    for contrast in tmp_dwi tmp_b0 tmp_mask;do
	file=${work_dir}/${id}_${contrast}.${ext};
	if [[ -e $file ]];then
	    /Users/alex/bass/gitfolder/gunnies/nifti_header_splicer.bash ${md} ${file} ${file};
	fi
    
    done

    file=${work_dir}/${id}_tmp_mask.${ext};
    orient_test=$(/Users/alex/bass/gitfolder/gunnies/find_relative_orientation_by_CoM.bash $md $file);
    
    orientation_in=$(echo ${orient_test} | cut -d ',' -f2 | cut -d ':' -f2);
    orientation_out=$(echo ${orient_test} | cut -d ',' -f1 | cut -d ':' -f2);
    echo "flexible orientation: ${orientation_in}";
    echo "reference orientation: ${orientation_out}";
    
    for contrast in dwi b0 mask;do
	img_in=${work_dir}/${id}_tmp_${contrast}.${ext};
	img_out=${work_dir}/${id}_${contrast}.${ext};

	if [[ "${orientation_out}" != "${orientation_in}" ]];then
	    echo "TRYING TO REORIENT...b0 and dwi and mask";
	    if [[ -e $img_in ]] && [[ ! -e $img_out ]];then 
		reorient_cmd="${img_xform_exec} ${mat_library} ${img_in} ${orientation_in} ${orientation_out} ${img_out}";
		${reorient_cmd};
		if [[ -e ${img_out} ]];then
		    rm $img_in;
		fi
	    elif [[ -e ${img_out} ]];then
		rm $img_in;
	    fi
	    
	else
	    mv $img_in $img_out;
	fi
    done
fi

# Turning off the following code block for now...
if (( 0 ));then
    file=${work_dir}/${id}_mask.${ext};
    orient_test=$(/Users/alex/bass/gitfolder/gunnies/find_relative_orientation_by_CoM.bash $md $file);

    orientation_in=$(echo ${orient_test} | cut -d ',' -f2 | cut -d ':' -f2);
    orientation_out=$(echo ${orient_test} | cut -d ',' -f1 | cut -d ':' -f2);
    echo "flexible orientation: ${orientation_in}";
    echo "reference orientation: ${orientation_out}";
    if [[ "${orientation_out}" != "${orientation_in}" ]];then
	echo "TRYING TO REORIENT...MASK!";
	img_in=${work_dir}/${id}_tmp2_mask.${ext};
	img_out=${work_dir}/${id}_mask.${ext};
    
	mv $img_out $img_in;
	reorient_cmd="${img_xform_exec} ${mat_library} ${img_in} ${orientation_in} ${orientation_out} ${img_out}";
	${reorient_cmd};
	if [[ -e ${img_out} ]];then
	    rm $img_in;
	fi
    fi
fi

for contrast in dwi b0 mask;do
    file=${work_dir}/${id}_${contrast}.${ext};
    /Users/alex/bass/gitfolder/gunnies/nifti_header_splicer.bash ${md} ${file} ${file};
    
done


# Apply updated mask to dwi and b0.
mask=${work_dir}/${id}_mask.${ext};

#for contrast in dwi b0;do
#    file=${work_dir}/${id}_${contrast}.${ext};
#    fslmaths ${file} -mas ${mask} ${file} -odt "input";
#done

if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_mask} ]] &&  [[ -f ${mask} ]];then
    rm ${tmp_mask};
fi


if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_dwi_out} ]] &&  [[ -f ${dwi_out} ]];then
    rm ${tmp_dwi_out};
fi

if [[ "x${cleanup}x" == "x1x" ]] && [[ -f ${tmp_b0_out} ]] &&  [[ -f ${b0_out} ]];then
    rm ${tmp_b0_out};
fi
