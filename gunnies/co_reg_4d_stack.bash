#!/bin/bash

## RJA, Badea Lab, 22 April 2020
# This function uses antsRegistration from the ANTs toolbox to affinely register all images in a 4D nifti stack to the first volume.
# In the '-work' folder are the individual images after the calculated transforms have been applied to the corresponding input image.
# These images are then reconcatonated into a registered 4D nift stack and placed in the '-results' folder. 
# Also in the '-results' folder are the individual transforms for each volume.
# Currently these transforms are NOT being applied to any DWI b-vector table.

# Though it doesn't belong here on a long-term basis, some minimal DWI processing has been included here.
# Right now, this will extract the b-table from a .bxh header, if there is one immediately next to the 4D nifti stack, and has the same name, save for the extension (.bxh instead of .nii or .nii.gz)

# There is hope that it will also spit out the dwi image contrast while it's at it...but that's a WIP.

nii4D=$1;
identifier=$2
dti=$3

if [[ ! -f ${nii4D} ]];then
    nii4D=${PWD}/${nii4D};
fi

job_desc="co_reg"; #e.g. co_reg
job_shorthand="Reg";#"Reg" for co_reg
ext="nii.gz";


sbatch_file='';


if [[ "x${identifier}x" == "xx" ]];then
    base=${nii4D##*/};
    runno=$(echo $base | cut -d '.' -f1);
else
    runno=$identifier;
fi

echo "Processing runno: ${runno}";

if [[ ! -f $nii4D ]];then    
    echo "ABORTING: Input file does not exist: ${nii4D}" && exit 1;
fi

YYY=$(PrintHeader $nii4D 2 | cut -d 'x' -f4);
XXX=$(expr $YYY - 1);

declare -i XXX; # Need to 

echo "Total number of volumes: $YYY";
echo "Number of independently oriented volumes: $XXX";

if [[ $XXX -lt 10 ]];then
    zeros='0';
elif [[ $XXX -lt 100 ]];then
    zeros='00';
else
    zeros='000';
fi

discard='1000';
discard=${discard/${zeros}/};
zero='0';
zero_pad=${zeros/${zero}/};

inputs="${BIGGUS_DISKUS}/${job_desc}_${runno}_m${zeros}-inputs/";
work="${BIGGUS_DISKUS}/${job_desc}_${runno}_m${zeros}-work/";
results="${BIGGUS_DISKUS}/${job_desc}_${runno}_m${zeros}-results/";

vol_zero="${inputs}${runno}_m${zeros}.${ext}";

echo "Target for coregistration: ${vol_zero}";


if [[ ! -d ${inputs} ]];then
    mkdir -p -m 775 ${inputs};
fi

if [[ ! -d ${work} ]];then
    mkdir -p -m 775 ${work};
fi

sbatch_folder="${work}/sbatch/";
if [[ ! -d ${sbatch_folder} ]];then
    mkdir -p -m 775 ${sbatch_folder};
fi

if [[ ! -d ${results} ]];then
    mkdir -p -m 775 $results;
fi

prefix="${BIGGUS_DISKUS}/${job_desc}_${runno}_m${zeros}-inputs/${runno}_m.nii.gz";

## To-do (2 June 2020 Tues): Make this a cluster job; use for later jobs: echo "#\$ -hold_jid ${jid_list}" >> ${sbatch_file};

if [[ ! -e ${vol_zero} ]];then
    if [[ ! -e ${prefix/_m/_m1000} ]];then
	echo "Splitting up nii4D volume...";
	${ANTSPATH}/ImageMath 4 ${prefix} TimeSeriesDisassemble ${nii4D};
    fi

    for file in $(ls ${inputs});do
	new_file=${inputs}/${file/_m${discard}/_m};
	if [[ ! -e ${new_file} ]];then
	    ln -s $file ${new_file};
	fi
    done
fi

work_vol_zero="${work}/${job_shorthand}_${runno}_m${zeros}.${ext}";
reassemble_list="${work_vol_zero} ";
jid_list='';


#for nn in $(seq 1 $XXX);do

if [[ ! -e ${work_vol_zero} ]];then
    ln -s ${vol_zero} ${work_vol_zero};
fi
echo "Dispatching co-registration jobs to the cluster:";

# Note the following line is necessarily complicated, as...
# the common sense line ('for nn in {01..$XXX}') does not work...
# https://stackoverflow.com/questions/169511/how-do-i-iterate-over-a-range-of-numbers-defined-by-variables-in-bash
for nn in $(eval echo "{${zero_pad}1..$XXX}");do
    num_string=$nn;
    if [[ $nn -lt 10 ]];then
        num_string="0${nn}";
    fi
   	#num_string=$nn;
	vol_xxx="${inputs}${runno}_m${num_string}.${ext}";
	out_prefix="${results}xform_${runno}_m${num_string}.${ext}";
	xform_xxx="${out_prefix}0GenericAffine.mat";
	vol_xxx_out="${work}/${job_shorthand}_${runno}_m${num_string}.${ext}";
	reassemble_list="${reassemble_list} ${vol_xxx_out} ";

	name="${job_desc}_${runno}_m${num_string}";
	sbatch_file="${sbatch_folder}/${name}.bash";
	#source_sbatch="${BIGGUS_DISKUS}/sinha_co_reg_nii4D_qsub_master.bash";
	#cp ${source_sbatch} ${sbatch_file};
	if [[ ! -e ${xform_xxx} ]] ||  [[ ! -e ${vol_xxx_out} ]];then

	    echo "#!/bin/bash" > ${sbatch_file};
	    echo "#\$ -l h_vmem=8000M,vf=8000M" >> ${sbatch_file};
	    echo "#\$ -M ${USER}@duke.edu" >> ${sbatch_file};
	    echo "#\$ -m ea" >> ${sbatch_file}; 
	    echo "#\$ -o ${sbatch_folder}"'/slurm-$JOB_ID.out' >> ${sbatch_file};
	    echo "#\$ -e ${sbatch_folder}"'/slurm-$JOB_ID.out' >> ${sbatch_file};
	    echo "#\$ -N ${name}" >> ${sbatch_file};

	    reg_cmd="if [[ ! -e ${xform_xxx} ]];then ${ANTSPATH}/antsRegistration  --float -d 3 -v  -m Mattes[ ${vol_zero},${vol_xxx},1,32,regular,0.3 ] -t Affine[0.05] -c [ 100x100x100,1.e-5,15 ] -s 0x0x0vox -f 4x2x1 -u 1 -z 1 -o ${out_prefix};fi";
	    apply_cmd="if [[ ! -e ${vol_xxx_out} ]];then ${ANTSPATH}/antsApplyTransforms -d 3 -e 0 -i ${vol_xxx} -r ${vol_zero} -o ${vol_xxx_out} -n Linear -t ${xform_xxx}  -v 0 --float;fi";

	    echo "${reg_cmd}" >> ${sbatch_file};
	    echo "${apply_cmd}" >> ${sbatch_file};

	    if ! command -v qsub &> /dev/null; then
	    	bash ${sbatch_file};
	    else
	    	# cmd="qsub -terse  -b y -V ${sbatch_file}";
	    	cmd="qsub -terse -V ${sbatch_file}";
	    	echo $cmd;
	    	job_id=$($cmd | tail -1);
	    	echo "JOB ID = ${job_id}; Job Name = ${name}";
		new_sbatch_file=${sbatch_file/${name}/${job_id}_${name}};
		mv ${sbatch_file} ${new_sbatch_file};
	    	jid_list="${jid_list}${job_id},";
	    fi
	fi
done

#Trim trailing comma from job id list:
jid_list=${jid_list%,};

reg_nii4D="${results}/${job_shorthand}_${runno}_nii4D.${ext}";
assemble_cmd="${ANTSPATH}/ImageMath 4 ${reg_nii4D} TimeSeriesAssemble 1 0 ${reassemble_list}";
#if [[ 1 -eq 2 ]];then # Uncomment when we want to short-circuit this to OFF
if [[ ! -f ${reg_nii4D} ]];then
    name="assemble_nii4D_${job_desc}_${runno}_m${zeros}";
    sbatch_file="${sbatch_folder}/${name}.bash";
    
    echo "#!/bin/bash" > ${sbatch_file};
    echo "#\$ -l h_vmem=32000M,vf=32000M" >> ${sbatch_file};
    echo "#\$ -M ${USER}@duke.edu" >> ${sbatch_file};
    echo "#\$ -m ea" >> ${sbatch_file}; 
    echo "#\$ -o ${sbatch_folder}"'/slurm-$JOB_ID.out' >> ${sbatch_file};
    echo "#\$ -e ${sbatch_folder}"'/slurm-$JOB_ID.out' >> ${sbatch_file};
    echo "#\$ -N ${name}" >> ${sbatch_file};
    if [[ "x${jid_list}x" != "xx" ]];then
	echo "#\$ -hold_jid ${jid_list}" >> ${sbatch_file};
    fi
    echo "${assemble_cmd}" >> ${sbatch_file};

    if ! command -v qsub &> /dev/null; then
        bash ${sbatch_file};
    else
        ass_cmd="qsub -terse -V ${sbatch_file}";
        echo $ass_cmd;
        job_id=$($ass_cmd | tail -1);
        echo "JOB ID = ${job_id}; Job Name = ${name}";
        new_sbatch_file=${sbatch_file/${name}/${job_id}_${name}};
        mv ${sbatch_file} ${new_sbatch_file};
    fi
fi 

if [[ "x${dti}x" == "x1x" ]];then
    bvecs=${reg_nii4D/\.${ext}/_fsl_bvecs.txt};
    bvecs=${bvecs/${job_shorthand}_/};
    bvals=${bvecs/bvecs/bvals};

    dsi_btable=${bvecs/fsl/dsi};
    dsi_btable=${dsi_btable/bvecs/btable};

    if [[ ! -f ${bvecs} ]];then
	bxh=${nii4D/${ext}/bxh};
	bvec_cmd="extractdiffdirs --fsl ${bxh} ${bvecs} ${bvals}";
	dsi_bvec_cmd="extractdiffdirs --dsistudio ${bxh} ${dsi_btable}";
	echo $bvec_cmd;
	$bvec_cmd;

	echo $dsi_bvec_cmd;
	$dsi_bvec_cmd;
    fi

    

fi
