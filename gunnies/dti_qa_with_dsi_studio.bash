#! /bin/bash
echo "Attempting to run DTI QA with DSI Studio (via command line, nonetheless!)"
nii4D=$1;
btable=$2; # This can also be a bvec OR a bval file, but the unspecified one needs have an otherwise identical path.
c_mask=$3;
out_folder=$4; # Optional; otherwise will use same folder as input nii4D.
cleanup=$5; # Optional; '--cleanup' : cleans up all DSI Studio files except stats, fa0, rd, ad, md


btest=${btable%\.txt};
btest=$(echo $btest | grep -cE b?tab);
b_string=" --b_table=${btable} ";

if [[ "xx${btest}xx" == "xx0xx" ]];then
    bvecs=$btable;
    bvals=${bvecs/bvec/bval};
    bvecs=${bvals/bval/bvec};
    b_string="--bvec=${bvecs}  --bval=${bvals} ";
fi

if [[ ! -d "${out_folder}" ]];then
    if [[ "x${out_folder}x" == 'x--cleanupx' ]];then
	cleanup='--cleanup';
    fi   
    out_folder=${nii4D%\/*}/;
else
    out_folder=${out_folder%\/}/;
fi

if [[ "x${cleanup}x" == 'x--cleanupx' ]];then
    cleanup=1;
else
    cleanup=0;
fi

echo "Cleanup = ${cleanup}";

prefix=${nii4D##*\/};
pid="9A99993E9A99193Fba3FCDCCCC3D9A99993Eb2041b96434aD03dcaCDCC4C3Ec"; #DSI Studio tracking parameters.

out_src="${out_folder}${prefix}.src.gz";
fib_file="${out_src}.dti.fib.gz";
fa_file="${fib_file}.fa0.nii.gz";
out_trk="${fib_file}_WB_0.3fa_thresh.trk.gz"; #WB => whole brain, 0.3fa_thresh => Using an fa threshold of 0.3
stat_file="${out_trk}.stat.txt";

echo "fa_file = ${fa_file}";


src_cmd='';
c_cmd='';
exp_cmd='';
t_cmd='';


if [[ ! -f "${stat_file}" ]];then
    if [[ ! -f "${out_trk}" ]];then
	t_cmd="dsi_studio --action=trk --source=${fib_file} --export=stat --output=${out_trk} --parameter_id=${pid}";
    fi 
fi

if [[ ! -f "${fa_file}" ]];then
    exp_cmd="dsi_studio --action=exp --source=${fib_file} --export=fa0,md,ad,rd";

    if [[ ! -f "${fib_file}" ]];then
	c_cmd="dsi_studio --action=rec --source=${out_src} --method=1 --output_dif=0 --thread_count=32 --mask=${c_mask} --check_btable=0";

	if [[ ! -f "${out_src}" ]];then
	    src_cmd="dsi_studio --action=src --source=${nii4D}  ${b_string} --output=${out_src}";
	fi
    fi
fi



if [[ ! -f "${fib_file}" ]] && [[ -f "${fa_file}" ]] && (( ! ${cleanup} ));then
    c_cmd="dsi_studio --action=rec --source=${out_src} --method=1 --output_dif=0 --thread_count=32 --mask=${c_mask} --check_btable=0";

    if [[ ! -f "${out_src}" ]];then
	src_cmd="dsi_studio --action=src --source=${nii4D}  ${b_string} --output=${out_src}";
    fi

fi

if [[ ! -f "${out_src}" ]] && (( ! ${cleanup} ));then
    src_cmd="dsi_studio --action=src --source=${nii4D} ${b_string}  --output=${out_src}";
fi



$src_cmd;
$c_cmd;
$exp_cmd;
$t_cmd;


if [[ -f $fa_file ]] && [[ -f $stat_file ]] && (( $cleanup ));then

    if [[ -f $out_src ]] || [[ -f $fib_file ]] || [[ -f $out_trk ]];then
	echo "Cleaning up large DSI Studio files:";

	if [[ -f $out_src ]];then
	    echo "     ${out_src}";
	    rm $out_src;
	fi
	if [[ -f $fib_file ]];then
	    echo "     ${fib_file}";
	    rm $fib_file;
	fi
	if [[ -f $out_trk ]];then
	    echo "    ${out_trk}";
	    rm $out_trk;
	fi
    fi
fi
