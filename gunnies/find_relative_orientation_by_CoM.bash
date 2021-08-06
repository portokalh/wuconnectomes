#! /bin/bash
# Takes in two images, and approximates in which octant their centers of mass are. For a given or assumed orientation of the reference image, this will output the relative orientation of the the moving image.

# Note that as of 25 February 2021, this only supports dimension flips. Dimension swaps may be added at a later time.

ref=$1
moving=$2

ref_error=0;
data_error=0;

if [[ ! -e ${ref} ]]; then
    ref_error=1;
fi

if [[ ! -e ${moving} ]]; then
    data_error=1;
fi

if (( $ref_error && $data_error )); then
    echo "ERROR: Reference image does not appear to exist: ${ref}."
    echo "ERROR: Moving image does not appear to exist: ${moving}."
    exit
elif (( $ref_error )); then
    echo "ERROR: Reference image does not appear to exist: ${ref}."
    exit;
elif (( $data_error )); then
    echo "ERROR: Moving image does not appear to exist: ${moving}."
    exit;
fi



# This works much more consistently if we binarize any input images...so please don't expect this to  work. Including the -u option to turn off in case one would like to try it on regular unmasked images (it could possible work, but we make no guarantees).

binarize=1;

ref_orientation=$3;
# Will default to 'RAS' if not specified. Currently doing little checking for non-sense input.

binarize_option=$4;

if [[ ${ref_orientation} == '-u' ]];then
    binarize=0;
    ref_orientation=$binarize_option;    
fi

if [[ "x${ref_orientation}x" == "xx" ]];then
    ref_orientation='RAS';
    binarize=1;
fi

if [[ "x${binarize_option}x" == "xx" ]];then
    binarize=1;
fi


# Binarize!

if (( $binarize ));then
    #new_ref=${TMPDIR}/"bin_$(date -r ${ref} --iso-8601=seconds)_${ref##*\/}";
    new_ref=${TMPDIR}/"bin_mytempfile_${ref##*\/}";
    fslmaths $ref -bin $new_ref  -odt char;
    ref=${new_ref};

    #new_moving=${TMPDIR}/"bin_$(date -r ${moving} --iso-8601=seconds)_${moving##*\/}";
    new_moving=${TMPDIR}/"bin_mytempfile_${moving##*\/}";
    fslmaths $moving -bin $new_moving  -odt char;
    moving=${new_moving};

fi


rox=${ref_orientation:0:1}
roy=${ref_orientation:1:1}
roz=${ref_orientation:2:1}

# The base assumption is that they are the same:

mox=$rox;
moy=$roy;
moz=$roz;


# Get array size (though really they should be the same, since this only really for co-sourced data, i.e. various diffusion contrasts.
ref_size_array=$(PrintHeader $ref 2);

ref_x_size=$(echo ${ref_size_array} | cut -d 'x' -f1 );
ref_y_size=$(echo ${ref_size_array} | cut -d 'x' -f2);
ref_z_size=$(echo ${ref_size_array} | cut -d 'x' -f3);

moving_size_array=$(PrintHeader $moving 2);

moving_x_size=$(echo ${moving_size_array} | cut -d 'x' -f1);
moving_y_size=$(echo ${moving_size_array} | cut -d 'x' -f2);
moving_z_size=$(echo ${moving_size_array} | cut -d 'x' -f3);

# Now find the CoM for each:
rx=$(fslstats $ref -C -l 0 | cut -d ' ' -f1);
ry=$(fslstats $ref -C -l 0 | cut -d ' ' -f2);
rz=$(fslstats $ref -C -l 0 | cut -d ' ' -f3);

mx=$(fslstats $moving -C -l 0 | cut -d ' ' -f1);
my=$(fslstats $moving -C -l 0 | cut -d ' ' -f2);
mz=$(fslstats $moving -C -l 0 | cut -d ' ' -f3);

# Calculate center index of image, counting from zero:
rcx=$(echo "( ${ref_x_size} - 1 ) / 2" | bc -l);
rcy=$(echo "( ${ref_y_size} - 1 ) / 2" | bc -l);
rcz=$(echo "( ${ref_z_size} - 1 ) / 2" | bc -l);

mcx=$(echo "( ${moving_x_size} - 1 ) / 2" | bc -l);
mcy=$(echo "( ${moving_y_size} - 1 ) / 2" | bc -l);
mcz=$(echo "( ${moving_z_size} - 1 ) / 2" | bc -l);

####
# Flip x?
if [[ $rx > $rcx && $mx < $mcx ]] || [[ $rx < $rcx && $mx > $mcx ]];then

test_letter=$rox;
new_val='';
if [[ $test_letter == 'R' ]];then
    new_val='L';
elif [[ $test_letter == 'L' ]];then
    new_val='R';
elif [[ $test_letter == 'A' ]];then
    new_val='P';
elif [[ $test_letter == 'P' ]];then
    new_val='A';
elif [[ $test_letter == 'S' ]];then
    new_val='I';
else
 ERROROROR;
fi

mox=$new_val;
fi
# Flip y?
if [[ $ry > $rcy && $my < $mcy ]] || [[ $ry < $rcy && $my > $mcy ]];then
test_letter=$roy;
new_val='';
if [[ $test_letter == 'R' ]];then
    new_val='L';
elif [[ $test_letter == 'L' ]];then
    new_val='R';
elif [[ $test_letter == 'A' ]];then
    new_val='P';
elif [[ $test_letter == 'P' ]];then
    new_val='A';
elif [[ $test_letter == 'S' ]];then
    new_val='I';
else
 ERROROROR;
fi

moy=$new_val;
fi 
# Flip z?
if [[ $rz > $rcz && $mz < $mcz ]] || [[ $rz < $rcz && $mz > $mcz ]];then
test_letter=$roz;
new_val='';
if [[ $test_letter == 'R' ]];then
    new_val='L';
elif [[ $test_letter == 'L' ]];then
    new_val='R';
elif [[ $test_letter == 'A' ]];then
    new_val='P';
elif [[ $test_letter == 'P' ]];then
    new_val='A';
elif [[ $test_letter == 'S' ]];then
    new_val='I';
else
 ERROROROR;
fi

moz=$new_val;
fi 

echo "Ref:${rox}${roy}${roz},Moving:${mox}${moy}${moz}";

#return "${mox}${moy}${moz}";
