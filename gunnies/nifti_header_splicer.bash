#! /bin/bash

## Function:
# This script takes the spatial components from the header of a reference nifti image,
# and splices them directly into the header of a nifti that is aligned to the reference in 
# the space of the data array (i.e. they would match in the first three dimensions
# if they were loaded in MATLAB without any header info.
#
# Reference: https://brainder.org/2012/09/23/the-nifti-file-format/

## Usage:
# ./nifti_header_splicer.bash $reference_nifti $data_array_nifti [$optional_output_nifti]

## Notes:
# The size of the reference nifti has no impact on the speed of this function, whether or not it is compressed (gzip'd), as it is only
# decompressed just enough to read the header.
# If the data_array nifti is uncompressed AND the same as the output nifti, then the entire script is a nearly instantaneous operation,
# irregardless of the size of the nifti.
# If the data_array nifti is large (~1GB+) and compressed, this may take an unreasonably long amount of time as it needs to 
# be decompressed and then recompressed on the disk.
# If the output file is different than the input data_array, it will first be copied in its original format (compressed or uncomppressed),
# so large files might experience noticable overhead from the copy operation (in addition to decompression/recompression if needed).  

## Requirements:
# This script implements the following common linux/unix functions, and will fail if not available: dd, xxd, cp, gzip, gunzip, zless, and fprintf
# It is assumed that all these are present.
# dd is slightly different between unix and linux, and this is tested for:
stat='';
if [[ $(echo HI | dd status=none  2>/dev/null | grep 'HI') ]];then
    stat=" status=none ";
fi

## Parse inputs:
ref=$1;
data=$2;
output=$3;

## Set output and copy_data flag:
copy_data=1;

if [[ "x${output}x" == "xx" ]];then
    output=$data;
fi

if [[ "x${data}x" == "x${output}x" ]];then
    copy_data=0;
fi

o_data=$data;

## Set input_gz flag:
input_gz=1;
if (( $(gzip -t $data 2>&1 | wc -l) > 0 ));then
    input_gz=0;
fi

## Set output_gz flag:
output_gz=0;
if [[ "${output##*.}" == "gz" ]];then
    output_gz=1;
fi

## Test for "large" uncompressed data array and throw performance warning if need be:
test=0;

if (( $output_gz )) || ((  $copy_data == 0 ));then
    if (( $input_gz ));then
	test=$(gzip -lv $data | cut -d '%' -f1 );test=${test%\ *};test=${test%%\ };test=${test##*\ };
    elif (( $output_gz ));then
	test=$(ls -s --block-size=1 $data | cut -d ' ' -f1);
    fi

    limit=1000000000; # Arbitrarily set to ~1GBi (but can be changed to your liking here)
    human_readable_limit="1 GBi";
    if (( $test > $limit ));then 
	echo "WARNING: Your uncompressed data size is greater than ${human_readable_limit}; this process may take an inconveniently amount of time to complete!";
	echo "         You may want to consider aborting.";
	sleep 5;
    fi
fi


## Create temporary header file from reference:
tmp_header="${ref}.tmp.hdr";
zless ${ref} | head -c 352 > ${tmp_header};

## Copy input data if needed:
if (( $copy_data ));then
    copied_file=$output;
    if (( $input_gz )) && (( $output_gz == 0 )); then
	copied_file="${copied_file}.gz";
    fi

    cp $data $copied_file;

    data=$copied_file;

fi

## Decompress (gunzip) input data if needed:
if (( $input_gz ));then
    gunzip $data;
    data=${data%.gz};
fi

## Run dd commands to splice in info from the ref header:
common_params=" if=${tmp_header} of=${data} conv=notrunc $stat"

# Block 1: The first 40 bytes, no offset. Fields: sizeof_hdr to dim_info.
dd $common_params skip=0 ibs=40 count=1 obs=40 seek=0;

# Block 2: 18 bytes with 74 byte offset. Fields: slice_start and the first 8 bytes of pix_dim.
dd $common_params ibs=2 skip=37 count=9 obs=2 seek=37;

# Block 3: 3 bytes with 120 byte offset. Fields: slice_end & slice_code.
dd $common_params ibs=3 skip=40 count=1 obs=3 seek=40;

# Block 4: The first 5 BITS of xyzt_units from data, and last 3 BITS from ref.
# In this case we need to split our dd command into input and output, with some xxd in between. (Offset of 123 bytes)
p1=$(dd if=${tmp_header} bs=1 skip=123 count=1 $stat | xxd -b | cut -d ' ' -f2 | head --bytes=5);
p2=$(dd if=${data} bs=1 skip=123 count=1 $stat | xxd -b | cut -d ' ' -f2 | tail --bytes=4);
printf '%x\n' "$((2#${p1}${p2}))" | dd of=${data} bs=1 skip=123 count=1 conv=notrunc ${stat};

# Block 5: 76 bytes with 252 byte offset. Fields: qform_code thru srow_z.
dd $common_params ibs=4 skip=63 count=19 obs=4 seek=63;

## Compress/gzip output if needed:
if (( $output_gz ));then
    gzip $data;
fi

## Clean up tmp header:
rm ${tmp_header};

## Pat yourself on the mastersciecing back:
echo "Spatial header info successfully copied from ${ref} onto image data from ${o_data}."
echo "Final result saved in: $output."

# This last line could cause issues when trying to debug, since we're assuming that our code
# works correctly and the end product is named the same as the specified output (in reality,
# the result is the final value of $data, with an optional '.gz'.
