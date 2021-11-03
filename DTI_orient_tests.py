import numpy as np
from bvec_handler import writebvec, find_bval_bvecs, orient_to_str, reorient_bvecs, read_bval_bvec
from file_tools import mkcdir
import itertools
import os

mainpath = '/Users/alex/jacques/APOE/bvec_testing/'
input_diifile = '/Users/alex/jacques/APOE/bvec_testing/N58794_masked.nii.gz'
mask_path = '/Users/alex/jacques/APOE/bvec_testing/N58794_masked_mask.nii.gz'
bval_path = '/Users/alex/jacques/APOE/bvec_testing/N58794_bvals_fix.txt'
bvec_path = '/Users/alex/jacques/APOE/bvec_testing/N58794_bvecs_fix.txt'
bvals, bvecs = read_bval_bvec(bval_path, bvec_path)

bvec_orient1 = (np.array(list(itertools.permutations([1, 2, 3]))))
bvec_orient2 = [elm*[-1, 1, 1] for elm in bvec_orient1]
bvec_orient3 = [elm*[1, -1, 1] for elm in bvec_orient1]
bvec_orient4 = [elm*[1, 1, -1] for elm in bvec_orient1]
bvec_orient_list = np.concatenate((bvec_orient1, bvec_orient2, bvec_orient3, bvec_orient4))

dtis='#! /bin/bash \n'
for bvec_orient in bvec_orient_list:
    strproperty = orient_to_str(bvec_orient)
    results_folder = os.path.join(mainpath, strproperty.replace('_',''))
    mkcdir(results_folder)
    bvecs_reoriented = reorient_bvecs(bvecs, bvec_orient)
    new_bvec_file = os.path.join(results_folder, strproperty.replace('_','')+'_bvecs.txt')
    writebvec(bvecs_reoriented, new_bvec_file, '',writeformat = "dsi", overwrite=False)
    output_path = os.path.join(results_folder, f'dtifit{strproperty}.nii.gz')
    dti_cmd = f'dtifit -k {input_diifile} -o {output_path} -m {mask_path} -r {new_bvec_file} -b {bval_path}\n'
    dtis = dtis + dti_cmd

bash_file_path = os.path.join(mainpath, 'dti_bash.sh')
with open (bash_file_path, 'w') as rsh:
    rsh.write(dtis)



