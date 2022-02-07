from dipy.io.streamline import load_trk
import os
import shutil
from tract_handler import trk_fixer
from file_tools import mkcdir
from tract_handler import reducetractnumber, ratio_to_str, get_ratio

TRK_folder = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT'
TRK_output = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT_fixed'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT'
TRK_output = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT_fixed'

TRK_folder = '/mnt/paros_MRI/jacques/AD_Decode/Analysis/TRK_MPCA_neworient/'
TRK_output = '/mnt/paros_MRI/jacques/AD_Decode/Analysis/TRK_MPCA_neworient_fixed/'

mkcdir(TRK_output)
for trk_name in os.listdir(TRK_folder):
    if trk_name.endswith('.trk'):
        trk_path = os.path.join(TRK_folder,trk_name)
        trk_newpath = os.path.join(TRK_output, trk_name)
        if not os.path.exists(trk_newpath):
            print(f'Beginning the fix of {trk_path}')
            trk_fixer(trk_path, trk_newpath, verbose=True)
        else:
            print(f'Already fixed {trk_name}')

trk_folder = TRK_output
new_trk_folder = '/mnt/paros_MRI/jacques/AD_Decode/Analysis/TRK_MPCA_neworient_fixed_100/'
mkcdir(new_trk_folder)
ratio = 100
filelist = os.listdir(trk_folder)
#filelist = sorted(filelist)
filelist.reverse()
for filepath in filelist:
    filename, f_ext = os.path.splitext(filepath)
    ratio1 =get_ratio(filepath)
    ratio2 = int(ratio1*ratio)
    filename2 = filename.replace(ratio_to_str(ratio1),ratio_to_str(ratio2))
    newfilepath = os.path.join(new_trk_folder, filename2 + f_ext)
    if not os.path.exists(newfilepath):
        print(f'Downsampling from {os.path.join(trk_folder,filepath)} to {newfilepath} with ratio {str(ratio)}')
        #print(f'{filename} converted to {filename2 + f_ext}')
        reducetractnumber(os.path.join(trk_folder,filepath), newfilepath, getdata=False, ratio=100, return_affine= False, verbose=False)
        print(f'succesfully created {newfilepath}')
    else:
        print(f'already created {newfilepath}')

"""
if testing:
    trk_test = '/Users/alex/jacques/AMD_TRK_testing/TRK/H21850_stepsize_2_ratio_100_wholebrain_pruned.trk'
    trk_new = '/Users/alex/jacques/AMD_TRK_testing/TRK/H21850_ratio100_fixed.trk'
    trk_fixer(trk_test, trk_new, verbose=True)
else:
"""
