from dipy.io.streamline import load_trk
import os
import shutil
from tract_handler import trk_fixer
from file_tools import mkcdir
from tract_handler import reducetractnumber, ratio_to_str, get_ratio
import random

TRK_folder = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT'
TRK_output = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT_fixed'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT'
TRK_output = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT_fixed'

TRK_folder = '/mnt/munin6/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT/'
TRK_output = '/mnt/munin6/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed_refixed/'

filelist = os.listdir(TRK_folder)
filelist = sorted(filelist)
#filelist.reverse()
random.shuffle(filelist)

mkcdir(TRK_output)
for trk_name in filelist:
    if trk_name.endswith('.trk'):
        trk_path = os.path.join(TRK_folder,trk_name)
        trk_newpath = os.path.join(TRK_output, trk_name)
        if not os.path.exists(trk_newpath):
            print(f'Beginning the fix of {trk_path}')
            trk_fixer(trk_path, trk_newpath, verbose=True)
        else:
            print(f'Already fixed {trk_name}')
