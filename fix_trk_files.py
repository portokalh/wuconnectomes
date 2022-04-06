from dipy.io.streamline import load_trk
import os
import shutil
from tract_handler import trk_fixer
from file_tools import mkcdir, getfromfile
from tract_handler import reducetractnumber, ratio_to_str, get_ratio
import random
from computer_nav import get_sftp, glob_remote
"""
TRK_folder = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT'
TRK_output = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT_fixed'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT'
TRK_output = '/Volumes/Data/Badea/Lab/human/AMD/TRK_MDT_fixed'

TRK_folder = '/mnt/munin6/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT/'
TRK_output = '/mnt/munin6/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Centroids_MDT_non_inclusive_100/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Centroids_MDT_non_inclusive_100_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed/'

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_fixed/'
"""

TRK_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_/'
TRK_output = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_fixed/'

remote=True
if remote:
    username, passwd = getfromfile('/Users/jas/samos_connect.rtf')

inpath, outpath, atlas_folder, sftp = get_sftp(remote, username=username,password=passwd)

if not remote:
    filelist = os.listdir(TRK_folder)
else:
    filelist = glob_remote(TRK_folder,sftp)

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
