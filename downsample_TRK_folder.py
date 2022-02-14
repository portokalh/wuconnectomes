import os
from tract_handler import reducetractnumber
from file_tools import mkcdir

trk_folder = '/mnt/paros_MRI/jacques/AMD/TRK_MDT_fixed'
new_trk_folder = '/mnt/paros_MRI/jacques/AMD/TRK_MDT_fixed_100'
trk_folder = '/Volumes/Data/Badea/Lab/jacques/ADDecode_test_connectomes/TRK_MPCA_fixed'
new_trk_folder = '/Volumes/Data/Badea/Lab/jacques/ADDecode_test_connectomes/TRK_MPCA_fixed_100'

trk_folder = '/Users/alex/jacques/AD_decode_comparison_temp/'
new_trk_folder = '/Users/alex/jacques/AD_decode_comparison_100_temp/'

trk_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA'
new_trk_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_100'

trk_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed'
new_trk_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT_fixed_100'

mkcdir(new_trk_folder)
ratio = 100
filelist = os.listdir(trk_folder)
filelist = sorted(filelist)
filelist.reverse()

for filepath in filelist:
    filename, f_ext = os.path.splitext(filepath)
    if f'ratio_{str(ratio)}' not in filename:
        newfilepath = os.path.join(new_trk_folder, filename + f'_ratio_{str(ratio)}'+f_ext)
        if not os.path.exists(newfilepath):
            print(f'Downsampling from {os.path.join(trk_folder,filepath)} to {newfilepath} with ratio {str(ratio)}')
            reducetractnumber(os.path.join(trk_folder,filepath), newfilepath, getdata=False, ratio=100, return_affine= False, verbose=False)
            print(f'succesfully created {newfilepath}')
        else:
            print(f'already created {newfilepath}')
