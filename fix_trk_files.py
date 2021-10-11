from dipy.io.streamline import load_trk
import os
import shutil
from tract_handler import trk_fixer
from file_tools import mkcdir

TRK_folder = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT'
TRK_output = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT_fixed'

hostname = 'samos'
if hostname == 'samos':
    TRK_folder = '/mnt/paros_MRI/jacques/AMD/TRK_MDT/'
    TRK_output = '/mnt/paros_MRI/jacques/AMD/TRK_MDT_fixed'

mkcdir(TRK_output)
for trk_name in os.listdir(TRK_folder):
    if trk_name.endswith('.trk'):
        trk_path = os.path.join(TRK_folder,trk_name)
        trk_newpath = os.path.join(TRK_output, trk_name)
        trk_fixer(trk_path, trk_newpath, verbose=True)

"""
if testing:
    trk_test = '/Users/alex/jacques/AMD_TRK_testing/TRK/H21850_stepsize_2_ratio_100_wholebrain_pruned.trk'
    trk_new = '/Users/alex/jacques/AMD_TRK_testing/TRK/H21850_ratio100_fixed.trk'
    trk_fixer(trk_test, trk_new, verbose=True)
else:
"""