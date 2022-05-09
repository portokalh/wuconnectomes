import numpy as np
from dipy.io.streamline import load_trk
import warnings
from dipy.tracking.streamline import transform_streamlines
import os, glob
from nifti_handler import getlabeltypemask
from file_tools import mkcdir, getfromfile, check_files
from tract_handler import ratio_to_str, gettrkpath, gettrkpath_testsftp
from convert_atlas_mask import atlas_converter
import socket
from excel_management import M_grouping_excel_save
import sys
from argument_tools import parse_arguments_function
from connectome_handler import connectivity_matrix_custom, connectivity_matrix_func
import random
from time import time
import getpass
from connectome_handler import _to_voxel_coordinates_warning, retweak_points
from dipy.viz import window, actor
from time import sleep
from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines
from dipy.segment.clustering import ClusterCentroid
from dipy.tracking.streamline import Streamlines
from tract_visualize import show_bundles, setup_view
from tract_save import save_trk_header
from tract_save import unload_trk
import errno
import pickle
from dipy.segment.clustering import QuickBundles
from dipy.io.image import load_nifti
from computer_nav import get_mainpaths, get_atlas, load_trk_remote, checkfile_exists_remote
from streamline_nocheck import load_trk as load_trk_spe

#def get_grouping(grouping_xlsx):
#    print('not done yet')

project = 'AMD'

remote=True
if remote:
    username, passwd = getfromfile('/Users/jas/samos_connect.rtf')

inpath, outpath, atlas_folder, sftp = get_mainpaths(remote,project = project, username=username,password=passwd)

if project=='AMD' or project=='AD_Decode':
    atlas_legends = get_atlas(atlas_folder, 'IIT')

# Setting identification parameters for ratio, labeling type, etc
ratio = 1
ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio', '')

inclusive = False
symmetric = True
fixed = False
overwrite = False

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

if fixed:
    fixed_str = '_fixed'
else:
    fixed_str = ''

labeltype = 'lrordered'
verbose = True
picklesave = True

function_processes = parse_arguments_function(sys.argv)
print(f'there are {function_processes} function processes')

if project=='AD_Decode':
    outpath = os.path.join(outpath,'Analysis')
    inpath = os.path.join(inpath, 'Analysis')


TRK_folder = os.path.join(inpath, f'TRK_MPCA_MDT{fixed_str}{folder_ratio_str}')
TRK_folder = os.path.join(inpath, f'TRK_rigidaff{fixed_str}{folder_ratio_str}')
label_folder = os.path.join(inpath, 'DWI')
#trkpaths = glob.glob(os.path.join(TRK_folder, '*trk'))
excel_folder = os.path.join(outpath, f'Excels_affinerigid{inclusive_str}{symmetric_str}{folder_ratio_str}')

mkcdir(excel_folder,sftp)

if not remote and os.path.exists(TRK_folder):
    raise Exception(f'cannot find TRK folder at {TRK_folder}')

# Initializing dictionaries to be filled
stream_point = {}
stream = {}
groupstreamlines = {}
groupLines = {}
groupPoints = {}
group_qb = {}
group_clusters = {}
groups_subjects = {}

if project == 'AD_Decode':
    subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361',
                'S02363',
                'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02469',
                'S02473',
                'S02485', 'S02491', 'S02490', 'S02506', 'S02523', 'S02524', 'S02535', 'S02654', 'S02666', 'S02670',
                'S02686',
                'S02690', 'S02695', 'S02715', 'S02720', 'S02737', 'S02745', 'S02753', 'S02765', 'S02771', 'S02781',
                'S02802',
                'S02804', 'S02813', 'S02812', 'S02817', 'S02840', 'S02842', 'S02871', 'S02877', 'S02898', 'S02926',
                'S02938',
                'S02939', 'S02954', 'S02967', 'S02987', 'S03010', 'S03017', 'S03028', 'S03033', 'S03034', 'S03045',
                'S03048',
                'S03069', 'S03225', 'S03265', 'S03293', 'S03308', 'S03321', 'S03343', 'S03350', 'S03378', 'S03391',
                'S03394']
    removed_list = ['S02523']
    str_identifier = f'_stepsize_2{ratio_str}_wholebrain_pruned'

elif project == 'AMD':
    groups_subjects['testing'] = ['H22825']
    groups_subjects['Initial AMD'] = ['H27640', 'H27778', 'H29020', 'H26637', 'H27680', 'H26765', 'H27017',
                                      'H26880', 'H28308', 'H28433', 'H28338', 'H26660', 'H28809', 'H27610',
                                      'H26745', 'H27111', 'H26974', 'H27391', 'H28748', 'H29025', 'H29013',
                                      'H27381', 'H26958', 'H28662', 'H26578', 'H28698', 'H27495', 'H28861',
                                      'H28115', 'H28437', 'H26850', 'H28532', 'H28377', 'H28463', 'H26890',
                                      'H28373', 'H28857', 'H27164', 'H27982']
    groups_subjects['Paired 2-YR AMD'] = ['H22825', 'H21850', 'H29225', 'H29304', 'H29060', 'H23210', 'H21836',
                                          'H29618', 'H22644', 'H22574', 'H22369', 'H29627', 'H29056', 'H22536',
                                          'H23143', 'H22320', 'H22898', 'H22864', 'H29264', 'H22683']
    groups_subjects['Initial Control'] = ['H26949', 'H27852', 'H28029', 'H26966', 'H27126', 'H28068', 'H29161',
                                          'H28955', 'H26862', 'H28262', 'H28856', 'H27842', 'H27246', 'H27869',
                                          'H27999', 'H29127', 'H28325', 'H26841', 'H29044', 'H27719', 'H27100',
                                          'H29254', 'H27682', 'H29002', 'H29089', 'H29242', 'H27488', 'H27841',
                                          'H28820', 'H27163', 'H28869', 'H28208', 'H27686']
    groups_subjects['Paired 2-YR Control'] = ['H29403', 'H22102', 'H29502', 'H22276', 'H29878', 'H29410', 'H22331',
                                              'H22368', 'H21729', 'H29556', 'H21956', 'H22140', 'H23309', 'H22101',
                                              'H23157', 'H21593', 'H21990', 'H22228', 'H23028', 'H21915']
    groups_subjects['Paired Initial Control'] = ['H27852', 'H28029', 'H26966', 'H27126', 'H29161', 'H28955',
                                                 'H26862', 'H27842', 'H27999', 'H28325', 'H26841', 'H27719',
                                                 'H27100', 'H27682', 'H29002', 'H27488', 'H27841', 'H28820',
                                                 'H28208', 'H27686']
    groups_subjects['Paired Initial AMD'] = ['H29020', 'H26637', 'H27111', 'H26765', 'H28308', 'H28433', 'H26660',
                                             'H28182', 'H27391', 'H28748', 'H28662', 'H26578', 'H28698',
                                             'H27495', 'H28861', 'H28115', 'H28377', 'H26890', 'H28373', 'H27164']

    # groups to go through
    groups_all = ['Paired 2-YR AMD','Initial AMD','Initial Control','Paired 2-YR Control','Paired Initial Control','Paired Initial AMD']
    groups = ['Paired Initial Control', 'Paired Initial AMD']
    groups = ['Paired 2-YR Control', 'Paired 2-YR AMD']
    groups = ['Initial AMD','Initial Control']
    removed_list=[]
    # groups = ['Paired 2-YR AMD']
    # groups = ['Paired 2-YR Control']
    # groups=[groups[0]]
    subjects = []
    str_identifier = f'*'

    for group in groups:
        subjects = subjects + groups_subjects[group]

elif project == 'APOE':
    raise Exception('not implemented')
else:
    txt = f'{project} not implemented'
    raise Exception(txt)

random.shuffle(subjects)
# removed_list = ['S02266']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

_, _, index_to_struct, _ = atlas_converter(atlas_legends)
labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', atlas_legends,
                                                     labeltype=labeltype, verbose=verbose, sftp=sftp)

print(f'Beginning streamline_prep run from {TRK_folder} for folder {excel_folder}')

for subject in subjects:

    trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned = False, verbose = verbose, sftp = sftp)

    if not exists:
        txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
        warnings.warn(txt)
        continue

    M_xlsxpath = os.path.join(excel_folder, subject + "_connectomes.xlsx")
    grouping_xlsxpath = os.path.join(excel_folder, subject + "_grouping.xlsx")

    _, exists = check_files([M_xlsxpath,grouping_xlsxpath], sftp=sftp)
    if np.all(exists) and not overwrite:
        if verbose:
            print(f'Found written file for subject {subject} at {M_xlsxpath} and {grouping_xlsxpath}')
        continue
    else:
        t1 = time()

        trkdata = load_trk_remote(trkpath, 'same', sftp)

        if verbose:
            print(f"Time taken for loading the trk file {trkpath} set was {str((- t1 + time()) / 60)} minutes")
        t2 = time()
        header = trkdata.space_attributes

        streamlines_world = transform_streamlines(trkdata.streamlines, np.linalg.inv(labelaffine))
        if function_processes == 1:

            M, _, _, _, grouping = connectivity_matrix_custom(streamlines_world, np.eye(4), labelmask,
                                                              inclusive=inclusive, symmetric=symmetric,
                                                              return_mapping=True,
                                                              mapping_as_streamlines=False, reference_weighting=None,
                                                              volume_weighting=False)
        else:
            M, _, _, _, grouping = connectivity_matrix_func(streamlines_world, np.eye(4), labelmask,
                                                            inclusive=inclusive,
                                                            symmetric=symmetric, return_mapping=True,
                                                            mapping_as_streamlines=False, reference_weighting=None,
                                                            volume_weighting=False,
                                                            function_processes = function_processes, verbose=False)

        M_grouping_excel_save(M, grouping, M_xlsxpath, grouping_xlsxpath, index_to_struct, verbose=False, sftp=sftp)

        del (trkdata)
        if verbose:
            print(f"Time taken for creating this connectome was set at {str((- t2 + time()) / 60)} minutes")
        if checkfile_exists_remote(grouping_xlsxpath,sftp):
            if verbose:
                print(f'Saved grouping for subject {subject} at {grouping_xlsxpath}')
        # grouping = extract_grouping(grouping_xlsxpath, index_to_struct, np.shape(M), verbose=verbose)
        else:
            raise Exception(f'saving of the excel at {grouping_xlsxpath} did not work')
