import numpy as np
from dipy.io.streamline import load_trk
import warnings
from dipy.tracking.streamline import transform_streamlines
import os, glob
from nifti_handler import getlabeltypemask
from file_tools import mkcdir
from tract_handler import ratio_to_str, gettrkpath
from convert_atlas_mask import atlas_converter
import socket
from excel_management import M_grouping_excel_save
import sys
from argument_tools import parse_arguments_function
from connectome_handler import connectivity_matrix_custom, connectivity_matrix_func
import random
from time import time

from nibabel.streamlines.array_sequence import ArraySequence
from dipy.segment.clustering import ClusterCentroid
from dipy.tracking.streamline import Streamlines
from tract_visualize import show_bundles, setup_view
from tract_save import save_trk_header
from tract_save import unload_trk
import pickle
from dipy.tracking.utils import connectivity_matrix
import errno
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric,mdf
from dipy.io.image import load_nifti
from dipy.viz import window, actor
from time import sleep
from dipy.tracking.streamline import set_number_of_points
from dipy.segment.clustering import QuickBundles

#def get_grouping(grouping_xlsx):
#    print('not done yet')

project = 'AD_Decode'

computer_name = socket.gethostname()

samos = False
if 'samos' in computer_name:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in computer_name or 'hydra' in computer_name:
    #mainpath = '/Users/alex/jacques/'
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    raise Exception('No other computer name yet')

#Setting identification parameters for ratio, labeling type, etc
ratio = 1
ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio','')

inclusive = True
symmetric = True
fixed = True
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

str_identifier = f'_stepsize_2{ratio_str}_wholebrain_pruned'
labeltype = 'lrordered'
verbose=True
picklesave=True

function_processes = parse_arguments_function(sys.argv)
print(f'there are {function_processes} function processes')

if project=='AD_Decode':
    mainpath=os.path.join(mainpath,project,'Analysis')
else:
    mainpath = os.path.join(mainpath, project)
TRK_folder = os.path.join(mainpath, f'TRK_MPCA_MDT{fixed_str}{folder_ratio_str}')
label_folder = os.path.join(mainpath, 'DWI')
trkpaths = glob.glob(os.path.join(TRK_folder, '*trk'))
excel_folder = os.path.join(mainpath, f'Excels_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')

print(excel_folder)
mkcdir(excel_folder)
if not os.path.exists(TRK_folder):
    raise Exception(f'cannot find TRK folder at {TRK_folder}')

#reference_img refers to statistical values that we want to compare to the streamlines, say fa, rd, etc

references = ['fa', 'md', 'rd', 'ad', 'b0']
references = []

#Initializing dictionaries to be filled
stream_point = {}
stream = {}
groupstreamlines={}
groupLines = {}
groupPoints = {}
group_qb = {}
group_clusters = {}
groups_subjects = {}

subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361', 'S02363',
        'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02469', 'S02473',
        'S02485', 'S02491', 'S02490', 'S02506', 'S02523', 'S02524', 'S02535', 'S02654', 'S02666', 'S02670', 'S02686',
        'S02690', 'S02695', 'S02715', 'S02720', 'S02737', 'S02745', 'S02753', 'S02765', 'S02771', 'S02781', 'S02802',
        'S02804', 'S02813', 'S02812', 'S02817', 'S02840', 'S02842', 'S02871', 'S02877', 'S02898', 'S02926', 'S02938',
        'S02939', 'S02954', 'S02967', 'S02987', 'S03010', 'S03017', 'S03028', 'S03033', 'S03034', 'S03045', 'S03048',
        'S03069', 'S03225', 'S03265', 'S03293', 'S03308', 'S03321', 'S03343', 'S03350', 'S03378', 'S03391', 'S03394']

random.shuffle(subjects)
#removed_list = ['S02266']
removed_list = ['S02523']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)


_, _, index_to_struct, _ = atlas_converter(ROI_legends)
labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends,
                                                                             labeltype=labeltype, verbose=verbose)
for subject in subjects:
    trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned=False, verbose=verbose)
    if not exists:
        txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
        warnings.warn(txt)
        continue

    M_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_connectomes.xlsx")
    grouping_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_grouping.xlsx")

    if (os.path.exists(M_xlsxpath) or os.path.exists(grouping_xlsxpath)) and not overwrite:
        print(f'Found written file for subject {subject} at {M_xlsxpath} and {grouping_xlsxpath}')
        continue
    else:
        t1 = time()
        trkdata = load_trk(trkpath, 'same')
        if verbose:
            print(f"Time taken for loading the trk file {trkpath} set was {str((- t1 + time())/60)} minutes")
        t2 = time()
        header = trkdata.space_attributes


        streamlines_world = transform_streamlines(trkdata.streamlines, np.linalg.inv(labelaffine))

        if function_processes == 1:

            M, _, _, _, grouping = connectivity_matrix_custom(streamlines_world, np.eye(4), labelmask,
                                              inclusive=inclusive, symmetric=symmetric, return_mapping=True,
                                              mapping_as_streamlines=False, reference_weighting = None, volume_weighting=False)
        else:
            M, _, _, _, grouping = connectivity_matrix_func(streamlines_world, np.eye(4), labelmask, inclusive = inclusive,
                                                   symmetric=symmetric, return_mapping=True, mapping_as_streamlines=False, reference_weighting = None,
                                                   volume_weighting=False, verbose = False)

        M_grouping_excel_save(M,grouping,M_xlsxpath, grouping_xlsxpath, index_to_struct, verbose=False)

        del(trkdata)
        if verbose:
            print(f"Time taken for creating this connectome was set at {str((- t2 + time())/60)} minutes")
        if os.path.exists(grouping_xlsxpath) and verbose:
            print(f'Saved grouping for subject {subject} at {grouping_xlsxpath}')
        #grouping = extract_grouping(grouping_xlsxpath, index_to_struct, np.shape(M), verbose=verbose)
        else:
            raise Exception(f'saving of the excel at {grouping_xlsxpath} did not work')
