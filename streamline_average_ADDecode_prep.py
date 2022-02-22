import numpy as np
from dipy.segment.clustering import QuickBundles
from dipy.io.streamline import load_trk, save_trk
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric,mdf
from dipy.io.image import load_nifti
import warnings
from dipy.viz import window, actor
from time import sleep
from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines
import os, glob
from tract_save import unload_trk
import pickle
from dipy.tracking.utils import connectivity_matrix
from nifti_handler import getlabeltypemask
from file_tools import mkcdir
from tract_handler import ratio_to_str, gettrkpath
from convert_atlas_mask import convert_labelmask, atlas_converter
import errno
import socket
from dipy.segment.clustering import ClusterCentroid
from dipy.tracking.streamline import Streamlines
from tract_visualize import show_bundles, setup_view
from tract_save import save_trk_header
from excel_management import M_grouping_excel_save, extract_grouping
import sys
from argument_tools import parse_arguments_function
from connectome_handler import connectivity_matrix_custom, connectivity_matrix_func
import random
from nibabel.streamlines.array_sequence import ArraySequence
from time import time

def get_grouping(grouping_xlsx):
    print('not done yet')

def get_diff_ref(label_folder, subject, ref):
    diff_path = os.path.join(label_folder,f'{subject}_{ref}_to_MDT.nii.gz')
    if os.path.exists(diff_path):
        return diff_path
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), diff_path)

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
ratio = 100
ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio','')

inclusive = False
symmetric = True
if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'


str_identifier = f'_stepsize_2{ratio_str}_wholebrain_pruned'
labeltype = 'lrordered'
verbose=True
picklesave=True

function_processes = parse_arguments_function(sys.argv)
print(f'there are {function_processes} function processes')
overwrite=False

if project=='AD_Decode':
    mainpath=os.path.join(mainpath,project,'Analysis')
else:
    mainpath = os.path.join(mainpath, project)
TRK_folder = os.path.join(mainpath, 'TRK_MPCA_MDT_fixed'+folder_ratio_str)
label_folder = os.path.join(mainpath, 'DWI')
trkpaths = glob.glob(os.path.join(TRK_folder, '*trk'))
#pickle_folder = os.path.join(mainpath, 'Pickle_MDT'+folder_ratio_str)
#centroid_folder = os.path.join(mainpath, 'Centroids_MDT'+folder_ratio_str)
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
#subjects = ['S02363']
#subjects = ['S02686']
#subjects = ['S02987']
#subjects = ['S02386']
#subjects = ['S02813']
random.shuffle(subjects)
#removed_list = ['S02266']
removed_list = ['S02523']
for remove in removed_list:
    if remove in subjects:
        subjects.remove(remove)

overwrite = True

_, _, index_to_struct, _ = atlas_converter(ROI_legends)
labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends,
                                                                             labeltype=labeltype, verbose=verbose)
for subject in subjects:
    trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned=False, verbose=verbose)
    if not exists:
        txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
        warnings.warn(txt)
        continue

    #picklepath_connectome = os.path.join(pickle_folder, subject + str_identifier + '_connectome.p')
    #picklepath_grouping = os.path.join(pickle_folder, subject + str_identifier + '_grouping.p')
    M_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_connectome.xlsx")
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

        #endlabels_list_all = [7559, 8618, 10070, 10872, 13321, 14375, 14390, 14752, 15397, 15795, 15829, 15876, 16160,
        #                      16919, 18050, 18446, 18843, 19263, 19633, 19695, 20046, 20065, 20069, 20076, 20493, 20497,
        #                      20501, 20523, 20927, 21359, 21373, 21683, 22167, 22570, 22586, 22914, 23390, 23838, 23842,
        #                      23851, 23877, 24172, 24945, 24950, 25344, 25383, 25416, 25746, 26064, 26124, 26158, 26820,
        #                      27854, 30558, 30918, 32315, 32366, 33667, 33711, 35699, 35745]

        #streamlines_world = streamlines_world[endlabels_list_all]

        #endpoingts[7559] = [[167 100 135] [ 95  82 135]]
        #endpoints_realones[7559] = [167,100,23] [95 82 26]
        """
        ref_MDT_path = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_dwi.nii.gz'
        scene = None

        lut_cmap = actor.colormap_lookup_table(
            scale_range=(0.01, 0.55))
        scene = setup_view(trkdata.streamlines[:2000], colors=lut_cmap, ref=ref_MDT_path, world_coords=True,
                           objectvals=[None],
                           colorbar=True, record=None, scene=scene)
        scene = setup_view(streamlines_world[:2000], colors=lut_cmap, ref=ref_MDT_path, world_coords=False,
                           objectvals=[None],
                           colorbar=True, record=None, scene=scene)
        """
        if function_processes == 1:
            #M, grouping = connectivity_matrix(trkdata.streamlines, trkdata.space_attributes[0], labelmask, inclusive=inclusive, symmetric=True, return_mapping=True, mapping_as_streamlines=False)
            M, _, _, _, grouping = connectivity_matrix_custom(streamlines_world, np.eye(4), labelmask,
                                              inclusive=inclusive, symmetric=symmetric, return_mapping=True,
                                              mapping_as_streamlines=False, reference_weighting = None, volume_weighting=False)
        else:
            #M, grouping = connectivity_matrix_func(trkdata.streamlines, function_processes, labelmask, symmetric = True, mapping_as_streamlines = False, affine_streams = trkdata.space_attributes[0], inclusive= inclusive, verbose=False)
            M, _, _, _, grouping = connectivity_matrix_func(streamlines_world, np.eye(4), labelmask, inclusive = inclusive,
                                                   symmetric=symmetric, return_mapping=True, mapping_as_streamlines=False, reference_weighting = None,
                                                   volume_weighting=False, verbose = False)

        """
        lut_cmap = actor.colormap_lookup_table(
            scale_range=(0.01, 0.55))
        scene = None
        ref_MDT_path = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_dwi.nii.gz'
        i = 1
        j = 9
        print(i, j)
        target_streamlines_list = grouping[i, j]
        #print(index_to_struct[i] + '_to_' + index_to_struct[j])

        if len(target_streamlines_list) > 0:
            target_streamlines = trkdata.streamlines[np.array(target_streamlines_list)]
            print(f'viewing {index_to_struct[i]} to {index_to_struct[j]}')
            scene = setup_view(target_streamlines, colors=lut_cmap, ref=ref_MDT_path, world_coords=True, objectvals=[None],
                      colorbar=False, record=None, scene=scene, interactive=True)
        else:
            print('hi')
        """
        M_grouping_excel_save(M,grouping,M_xlsxpath, grouping_xlsxpath, index_to_struct, verbose=False)

        del(trkdata)
        if verbose:
            print(f"Time taken for creating this connectome was set at {str((- t2 + time())/60)} minutes")
        if os.path.exists(grouping_xlsxpath) and verbose:
            print(f'Saved grouping for subject {subject} at {grouping_xlsxpath}')
        #grouping = extract_grouping(grouping_xlsxpath, index_to_struct, np.shape(M), verbose=verbose)
        else:
            raise Exception(f'saving of the excel at {grouping_xlsxpath} did not work')


        """
        from connectome_handler import _mapping_to_voxel, _to_voxel_coordinates_warning, retweak_points

        endpoints = np.array([sl[0::len(sl) - 1] for sl in streamlines_world[grouping[1,9]]])
        lin_T, offset = _mapping_to_voxel(np.eye(4))
        endpoints = _to_voxel_coordinates_warning(endpoints, lin_T, offset)
        endpoints = retweak_points(endpoints, np.shape(labelmask))
        i, j, k = endpoints.T
        endlabels = labelmask[i, j, k]
        print(endlabels)
        endlabels_list = []
        for i in np.arange(np.shape(endlabels)[1]):
            if (endlabels[0,i] == 1 and endlabels[1,i] == 9) or (endlabels[0,i] == 9 and endlabels[1,i] == 1):
                endlabels_list.append(i)
        'S02813'
        endlabels_list = [81, 113, 162, 198, 313, 373, 374, 391, 433, 460, 465, 471, 483, 532, 601, 629, 661, 705, 737, 758, 771, 775, 777, 779, 815, 817, 818, 825, 856, 900, 906, 918, 964, 996, 1002, 1023, 1069, 1108, 1109, 1110, 1115, 1125, 1176, 1178, 1206, 1213, 1219, 1237, 1258, 1271, 1277, 1313, 1389, 1532, 1547, 1610, 1615, 1663, 1673, 1734, 1736]
        endlabels_list_all = []
        for endlabels_sl in endlabels_list:
            np.append(grouping[1,9][endlabels_sl])
        endlabels_list_all = [7559, 8618, 10070, 10872, 13321, 14375, 14390, 14752, 15397, 15795, 15829, 15876, 16160, 16919, 18050, 18446, 18843, 19263, 19633, 19695, 20046, 20065, 20069, 20076, 20493, 20497, 20501, 20523, 20927, 21359, 21373, 21683, 22167, 22570, 22586, 22914, 23390, 23838, 23842, 23851, 23877, 24172, 24945, 24950, 25344, 25383, 25416, 25746, 26064, 26124, 26158, 26820, 27854, 30558, 30918, 32315, 32366, 33667, 33711, 35699, 35745]
            
        
        scene = setup_view(streamlines[11], colors=lut_cmap, ref=label_volume, world_coords=False,
                   objectvals=[None],
                   colorbar=True, record=None, scene=scene)
                   
                   
        from connectome_handler import _mapping_to_voxel, _to_voxel_coordinates_warning, retweak_points
from collections import defaultdict, OrderedDict

for streamline in streamlines_world[target_streamlines_list]:
    entire = _to_voxel_coordinates_warning(streamline, lin_T, offset)
    entire = retweak_points(entire, np.shape(labelmask))
    i, j, k = entire.T
    entirelabels = list(OrderedDict.fromkeys(labelmask[i, j, k]))
    print(entirelabels)
    print('\n\n')
        """

#old testing grounds, good visualization tools
"""
#ref_MDT_path = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_dwi.nii.gz'

from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
scene = None

streamlines = trkdata.streamlines
lin_T, offset = _mapping_to_voxel(header[0])
streamlines_world_connectome = []
for i in np.arange(10000):
    entire = _to_voxel_coordinates(streamlines[i], lin_T, offset)
    streamlines_world_connectome.append(entire)
streamlines_world_connectome = ArraySequence(streamlines_world_connectome)

target_isocenter = np.array([[-1., 0., 0., 128.], [0., -1., 0., 128.], [0., 0., 1., -67.], [0., 0., 0., 1.]])
target_isocenter[2, 3] = 67
MDT_affine = load_nifti(ref_MDT_path)[1]
streamlines_world_affine = transform_streamlines(streamlines, MDT_affine)

scene = setup_view(streamlines_world_affine[:2000], colors=lut_cmap, ref=ref_MDT_path, world_coords=False,
                   objectvals=[None],
                   colorbar=True, record=None, scene=scene)
scene = setup_view(streamlines_world_connectome[:2000], colors=lut_cmap, ref=ref_MDT_path, world_coords=False,
                   objectvals=[None],
                   colorbar=True, record=None, scene=scene)


for i in np.arange(1,85):
    for j in np.arange(1,85):
        i = 1
        j = 9
        print(i, j)
        target_streamlines_list = grouping[i, j]
        #print(index_to_struct[i] + '_to_' + index_to_struct[j])
        if len(target_streamlines_list) > 0:
            target_streamlines = trkdata.streamlines[np.array(target_streamlines_list)]
            print(f'viewing {index_to_struct[i]} to {index_to_struct[j]}')
            scene = setup_view(target_streamlines, colors=lut_cmap, ref=ref_MDT_path, world_coords=True, objectvals=[None],
                      colorbar=False, record=None, scene=scene, interactive=True)
        else:
            print(f'skipping {index_to_struct[i]} to {index_to_struct[j]}')

        #setup_view(target_streamlines, colors=lut_cmap, ref=ref_MDT_path, world_coords=True, objectvals=[None],
         #          colorbar=False, record=None, scene=None, interactive=True)
#setup_view(target_streamlines, colors=lut_cmap, ref=labeloutpath, world_coords=True, objectvals=[None],
#           colorbar=False, record=None, scene=None, interactive=True)
"""