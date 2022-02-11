
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
from file_tools import mkcdir, check_files
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
from tract_manager import connectivity_matrix_func

def get_grouping(grouping_xlsx):
    print('not done yet')

def get_diff_ref(label_folder, subject, ref):
    diff_path = os.path.join(label_folder,f'{subject}_{ref}_to_MDT.nii.gz')
    if os.path.exists(diff_path):
        return diff_path
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), diff_path)

#set parameter
num_points1 = 50
distance1 = 1
feature1 = ResampleFeature(nb_points=num_points1)
metric1 = AveragePointwiseEuclideanMetric(feature=feature1)

#group cluster parameter
num_points2 = 50
distance2 = 2
feature2 = ResampleFeature(nb_points=num_points2)
metric2 = AveragePointwiseEuclideanMetric(feature=feature2)

project = 'AD_Decode'

huma_projects = ''
hostname = socket.gethostname()


samos = False
if 'samos' in hostname:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in hostname:
    #mainpath = '/Users/alex/jacques/'
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    raise Exception('No other computer name yet')

if project=='AD_Decode':
    mainpath=os.path.join(mainpath,project,'Analysis')
else:
    mainpath = os.path.join(mainpath, project)
TRK_folder = os.path.join(mainpath, 'TRK_MPCA_MDT_fixed')
label_folder = os.path.join(mainpath, 'DWI')
trkpaths = glob.glob(os.path.join(TRK_folder, '*trk'))
figures_folder = os.path.join(mainpath, 'Figures_MDT')
pickle_folder = os.path.join(mainpath, 'Pickle_MDT')
centroid_folder = os.path.join(mainpath, 'Centroids_MDT')
excel_folder = os.path.join(mainpath, 'Excels_MDT')
mkcdir([figures_folder, pickle_folder, centroid_folder, excel_folder])
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

"""
if project == 'AD_Decode':
    groups_subjects['APOE3'] = ['S02402','S02266','S02720','S02812','S02373','S02231','S02410','S01912','S02451','S02485','S02473','S02506','S02524','S02535','S02686','S02695','S02753','S02765','S02804','S02817','S02842','S02871','S02926','S02938','S02939','S02967','S02320','S02110','S02289','S03017','S03010','S02987','S02227','S03033','S03034','S03069']
    groups_subjects['APOE4']= ['S02363','S02386','S02421','S02424','S02446','S02491','S02654','S02666','S02690','S02715','S02737','S02771','S02781','S02802','S02813','S02840','S02224','S02877','S02898','S02954','S02361','S02390','S02670','S03045','S03048']
    #groups to go through
    groups = ['APOE4','APOE3']

group_toview = groups[0]

for group in groups:
    groupstreamlines[group]=[]
    for ref in references:
        groupLines[group, ref]=[]
        groupPoints[group, ref]=[]
"""

subjects = ['S01912', 'S02110', 'S02224', 'S02227', 'S02230', 'S02231', 'S02266', 'S02289', 'S02320', 'S02361', 'S02363',
        'S02373', 'S02386', 'S02390', 'S02402', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02469', 'S02473',
        'S02485', 'S02491', 'S02490', 'S02506', 'S02523', 'S02524', 'S02535', 'S02654', 'S02666', 'S02670', 'S02686',
        'S02690', 'S02695', 'S02715', 'S02720', 'S02737', 'S02745', 'S02753', 'S02765', 'S02771', 'S02781', 'S02802',
        'S02804', 'S02813', 'S02812', 'S02817', 'S02840', 'S02842', 'S02871', 'S02877', 'S02898', 'S02926', 'S02938',
        'S02939', 'S02954', 'S02967', 'S02987', 'S03010', 'S03017', 'S03028', 'S03033', 'S03034', 'S03045', 'S03048',
        'S03069', 'S03225', 'S03265', 'S03293', 'S03308', 'S03321', 'S03343', 'S03350', 'S03378', 'S03391', 'S03394']

#Setting identification parameters for ratio, labeling type, etc
ratio = 1
ratio_str = ratio_to_str(ratio)
str_identifier = '_MDT'+ratio_str
str_identifier = '_stepsize_2_all_wholebrain_pruned'
labeltype = 'lrordered'
verbose=True
picklesave=True

function_processes = parse_arguments_function(sys.argv)
print(f'there are {function_processes} function processes')
overwrite=False

"""
labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends,
                                                                         labeltype=labeltype, verbose=verbose)
trkpath = '/Users/alex/jacques/AMD/TRK_MDT/H29403_MDT_ratio_100.trk'
trkdata = load_trk(trkpath, 'same')
streamlines = trkdata.streamlines

setup_view(streamlines[0:6], ref=labeloutpath, world_coords=True)

labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends,
                                                                         labeltype=labeltype, verbose=verbose)
trkpath = '/Users/alex/jacques/AMD/TRK_MDT/H29403_MDT_ratio_100.trk'
trk_testpath = '/Users/alex/jacques/AMD/TRK_MDT/H29403_MDT_ratio_firststreamlines.trk'
trkdata = load_trk(trkpath, 'same')
header = trkdata.space_attributes
save_trk_header(filepath=trk_testpath, streamlines=trkdata.streamlines[0:6], header=header,
                    affine=np.eye(4), verbose=verbose)
setup_view(trkdata.streamlines[0:6], ref=labeloutpath, world_coords=True)
"""


_, _, index_to_struct, _ = atlas_converter(ROI_legends)
labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends,
                                                                             labeltype=labeltype, verbose=verbose)
for subject in subjects:
    trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned=False, verbose=verbose)
    if not exists:
        txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
        warnings.warn(txt)
        continue
    #streamlines, header, _ = unload_trk(trkpath)
    picklepath_connectome = os.path.join(pickle_folder, subject + str_identifier + '_connectome.p')
    picklepath_grouping = os.path.join(pickle_folder, subject + str_identifier + '_grouping.p')
    M_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_connectome.xlsx")
    grouping_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_grouping.xlsx")
    #if os.path.exists(picklepath_grouping) and not overwrite:
    #    with open(picklepath_grouping, 'rb') as f:
        #        grouping = pickle.load(f)
    if os.path.exists(picklepath_connectome) or os.path.exists(grouping_xlsxpath):
        print(f'Found written file for subject {subject}')
        continue
    else:
        trkdata = load_trk(trkpath, 'same')
        header = trkdata.space_attributes
        if function_processes == 1:
            M, grouping = connectivity_matrix(trkdata.streamlines, trkdata.space_attributes[0], labelmask, inclusive=True, symmetric=True,
                                    return_mapping=True,
                                    mapping_as_streamlines=False, verbose=verbose)
        else:
            M, grouping_temp = connectivity_matrix_func(trkdata.streamlines, function_processes, labelmask, symmetric = True, mapping_as_streamlines = False, affine_streams = trkdata.space_attributes[0], inclusive= True, verbose=verbose)
        M_grouping_excel_save(M,grouping_temp,M_xlsxpath, grouping_xlsxpath, index_to_struct, verbose=False)
    del(trkdata)
    if os.path.exists(grouping_xlsxpath):
        grouping = extract_grouping(grouping_xlsxpath, index_to_struct, np.shape(M), verbose=verbose)
    else:
        raise Exception(f'saving of the excel at {grouping_xlsxpath} did not work')

