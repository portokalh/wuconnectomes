
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

project = 'AMD'

huma_projects = ''
hostname = socket.gethostname()


samos = False
if 'samos' in hostname:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in hostname:
    mainpath = '/Users/alex/jacques/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    print(f'no option for {hostname}')

mainpath = os.path.join(mainpath, project)
TRK_folder = os.path.join(mainpath, 'TRK_MDT')
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

if project == 'AMD':
    groups_subjects['testing'] = ['H22825']
    groups_subjects['Initial AMD'] = ['H27640', 'H27778', 'H29020', 'H26637', 'H27680', 'H26765', 'H27017', 'H26880', 'H28308', 'H28433', 'H28338', 'H26660', 'H28809', 'H27610', 'H26745', 'H27111', 'H26974', 'H27391', 'H28748', 'H29025', 'H29013', 'H27381', 'H26958', 'H28662', 'H26578', 'H28698', 'H27495', 'H28861', 'H28115', 'H28437', 'H26850', 'H28532', 'H28377', 'H28463', 'H26890', 'H28373', 'H28857', 'H27164', 'H27982']
    groups_subjects['Paired 2-YR AMD'] = ['H22825', 'H21850', 'H29225', 'H29304', 'H29060', 'H23210', 'H21836', 'H29618', 'H22644', 'H22574', 'H22369', 'H29627', 'H29056', 'H22536', 'H23143', 'H22320', 'H22898', 'H22864', 'H29264', 'H22683']
    groups_subjects['Initial Control'] = ['H26949', 'H27852', 'H28029', 'H26966', 'H27126', 'H28068', 'H29161', 'H28955', 'H26862', 'H28262', 'H28856', 'H27842', 'H27246', 'H27869', 'H27999', 'H29127', 'H28325', 'H26841', 'H29044', 'H27719', 'H27100', 'H29254', 'H27682', 'H29002', 'H29089', 'H29242', 'H27488', 'H27841', 'H28820', 'H27163', 'H28869', 'H28208', 'H27686']
    groups_subjects['Paired 2-YR Control'] = ['H29403', 'H22102', 'H29502', 'H22276', 'H29878', 'H29410', 'H22331', 'H22368', 'H21729', 'H29556', 'H21956', 'H22140', 'H23309', 'H22101', 'H23157', 'H21593', 'H21990', 'H22228', 'H23028', 'H21915']
    groups_subjects['Paired Initial Control'] = ['H27852', 'H28029', 'H26966', 'H27126', 'H29161', 'H28955', 'H26862', 'H27842', 'H27999', 'H28325', 'H26841', 'H27719', 'H27100', 'H27682', 'H29002', 'H27488', 'H27841', 'H28820', 'H28208', 'H27686']
    groups_subjects['Paired Initial AMD'] = ['H29020', 'H26637', 'H27111', 'H26765', 'H28308', 'H28433', 'H26660', 'H28182', 'H27111', 'H27391', 'H28748', 'H28662', 'H26578', 'H28698', 'H27495', 'H28861', 'H28115', 'H28377', 'H26890', 'H28373', 'H27164']

    #groups to go through
    groups = ['Initial AMD','Initial Control','Paired 2-YR Control','Paired Initial Control','Paired Initial AMD',
'Paired 2-YR AMD']

#groups = ['Paired 2-YR AMD']
#groups = ['Paired 2-YR Control']
#groups=[groups[0]]
group_toview = groups[0]

if project == 'AD_Decode':
    raise Exception('not implemented')

if project == 'APOE':
    raise Exception('not implemented')

if project == 'AD_Decode':
    raise Exception('not implemented')

for group in groups:
    groupstreamlines[group]=[]
    for ref in references:
        groupLines[group, ref]=[]
        groupPoints[group, ref]=[]


#Setting identification parameters for ratio, labeling type, etc
ratio = 1
ratio_str = ratio_to_str(ratio)
str_identifier = '_MDT'+ratio_str
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

for group in groups:
    group_str = group.replace(' ', '_')
    _, _, index_to_struct, _ = atlas_converter(ROI_legends)
    subjects = groups_subjects[group]
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
