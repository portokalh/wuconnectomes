
import numpy as np
from dipy.segment.clustering import QuickBundles
from dipy.io.streamline import load_trk, save_trk
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric,mdf
from dipy.io.image import load_nifti
import warnings

from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines
import os, glob
import pickle
from nifti_handler import getlabeltypemask
from file_tools import mkcdir, check_files
from tract_handler import ratio_to_str, gettrkpath
from convert_atlas_mask import atlas_converter
import errno
import socket
from tract_save import save_trk_header
from excel_management import M_grouping_excel_save, extract_grouping
import sys
from argument_tools import parse_arguments_function
from connectome_handler import connectivity_matrix_func
from dipy.tracking.utils import length
from dipy.viz import window, actor
from time import sleep
from dipy.segment.clustering import ClusterCentroid
from dipy.tracking.streamline import Streamlines
from tract_visualize import show_bundles, setup_view
from dipy.tracking.utils import connectivity_matrix
from tract_save import unload_trk

def get_grouping(grouping_xlsx):
    print('not done yet')

def get_diff_ref(label_folder, subject, ref):
    diff_path = os.path.join(label_folder,f'{subject}_{ref}_to_MDT.nii.gz')
    if os.path.exists(diff_path):
        return diff_path
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), diff_path)

"""
'1 Cerebellum-Cortex_Right---Cerebellum-Cortex_Left 9 1 with weight of 3053.5005\n'
    '2 inferiortemporal_Left---Cerebellum-Cortex_Left 24 1 with weight of 463.1322\n'
    '3 inferiortemporal_Right---inferiorparietal_Right 58 57 with weight of 435.9886\n'
    '4 middletemporal_Right---inferiorparietal_Right 64 57 with weight of 434.9106\n'
    '5 fusiform_Left---Cerebellum-Cortex_Left 22 1 with weight of 402.0991\n'
"""
#target_tuple = (9, 1)
#target_tuple = (64, 57)
#target_tuple = (24, 1)
#target_tuple = (58, 57)
#target_tuple = (64, 57)
#target_tuple = (22, 1)
#target_tuple = (30, 50) #The connectomes to check up on and create groupings clusters for

#target_tuple = (39,32)

#set parameter
num_points1 = 50
distance1 = 1
#group cluster parameter
num_points2 = 50
distance2 = 2

ratio = 1
project = 'AD_Decode'
skip_subjects = True
write_streamlines = True
allow_preprun = False
verbose=True
picklesave=True
overwrite=False
inclusive = False
symmetric = True
write_stats = True
write_txt = True
constrain_groups = True

#target_tuples = [(9, 1), (24,1), (22, 1), (58, 57), (64, 57),(23,24),(24,30),(23,30)]
#target_tuples = [(9, 1), (77,43), (58,57), (24,1), (22,1)]
#target_tuples = [(58, 30), (58,45), (64,30), (58,24), (64,45)]
#target_tuples = [(64,57,(58,57),(64,58))]
#target_tuples = [(58,24), (58, 30), (64,30), (64,24), (58,48)]
#target_tuples = [(9,1), (57, 9), (61,23), (84,23), (80,9)]

#target_tuples.reverse()
#target_tuples = target_tuples[:3]
#target_tuples = [(9,1)]
#target_tuples = [(24,30),(23,30),(23,24)]


#genotype_noninclusive
target_tuples = [(9, 1), (24, 1), (58, 57), (64, 57), (22, 1)]
#target_tuples = [(24, 1)]
#genotype_noninclusive_volweighted_fa
#target_tuples = [(9, 1), (57, 9), (61, 23), (84, 23), (80, 9)]

#sex_noninclusive
#target_tuples = [(64, 57), (58, 57), (9, 1), (64, 58), (80,58)]
#sex_noninclusive_volweighted_fa
#target_tuples = [(58, 24), (58, 30), (64, 30), (64, 24), (58,48)]


labeltype = 'lrordered'
#reference_img refers to statistical values that we want to compare to the streamlines, say fa, rd, etc
references = ['fa', 'md', 'rd', 'ad', 'b0']
references = ['fa', 'md']
references = ['fa', 'md', 'ln', 'rd', 'ad']

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

computer_name = socket.gethostname()

samos = False
if 'samos' in computer_name:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in computer_name or 'hydra' in computer_name:
    #mainpath = '/Users/alex/jacques/'
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    ref_MDT_folder = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/reg_images/'
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
    ref_MDT_folder = '/mnt/munin6/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/reg_images/'
else:
    raise Exception('No other computer name yet')

# Setting identification parameters for ratio, labeling type, etc
ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio', '')

str_identifier = f'_stepsize_2{ratio_str}_wholebrain_pruned'
labeltype = 'lrordered'

function_processes = parse_arguments_function(sys.argv)
print(f'there are {function_processes} function processes')

if project=='AD_Decode':
    mainpath=os.path.join(mainpath,project,'Analysis')
else:
    mainpath = os.path.join(mainpath, project)
TRK_folder = os.path.join(mainpath, f'TRK_MPCA_MDT_fixed{folder_ratio_str}')

label_folder = os.path.join(mainpath, 'DWI')
if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'


trkpaths = glob.glob(os.path.join(TRK_folder, '*trk'))
pickle_folder = os.path.join(mainpath, f'Pickle_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
centroid_folder = os.path.join(mainpath, f'Centroids_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
stats_folder = os.path.join(mainpath, f'Statistics_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
excel_folder = os.path.join(mainpath, f'Excels_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
mkcdir([pickle_folder, centroid_folder, stats_folder, excel_folder])
if not os.path.exists(TRK_folder):
    raise Exception(f'cannot find TRK folder at {TRK_folder}')

#Initializing dictionaries to be filled
stream_point = {}
stream = {}
groupstreamlines={}
groupstreamlines_orig={}
groupLines = {}
groupPoints = {}
group_qb = {}
group_clusters = {}
groups_subjects = {}

if project == 'AD_Decode':
    groups_subjects['APOE3'] = ['S02402','S02720','S02812','S02373','S02231','S02410','S01912','S02451','S02485','S02473','S02506','S02524','S02535','S02686','S02695','S02753','S02765','S02804','S02817','S02842','S02871','S02926','S02938','S02939','S02967','S02320','S02110','S02289','S03017','S03010','S02987','S02227','S03033','S03034','S03069','S03308','S03321','S03350','S02266']
    groups_subjects['APOE4']= ['S02363','S02386','S02421','S02424','S02446','S02491','S02654','S02666','S02690','S02715','S02737','S02771','S02781','S02802','S02813','S02840','S02224','S02877','S02898','S02954','S02361','S02390','S02670','S03045','S03048','S03225','S03265','S03293','S03343','S03378','S03391']
    groups_subjects['APOEtestrun'] = ['S02386','S02363']

    groups_subjects['Male'] =['S01912', 'S02110', 'S02231', 'S02402', 'S02469', 'S02473', 'S02491', 'S02535', 'S02654', 'S02289', 'S02266', 'S02666', 'S02670', 'S02690', 'S02753', 'S02227', 'S02813', 'S02842', 'S02224', 'S02871', 'S02938', 'S02939', 'S02954', 'S02987', 'S03010', 'S02320', 'S03017', 'S03028', 'S03048', 'S03069', 'S03225', 'S03265', 'S03293', 'S03350', 'S03391']
    groups_subjects['Female'] = ['S02363', 'S02373', 'S02386', 'S02390', 'S02410', 'S02421', 'S02424', 'S02446', 'S02451', 'S02506', 'S02524', 'S02686', 'S02695', 'S02715', 'S02720', 'S02737', 'S02765', 'S02771', 'S02781', 'S02802', 'S02804', 'S02812', 'S02817', 'S02840', 'S02877', 'S02898', 'S02926', 'S02967', 'S03033', 'S03034', 'S03045', 'S02361', 'S03308', 'S03321', 'S03343', 'S03378']

    #groups to go through
    groups_all = ['APOE4','APOE3']
    groups= ['APOE3', 'APOE4']

    #groups = ['APOE3']
    #groups = ['Male','Female']
    #groups =['Female']
    #groups = ['APOEtestrun']

removed_list = ['S02654','S02523']

for group in groups:
    for remove in removed_list:
        if remove in groups_subjects[group]:
            groups_subjects[group].remove(remove)

if constrain_groups:
    group_sizes = []
    for group in groups_all:
        #group_sizes[group] = np.size(groups_subjects[group])
        group_sizes.append(np.size(groups_subjects[group]))
    group_min = np.min(group_sizes)
    for group in groups_all:
        groups_subjects[group] = groups_subjects[group][:group_min]
    print(group_sizes)

if project == 'APOE':
    raise Exception('not implemented')


feature1 = ResampleFeature(nb_points=num_points1)
metric1 = AveragePointwiseEuclideanMetric(feature=feature1)

feature2 = ResampleFeature(nb_points=num_points2)
metric2 = AveragePointwiseEuclideanMetric(feature=feature2)

overwrite=True

for target_tuple in target_tuples:

    for group in groups:
        groupstreamlines[group] = []
        groupstreamlines_orig[group] = []
        for ref in references:
            groupLines[group, ref] = []
            groupPoints[group, ref] = []

    _, _, index_to_struct, _ = atlas_converter(ROI_legends)
    print(f'Starting the run for {index_to_struct[target_tuple[0]]} to {index_to_struct[target_tuple[1]]}')

    for group in groups:
        print(f'Going through group {group}')

        group_str = group.replace(' ', '_')
        centroid_file_path = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_centroid.py')
        streamline_file_path = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_streamlines.trk')
        stats_path = os.path.join(stats_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_stats.xlsx')
        if write_stats:
            import xlsxwriter
            workbook = xlsxwriter.Workbook(stats_path)
            worksheet = workbook.add_worksheet()
            l=1
            for ref in references:
                worksheet.write(0,l, ref + ' mean')
                worksheet.write(0,l+1, ref + ' min')
                worksheet.write(0,l+2, ref + ' max')
                worksheet.write(0,l+3, ref + ' std')
                l=l+4
            #if verbose:
            #    print(f'Saved connectome at {output_path}')
            #streamline_file_path_orig = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_streamlines.trk')
        grouping_files = {}
        exists=True

        for ref in references:
            grouping_files[ref,'lines']=(os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_' + ref + '_lines.py'))
            grouping_files[ref, 'points'] = (os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_' + ref + '_points.py'))
            list_files, exists = check_files(grouping_files)
        if not os.path.exists(centroid_file_path) or not np.all(exists) or (not os.path.exists(streamline_file_path) and write_streamlines) or (not os.path.exists(stats_path) and write_stats) or overwrite:
            subjects = groups_subjects[group]
            subj = 1
            for subject in subjects:
                trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned=False, verbose=True)
                if not exists:
                    txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
                    warnings.warn(txt)
                    continue
                #streamlines, header, _ = unload_trk(trkpath)
                if np.shape(groupLines[group, ref])[0] != np.shape(groupstreamlines[group])[0]:
                    raise Exception('happened from there')
                trkdata = load_trk(trkpath, 'same')
                header = trkdata.space_attributes
                picklepath_connectome = os.path.join(pickle_folder, subject + str_identifier + '_connectomes.p')
                picklepath_grouping = os.path.join(pickle_folder, subject + str_identifier + '_grouping.p')
                M_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_connectomes.xlsx")
                grouping_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_grouping.xlsx")
                #if os.path.exists(picklepath_grouping) and not overwrite:
                #    with open(picklepath_grouping, 'rb') as f:
                #        grouping = pickle.load(f)
                if os.path.exists(picklepath_connectome):
                    with open(picklepath_connectome, 'rb') as f:
                        M = pickle.load(f)
                if os.path.exists(grouping_xlsxpath):
                    grouping = extract_grouping(grouping_xlsxpath, index_to_struct, None, verbose=verbose)
                else:
                    if allow_preprun:
                        labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT',
                                                                                                 ROI_legends,
                                                                                                 labeltype=labeltype,
                                                                                                 verbose=verbose)
                        streamlines_world = transform_streamlines(trkdata.streamlines, np.linalg.inv(labelaffine))

                        #M, grouping = connectivity_matrix_func(trkdata.streamlines, function_processes, labelmask,
                        #                                       symmetric=True, mapping_as_streamlines=False,
                        #                                       affine_streams=trkdata.space_attributes[0],
                        #                                       inclusive=inclusive)
                        M, grouping = connectivity_matrix_func(streamlines_world, np.eye(4), labelmask, inclusive=inclusive,
                                                 symmetric=symmetric, return_mapping=True, mapping_as_streamlines=False,
                                                 reference_weighting=None,
                                                 volume_weighting=False, verbose=False)
                        M_grouping_excel_save(M, grouping, M_xlsxpath, grouping_xlsxpath, index_to_struct,
                                              verbose=False)
                    else:
                        print(f'skipping subject {subject} for now as grouping file is not calculated. Best rerun it afterwards ^^')
                        continue

                target_streamlines_list = grouping[target_tuple[0], target_tuple[1]]
                if np.size(target_streamlines_list) == 0:
                    txt = f'Did not have any streamlines for {index_to_struct[target_tuple[0]]} to {index_to_struct[target_tuple[1]]} for subject {subject}'
                    warnings.warn(txt)
                    continue
                target_streamlines = trkdata.streamlines[np.array(target_streamlines_list)]
                target_streamlines_set = set_number_of_points(target_streamlines, nb_points=num_points2)
                #del(target_streamlines, trkdata)
                target_qb = QuickBundles(threshold=distance1, metric=metric1)

                if write_stats:
                    l = 1
                    worksheet.write(subj, 0, subject)
                for ref in references:

                    if ref != 'ln':
                        ref_img_path = get_diff_ref(ref_MDT_folder, subject, ref)
                        ref_data, ref_affine = load_nifti(ref_img_path)

                        from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
                        from collections import defaultdict, OrderedDict
                        from itertools import combinations, groupby

                        edges = np.ndarray(shape=(3, 0), dtype=int)
                        lin_T, offset = _mapping_to_voxel(trkdata.space_attributes[0])
                        stream_ref = []
                        stream_point_ref = []
                        for sl, _ in enumerate(target_streamlines_set):
                            # Convert streamline to voxel coordinates
                            entire = _to_voxel_coordinates(target_streamlines_set[sl], lin_T, offset)
                            i, j, k = entire.T
                            ref_values = ref_data[i, j, k]
                            stream_point_ref.append(ref_values)
                            stream_ref.append(np.mean(ref_values))
                    else:
                        stream_ref = list(length(target_streamlines))
                    """
                    from dipy.viz import window, actor
                    from tract_visualize import show_bundles, setup_view
                    import nibabel as nib

                    lut_cmap = actor.colormap_lookup_table(
                        scale_range=(0.05, 0.3))

                    scene = setup_view(nib.streamlines.ArraySequence(target_streamlines[33:34]), colors=lut_cmap,
                                       ref=ref_img_path, world_coords=True,
                                       objectvals=[None], colorbar=True, record=None, scene=None, interactive=True)
                    """

                    if write_stats:
                        worksheet.write(subj, l, np.mean(stream_ref))
                        worksheet.write(subj, l+1, np.min(stream_ref))
                        worksheet.write(subj, l+2, np.max(stream_ref))
                        worksheet.write(subj, l+3, np.std(stream_ref))
                        l=l+4
                    if not (group, ref) in groupLines.keys():
                        groupLines[group, ref]=(stream_ref)
                    else:
                        groupLines[group, ref].extend(stream_ref)
                        groupPoints[group, ref].extend(stream_point_ref)
                subj += 1
                groupstreamlines[group].extend(target_streamlines_set)


            if write_stats:
                worksheet.write(subj, 0, group)
                l=1
                for ref in references:
                    worksheet.write(subj, l, np.mean(groupLines[group, ref]))
                    worksheet.write(subj, l + 1, np.min(groupLines[group, ref]))
                    worksheet.write(subj, l + 2, np.max(groupLines[group, ref]))
                    worksheet.write(subj, l + 3, np.std(groupLines[group, ref]))
                    l=l+4
                workbook.close()

                #groupstreamlines_orig[group].extend(target_streamlines)


            group_qb[group] = QuickBundles(threshold=distance2, metric=metric2)
            group_clusters[group] = group_qb[group].cluster(groupstreamlines[group])
            if os.path.exists(centroid_file_path) and overwrite:
                os.remove(centroid_file_path)
            if not os.path.exists(centroid_file_path):
                if verbose:
                    print(f'Summarized the clusters for group {group} at {centroid_file_path}')
                pickle.dump(group_clusters[group], open(centroid_file_path, "wb"))


            if np.shape(groupLines[group, ref])[0] != np.shape(groupstreamlines[group])[0]:
                raise Exception('happened from there')

            if os.path.exists(streamline_file_path) and overwrite and write_streamlines:
                os.remove(streamline_file_path)
            if not os.path.exists(streamline_file_path) and write_streamlines:
                if verbose:
                    print(f'Summarized the streamlines for group {group} at {streamline_file_path}')
                pickle.dump(groupstreamlines[group], open(streamline_file_path, "wb"))
                save_trk_header(filepath= streamline_file_path, streamlines = groupstreamlines[group], header = header, affine=np.eye(4), verbose=verbose)

            """
            if not os.path.exists(streamline_file_path_orig) and write_streamlines:
                if verbose:
                    print(f'Summarized the streamlines for group {group} at {streamline_file_path}')
                pickle.dump(groupstreamlines_orig[group], open(streamline_file_path, "wb"))
                save_trk_header(filepath= streamline_file_path, streamlines = groupstreamlines_orig[group], header = header, affine=np.eye(4), verbose=verbose)
            """

            for ref in references:
                if overwrite:
                    if os.path.exists(grouping_files[ref,'lines']):
                        os.remove(grouping_files[ref,'lines'])
                    if os.path.exists(grouping_files[ref,'points']):
                        os.remove(grouping_files[ref,'points'])
                if not os.path.exists(grouping_files[ref,'lines']):
                    if verbose:
                        print(f"Summarized the clusters for group {group} and statistics {ref} at {grouping_files[ref,'lines']}")
                    pickle.dump(groupLines[group, ref], open(grouping_files[ref,'lines'], "wb"))
                if not os.path.exists(grouping_files[ref, 'points']):
                    if verbose:
                        print(f"Summarized the clusters for group {group} and statistics {ref} at {grouping_files[ref,'lines']}")
                    pickle.dump(groupPoints[group, ref], open(grouping_files[ref,'points'], "wb"))
            pickle.dump(groupLines[group, 'ln'], open(grouping_files['ln', 'lines'], "wb"))

        else:
            print(f'Centroid file was found at {centroid_file_path}, reference files for {references}')
            with open(centroid_file_path, 'rb') as f:
                group_clusters[group] = pickle.load(f)
            for ref in references:
                ref_path_lines = grouping_files[ref, 'lines']
                with open(ref_path_lines, 'rb') as f:
                    groupLines[group,ref] = pickle.load(f)
                #ref_path_points = grouping_files[ref, 'points']
                #groupPoints[group, ref] = grouping_files[ref, 'points']

ref_mean = {}
for reference in references:
    for group in groups:
        ref_mean[reference,group] = np.mean(groupLines[group,ref])

for group in groups:
    cluster = group_clusters[group]
    group_str = group.replace(' ', '_')
    idx_path = os.path.join(centroid_folder,
                            group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                            index_to_struct[target_tuple[1]] + '_idx.py')
    if os.path.exists(idx_path):
        continue
    else:
        top_idx_list = sorted(range(len(cluster.clusters_sizes())), key=lambda i: cluster.clusters_sizes()[i],
                              reverse=True)
        if verbose:
            print(f'Listed the biggest clusters for group {group} at {idx_path}')
        pickle.dump(top_idx_list, open(idx_path, "wb"))

"""
toview=False
if toview:
    #group_toview = 'Initial AMD'
    viz_top_bundle = True
    ref = None
    ref = '/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz'

    num_of_bundles = 5
    color_list = [(np.random.randint(0, 255), np.random.randint(0, 255), np.random.randint(0, 255))
                  for n in range(num_of_bundles)]
    color_list_dis_all = [window.colors.green, window.colors.yellow,
                          window.colors.red, window.colors.brown,
                          window.colors.orange, window.colors.blue]
    color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]

    if num_of_bundles <= 6:
        colors = color_list_dis
    else:
        colors = color_list

    cluster = group_clusters[group_toview]
    group_str = group_toview.replace(' ', '_')
    idx_path = os.path.join(centroid_folder,group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_idx.py')
    if os.path.exists(idx_path):
        with open(idx_path, 'rb') as f:
            top_idx_list = pickle.load(f)
    else:
        top_idx_list = sorted(range(len(cluster.clusters_sizes())), key=lambda i: cluster.clusters_sizes()[i],
                     reverse=True)
        pickle.dump(top_idx_list, open(idx_path, "wb"))
    top_idx = top_idx_list[:num_of_bundles]

    bundle_list = [cluster.clusters[idx] for idx in top_idx]

    setup_view(bundle_list, colors = colors,ref=ref, world_coords=True)

groups_toview = ['Paired 2-YR Control','Paired 2-YR AMD' ]
toview_multi = False
num_of_bundles = 10

if toview_multi:
    ref = '/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz'

    num_of_groups = np.size(groups_toview)

    color_list = [(np.random.randint(0, 255), np.random.randint(0, 255), np.random.randint(0, 255))
                  for n in range(num_of_groups)]
    color_list_dis_all = [window.colors.green, window.colors.yellow,
                          window.colors.red, window.colors.brown,
                          window.colors.orange, window.colors.blue]
    color_list_dis = [color_list_dis_all[i] for i in range(num_of_groups)]
    if num_of_groups <= 6:
        colors = color_list_dis
    else:
        colors = color_list
    num_of_bundles = 10
    bundle_superlist = []
    for group in groups_toview:
        cluster = group_clusters[group]
        group_str = group.replace(' ', '_')
        idx_path = os.path.join(centroid_folder,
                                group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                index_to_struct[target_tuple[1]] + '_idx.py')
        if os.path.exists(idx_path):
            with open(idx_path, 'rb') as f:
                top_idx_list = pickle.load(f)
        else:
            top_idx_list = sorted(range(len(cluster.clusters_sizes())), key=lambda i: cluster.clusters_sizes()[i],
                                  reverse=True)
            pickle.dump(top_idx_list, open(idx_path, "wb"))
        top_idx = top_idx_list[:num_of_bundles]

        bundle_list = [cluster.clusters[idx] for idx in top_idx]
        bundle_superlist.append(bundle_list)

    setup_view(bundle_superlist, colors = colors,ref=ref, world_coords=True)


if viz_top_bundle:
    np.random.seed(123)
    num_of_bundles = 5

    cluster = group_clusters[group_toview]
    name = f'Group_{group_toview}' + str(num_of_bundles)

    top_idx = sorted(range(len(cluster.clusters_sizes())), key=lambda i: cluster.clusters_sizes()[i],
                     reverse=True)[:num_of_bundles]

    bundle_list = [cluster.clusters[idx] for idx in top_idx]
    color_list = [(np.random.randint(0, 255), np.random.randint(0, 255), np.random.randint(0, 255))
                  for n in range(num_of_bundles)]
    color_list_dis_all = [window.colors.green, window.colors.yellow,
                          window.colors.red, window.colors.brown,
                          window.colors.orange, window.colors.blue]
    color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]

    if num_of_bundles <= 6:
        colors = color_list_dis
    else:
        colors = color_list

    show_bundles(bundle_list, colors, ref=ref)


#color by lines, select a bundle?
np.random.seed(123)
bundle_id = 40
ref_toview = ['fa']

if viz_top_bundle:

    clusters = group_clusters[group_toview]
    groupLines = groupLines[group_toview, ref_toview]
    name = f'Group_Gen3-Bundle {str(bundle_id)}'

    top_idx = sorted(range(len(clusters.clusters_sizes())), key=lambda i: clusters.clusters_sizes()[i],
                     reverse=True)[:num_of_bundles]

    k = clusters.clusters[bundle_id]
    bundle_ref = []
    for idx in k.indices:
        bundle_ref.append(groupLines[idx])

    #         cmap = actor.colormap_lookup_table(
    #         scale_range=(np.min(bundle_ref), np.max(bundle_ref)))
    cmap = actor.colormap_lookup_table(
        scale_range=(0.1, 0.5))

    # color by line-average fa
    renderer = window.Renderer()
    renderer.clear()
    renderer = window.Renderer()
    stream_actor3 = actor.line(clusters.clusters[bundle_id], np.array(bundle_ref), lookup_colormap=cmap)
    renderer.add(stream_actor3)
    bar = actor.scalar_bar(cmap)
    renderer.add(bar)
    # Uncomment the line below to show to display the window
    window.show(renderer, size=(600, 600), reset_camera=False)

# viz top bundle
np.random.seed(123)
num_of_bundles = 5

if viz_top_bundle:

    clusters = group_clusters[group_toview]
    name = f'Group_{group_toview}-Bundle top ' + str(num_of_bundles)

    top_idx = sorted(range(len(group_clusters.clusters_sizes())), key=lambda i: group_clusters.clusters_sizes()[i],
                     reverse=True)[:num_of_bundles]

    bundle_list = [group_clusters.clusters[idx] for idx in top_idx]
    color_list = [(np.random.randint(0, 255), np.random.randint(0, 255), np.random.randint(0, 255))
                  for n in range(num_of_bundles)]
    color_list_dis_all = [window.colors.green, window.colors.yellow,
                          window.colors.red, window.colors.brown,
                          window.colors.orange, window.colors.blue]
    color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]

    if num_of_bundles <= 6:
        colors = color_list_dis
    else:
        colors = color_list

    show_bundles(bundle_list, colors, fa=1)


group_qb = {}
group_clusters = {}
for group in groups:
    group_str = group.replace(' ', '_')
    centroid_file_path = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '.py')
    if not os.path.exists(centroid_file_path):
        group_qb[group] = QuickBundles(threshold=distance2, metric=metric2)
        group_clusters[group] = group_qb[group].cluster(groupstreamlines[group])
        pickle.dump(grouping, open(centroid_file_path, "wb"))
    else:
        if os.path.exists(picklepath_grouping):
            with open(picklepath_grouping, 'rb') as f:
                grouping = pickle.load(f)
"""

    #save_trk(group_qb[group].cluster(groupstreamlines[group]), centroid_file_path)
    #save_trk_heavy_duty(centroid_file_path, streamlines=group_clusters[group], affine=np.eye(4), header=header)


#print("Young Group Nb. clusters:", len(group3_clusters))


# registration
# srr = StreamlineLinearRegistration()
##srm = srr.optimize(static=target_clusters_control.centroids, moving=target_clusters.centroids)
# target_str_aligned = srm.transform(target_streamlines)
# native_target_strea
