
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
    mainpath = '/Users/alex/jacques/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    print(f'no option for {hostname}')

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

verbose = True

#Initializing dictionaries to be filled
stream_point = {}
stream = {}
groupstreamlines={}
groupLines = {}
groupPoints = {}
group_qb = {}
group_clusters = {}
groups_subjects = {}

if project == 'AD_Decode':
    groups_subjects['APO3'] = ['S02402','S02266','S02720','S02812','S02373','S02231','S02410','S01912','S02451','S02485','S02473','S02506','S02524','S02535','S02686','S02695','S02753','S02765','S02804','S02817','S02842','S02871','S02926','S02938','S02939','S02967','S02320','S02110','S02289','S03017','S03010','S02987','S02227','S03033','S03034','S03069']
    groups_subjects['APOE4']= ['S02363','S02386','S02421','S02424','S02446','S02491','S02654','S02666','S02690','S02715','S02737','S02771','S02781','S02802','S02813','S02840','S02224','S02877','S02898','S02954','S02361','S02390','S02670','S03045','S03048']
    #groups to go through
    groups = ['APO3','APOE4']

#groups = ['Paired 2-YR AMD']
#groups = ['Paired 2-YR Control']
#groups=[groups[0]]
group_toview = groups[0]


if project == 'APOE':
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

"""
'1 Cerebellum-Cortex_Right---Cerebellum-Cortex_Left 9 1 with weight of 3053.5005\n'
    '2 inferiortemporal_Left---Cerebellum-Cortex_Left 24 1 with weight of 463.1322\n'
    '3 inferiortemporal_Right---inferiorparietal_Right 58 57 with weight of 435.9886\n'
    '4 middletemporal_Right---inferiorparietal_Right 64 57 with weight of 434.9106\n'
    '5 fusiform_Left---Cerebellum-Cortex_Left 22 1 with weight of 402.0991\n'
"""
target_tuple = (9, 1)
#target_tuple = (24, 23)
#target_tuple = (30, 23)
#target_tuple = (24, 1)
#target_tuple = (58, 57)
#target_tuple = (30, 50) #The connectomes to check up on and create groupings clusters for

target_tuple = (39,32)

function_processes = parse_arguments_function(sys.argv)

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
    centroid_file_path = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '.py')
    grouping_files = {}
    exists=True
    for ref in references:
        grouping_files[ref,'lines']=(os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_' + ref + '_lines.py'))
        grouping_files[ref, 'points'] = (os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '_' + ref + '_points.py'))
        _, exists = check_files(grouping_files)
    if not os.path.exists(centroid_file_path) or np.any(exists) is False or overwrite:
    #if not os.path.exists(centroid_file_path):
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
            trkdata = load_trk(trkpath, 'same')
            header = trkdata.space_attributes
            picklepath_connectome = os.path.join(pickle_folder, subject + str_identifier + '_connectome.p')
            picklepath_grouping = os.path.join(pickle_folder, subject + str_identifier + '_grouping.p')
            M_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_connectome.xlsx")
            grouping_xlsxpath = os.path.join(excel_folder, subject + str_identifier + "_grouping.xlsx")
            #if os.path.exists(picklepath_grouping) and not overwrite:
            #    with open(picklepath_grouping, 'rb') as f:
            #        grouping = pickle.load(f)
            if os.path.exists(picklepath_connectome) and not overwrite:
                with open(picklepath_connectome, 'rb') as f:
                    M = pickle.load(f)
            if os.path.exists(grouping_xlsxpath):
                grouping = extract_grouping(grouping_xlsxpath, index_to_struct, None, verbose=verbose)
            else:
                if function_processes == 1:
                    M, grouping_temp = connectivity_matrix(trkdata.streamlines, trkdata.space_attributes[0], labelmask, inclusive=True, symmetric=True,
                                        return_mapping=True,
                                        mapping_as_streamlines=False)
                else:
                    M, grouping_temp = connectivity_matrix_func(trkdata.streamlines, function_processes, labelmask,
                                                                symmetric=True, mapping_as_streamlines=False,
                                                                affine_streams=trkdata.space_attributes[0],
                                                                inclusive=True, verbose=verbose)
                M_grouping_excel_save(M,grouping_temp,M_xlsxpath, grouping_xlsxpath, index_to_struct, verbose=False)

                if os.path.exists(grouping_xlsxpath):
                    grouping = extract_grouping(grouping_xlsxpath, index_to_struct, np.shape(M), verbose=verbose)
                else:
                    raise Exception(f'saving of the excel at {grouping_xlsxpath} did not work')

                if picklesave:
                    pickle.dump(grouping, open(picklepath_grouping, "wb"))
                    pickle.dump(M, open(picklepath_connectome, "wb"))
                    if verbose:
                        txt = ("The connectomes were saved at " + picklepath_grouping + "and " + picklepath_connectome)
                        #send_mail(txt, subject="Pickle save")
                        print(txt)

            target_streamlines_list = grouping[target_tuple[0], target_tuple[1]]
            target_streamlines = trkdata.streamlines[target_streamlines_list]
            target_streamlines_set = set_number_of_points(target_streamlines, nb_points=num_points2)
            del(target_streamlines, trkdata)
            target_qb = QuickBundles(threshold=distance1, metric=metric1)
            target_clusters = target_qb.cluster(target_streamlines_set)

            #srr = StreamlineLinearRegistration()
            #srm = srr.optimize(static=target_clusters_control.centroids, moving=target_clusters.centroids)
            #target_str_aligned = srm.transform(target_streamlines_set)

            # first clustering for transform matrix
            #target_qb = QuickBundles(threshold=distance1, metric=metric1)
            #target_clusters = target_qb.cluster(target_streamlines_set)
            #         print('NO.'+str(j+1)+' '+runno+" Nb. clusters:", len(target_clusters))

            for ref in references:
                ref_img_path = get_diff_ref(label_folder, subject, ref)
                ref_data, ref_affine = load_nifti(ref_img_path)
                stream_ref = []
                stream_point_ref = []
                for s in range(len(target_streamlines_set)):
                    point_ref = [ref_data[int(k[0]), int(k[1]), int(k[2])] for k in target_streamlines_set[s]]
                    stream_point_ref.append(point_ref)
                    stream_ref.append(np.mean(point_ref))

                groupLines[group, ref].extend(stream_ref)
                groupPoints[group, ref].extend(stream_point_ref)

            groupstreamlines[group].extend(target_streamlines_set)

        group_qb[group] = QuickBundles(threshold=distance2, metric=metric2)
        group_clusters[group] = group_qb[group].cluster(groupstreamlines[group])
        if os.path.exists(centroid_file_path) and overwrite:
            os.remove(centroid_file_path)
        if not os.path.exists(centroid_file_path):
            if verbose:
                print(f'Summarized the clusters for group {group} at {centroid_file_path}')
            pickle.dump(group_clusters[group], open(centroid_file_path, "wb"))

        for ref in references:
            if overwrite:
                os.remove([grouping_files[ref,'lines'], grouping_files[ref,'points']])
            if not os.path.exists(grouping_files[ref,'lines']):
                if verbose:
                    print(f"Summarized the clusters for group {group} and statistics {ref} at {grouping_files[ref,'lines']}")
                pickle.dump(groupLines[group, ref], grouping_files[ref,'lines'])
                pickle.dump(groupPoints[group, ref], grouping_files[ref,'points'])

    else:
        print(f'Centroid file was found at {centroid_file_path}')
        with open(centroid_file_path, 'rb') as f:
            group_clusters[group] = pickle.load(f)
        for ref in references:
            ref_path_lines = grouping_files[ref, 'lines']
            groupLines[group, ref] = grouping_files[ref, 'lines']
            ref_path_points = grouping_files[ref, 'points']
            groupPoints[group, ref] = grouping_files[ref, 'points']

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

toview=True
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

"""
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
