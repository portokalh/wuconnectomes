from dipy.io.streamline import load_trk, save_trk
from dipy.viz import window, actor
import os
import pickle
from tract_visualize import show_bundles, setup_view
from convert_atlas_mask import convert_labelmask, atlas_converter
from tract_handler import ratio_to_str, gettrkpath
from itertools import compress
import numpy as np
import nibabel as nib, socket
from file_tools import mkcdir
from streamline_nocheck import load_trk as load_trk_spe
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric
import warnings
from dipy.align.streamlinear import StreamlineLinearRegistration
import copy

"""
# ,(23,30)
target_tuples = [(9, 1), (24, 1), (22, 1), (58, 57), (64, 57)]
target_tuples = [(9, 1), (24, 1), (22, 1), (58, 57), (23, 24), (64, 57)]
target_tuples = [(58, 57), (9, 1), (24, 1), (22, 1), (64, 57), (23, 24), (24, 30), (23, 30)]
target_tuples = [(24, 30), (23, 24)]
target_tuples = [(80, 58)]
"""

computer_name = socket.gethostname()

project = 'AD_Decode'
#genotype_noninclusive
#target_tuples = [(9, 1), (24, 1), (58, 57), (64, 57), (22, 1)]
#genotype_noninclusive_volweighted_fa
#target_tuples = [(9, 1), (57, 9), (61, 23), (84, 23), (80, 9)]

#sex_noninclusive
#target_tuples = [(64, 57), (58, 57), (9, 1), (64, 58), (80,58)]
#target_tuples = [(64,57)]
#sex_noninclusive_volweighted_fa
#target_tuples = [(58, 24), (58, 30), (64, 30), (64, 24), (58,48)]

#target_tuples = [(9,1)]
#groups = ['APOE4', 'APOE3']


fixed = True
record = ''

inclusive = False
symmetric = True
write_txt = True
ratio = 100
top_percentile = 100
num_bundles = 20

selection = 'num_streams'
coloring = 'bundles_coloring'
references = ['fa','md']
references = ['fa']
cutoffref = 0
groups = ['Male','Female']
groups = ['APOE4', 'APOE3']
non_control = groups[0]
control = groups[1]
write_stats = False
registration = False
overwrite = True

changewindow_eachtarget = False

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

# if fixed:
#    fixed_str = '_fixed'
# else:
#    fixed_str = ''

samos = False
if 'samos' in computer_name:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in computer_name:
    # mainpath = '/Users/alex/jacques/'
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    raise Exception('No other computer name yet')


ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio', '')

_, _, index_to_struct, _ = atlas_converter(ROI_legends)

if project == 'AMD':
    mainpath = os.path.join(mainpath, project)
    groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
              'Paired Initial AMD']
    anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'

if project == 'AD_Decode':
    mainpath = os.path.join(mainpath, project, 'Analysis')
    anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_b0.nii.gz'

# figures_path = '/Volumes/Data/Badea/Lab/human/AMD/Figures_MDT_non_inclusive/'
# centroid_folder = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/'
figures_path = os.path.join(mainpath, f'Figures_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
centroid_folder = os.path.join(mainpath, f'Centroids_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
trk_folder = os.path.join(mainpath, f'Centroids_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
stats_folder = os.path.join(mainpath, f'Statistics_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')

mkcdir([figures_path, centroid_folder, stats_folder])

# groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
#          'Paired Initial AMD']

# anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'


# superior frontal right to cerebellum right

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

scene = None



for target_tuple in target_tuples:

    interactive = True

    print(target_tuple[0], target_tuple[1])
    region_connection = index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]]
    print(region_connection)

    if write_txt:
        text_path = os.path.join(figures_path, region_connection + '_stats.txt')
        testfile = open(text_path, "w")
        testfile.write("Parameters for groups\n")
        testfile.close()

    if changewindow_eachtarget:
        firstrun = True

    selected_bundles = {}
    selected_centroids = {}
    selected_sizes = {}
    streamlines = {}
    num_bundles_group = {}

    ref_lines = {}
    ref_points = {}

    for group in groups:

        selected_bundles[group] = []
        selected_centroids[group] = []
        selected_sizes[group] = []

        print(f'Setting up group {group}')
        group_str = group.replace(' ', '_')

        stats_path = os.path.join(stats_folder,
                                  group_str + '_MDT' + ratio_str + '_' + region_connection + '_bundle_stats.xlsx')

        if not os.path.exists(stats_path) or overwrite:
            if os.path.exists(stats_path):
                os.remove(stats_path)
            import xlsxwriter
            workbook = xlsxwriter.Workbook(stats_path)
            worksheet = workbook.add_worksheet()
            l=1
            worksheet.write(0,l,'Number streamlines')
            l+=1
            for ref in references:
                worksheet.write(0,l, ref + ' mean')
                worksheet.write(0,l+1, ref + ' min')
                worksheet.write(0,l+2, ref + ' max')
                worksheet.write(0,l+3, ref + ' std')
                l=l+4
        else:
            print(f'The file {stats_path} already exists and no overwrite enabled: skipping')
            #continue

        centroid_file_path = os.path.join(centroid_folder,
                                          group_str + '_MDT' + ratio_str + '_' + region_connection + '_centroid.py')

        trk_path = os.path.join(trk_folder,
                                group_str + '_MDT' + ratio_str + '_' + region_connection + '_streamlines.trk')

        # '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Centroids_MDT_non_inclusive_symmetric_100/APOE4_MDT_ratio_100_ctx-lh-inferiorparietal_left_to_ctx-lh-inferiortemporal_left_streamlines.trk'

        for ref in references:
            ref_path_lines = os.path.join(centroid_folder,
                                   group_str + '_MDT' + ratio_str + '_' + region_connection + f'_{ref}_lines.py')
            ref_path_points = os.path.join(centroid_folder,
                                   group_str + '_MDT' + ratio_str + '_' + region_connection + f'_{ref}_points.py')

            if os.path.exists(ref_path_points):
                with open(ref_path_points, 'rb') as f:
                    ref_points[group,ref] = pickle.load(f)
            else:
                txt = f'Could not find file {ref_path_points} for group {group} reference {ref}'
                raise Exception(txt)

            if os.path.exists(ref_path_lines):
                with open(ref_path_lines, 'rb') as f:
                    ref_lines[group,ref] = pickle.load(f)
            else:
                txt = f'Could not find file {ref_path_lines} for group {group} reference {ref}'
                raise Exception(txt)

        if os.path.exists(trk_path):
            try:
                streamlines_data = load_trk(trk_path, 'same')
            except:
                streamlines_data = load_trk_spe(trk_path, 'same')
        streamlines[group] = streamlines_data.streamlines

        if top_percentile<100:
            cutoff = np.percentile(ref_lines[group,references[cutoffref]], 100 - top_percentile)
            select_streams = ref_lines[group,references[cutoffref]] > cutoff
            streamlines[group] = list(compress(streamlines[group], select_streams))
            streamlines[group] = nib.streamlines.ArraySequence(streamlines[group])

            for ref in references:
                if np.shape(streamlines[group])[0] != np.shape(ref_lines[group][ref])[0]:
                    raise Exception('Inconsistency between streamlines and fa lines')
                ref_lines[group,ref] = list(compress(ref_lines[group,ref], select_streams))

        group_qb = QuickBundles(threshold=distance2, metric=metric2)
        group_clusters = group_qb.cluster(streamlines[group])
        #group2_qb = QuickBundles(threshold=distance2, metric=metric2)
        #group2_clusters = group2_qb.cluster(groupstreamlines2)

        num_bundles_group[group]=0
        if selection =='num_streams':
            num_streamlines = [np.shape(cluster)[0] for cluster in group_clusters.clusters]
            num_streamlines = group_clusters.clusters_sizes()
            top_bundles = sorted(range(len(num_streamlines)), key=lambda i: num_streamlines[i], reverse=True)[:]
        for bundle in top_bundles:
            selected_bundles[group].append(group_clusters.clusters[bundle])
            selected_centroids[group].append(group_clusters.centroids[bundle])
            selected_sizes[group].append(group_clusters.clusters_sizes()[bundle])
            num_bundles_group[group]+=1
        bun_num = 0

        bundles_ref = {}
        bundles_ref_mean = {}
        for ref in references:
            bundles_ref[ref] = []
            bundles_ref_mean[ref] = []

        for bundle in selected_bundles[group]:
            for ref in references:
                bundle_ref = []
                for idx in bundle.indices:
                    bundle_ref.append(ref_lines[group,ref][idx])
                bundles_ref[ref].append(bundle_ref)
                bundles_ref_mean[ref].append(np.mean(bundle_ref))

        empty_bundles = {}
        for ref in references:
            empty_bundles[ref] = 0

        if write_stats:
            bun_num=0
            for bundle in top_bundles:
                l=0
                worksheet.write(bun_num+1, l, bun_num+1)
                l+=1
                worksheet.write(bun_num + 1, l, np.shape(group_clusters.clusters[bundle])[0])
                l+=1
                for ref in references:
                    if np.mean(bundles_ref[ref][bun_num])==0:
                        empty_bundles[ref] += 1
                    worksheet.write(bun_num+1, l+0, np.mean(bundles_ref[ref][bun_num]))
                    worksheet.write(bun_num+1, l+1, np.min(bundles_ref[ref][bun_num]))
                    worksheet.write(bun_num+1, l+2, np.max(bundles_ref[ref][bun_num]))
                    worksheet.write(bun_num+1, l+3, np.std(bundles_ref[ref][bun_num]))
                    l = l + 4
                bun_num+=1
            workbook.close()
        for ref in references:
            if empty_bundles[ref]>0:
                print(f'Found {empty_bundles} empty bundles out of {np.size(top_bundles)} for {ref} in group {group} for {region_connection}')


    if registration:
        srr = StreamlineLinearRegistration()
        for streamline,i in enumerate(selected_centroids[non_control]):
            srm = srr.optimize(static=selected_centroids[control], moving=streamline)
            streamlines[control][i] = srm.transform(streamline)

    from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric, mdf

    #dist_all = np.zeros((np.size(selected_bundles[control]), np.size(selected_bundles[non_control])))
    dist_all = np.zeros((num_bundles, num_bundles))


    """
    top_idx_group_control = sorted(range(len(selected_sizes[control])),
                            key=lambda i: selected_sizes[group][i], reverse=True)[:num_bundles]
    top_idx_group_noncontrol = sorted(range(len(selected_sizes[non_control])),
                            key=lambda i: selected_sizes[group][i], reverse=True)[:num_bundles]
                            
    if not np.all(top_idx_group_control == np.arange(num_bundles)) or not np.all(top_idx_group_noncontrol == np.arange(num_bundles)):
        warnings.warn('There is indeed a difference between the two')
    else:
        print('no difference between the two')
    """

    for g3 in np.arange(num_bundles):
        for g4 in np.arange(num_bundles):
            dist_all[g3, g4] = (mdf(selected_centroids[control][g3], selected_centroids[non_control][g4]))

    dist_all_fix = copy.copy(dist_all)
    dist_all_idx = []
    #for i in range(len(selected_centroids[group])):
    for i in np.arange(num_bundles):
        idx = np.argmin(dist_all_fix[i, :])
        dist_all_idx.append([i, idx])
        dist_all_fix[:, idx] = 100000

    dist_group3_idx = [dist_all_idx[iii][0] for iii in range(num_bundles)]  # size id
    dist_group4_idx = [dist_all_idx[iii][1] for iii in range(num_bundles)]  # size id

    group_list = {}
    dist_idx = {}
    for j,group in enumerate(groups):
        dist_idx[group] = [dist_all_idx[iii][j] for iii in range(num_bundles)]
        group_list[group]=([np.arange(num_bundles)[dist_all_idx[i][j]] for i in range(num_bundles)])

    import pandas as pd
    from dipy.tracking import utils

    num_bundles_full_stats = 10

    for group in groups:
        groupcsv = np.zeros((1, 5+np.size(references)))
        references_string = "_".join(references)
        csv_summary = os.path.join(stats_folder, group + '_' + region_connection + ratio_str + f'_bundle_stats_{references_string}.csv')
        if not os.path.exists(csv_summary) or overwrite:
            for i in range(num_bundles_full_stats):
                idsize = dist_idx[group][i]
                idbundle = group_list[group][i]
                fa = []
                for s in selected_bundles[group][idbundle].indices:
                    #temp = np.hstack((idsize * np.ones((num_points2, 1)),
                    #                  idbundle * np.ones((num_points2, 1)),
                    #                  s * np.ones((num_points2, 1)),
                    #                  np.array(range(num_points2)).reshape(num_points2, 1),
                    #                  list(utils.length([streamlines[group][s]])) * np.ones((num_points2, 1)),
                    #                  np.array(ref_points[group, ref][s]).reshape(num_points2, 1)))
                    temp = np.hstack((idsize * np.ones((num_points2, 1)),
                                      idbundle * np.ones((num_points2, 1)),
                                      s * np.ones((num_points2, 1)),
                                      np.array(range(num_points2)).reshape(num_points2, 1),
                                      list(utils.length([streamlines[group][s]])) * np.ones((num_points2, 1))))
                    for ref in references:
                        temp = np.hstack((temp,np.array(ref_points[group, ref][s]).reshape(num_points2, 1)))
                    groupcsv = np.vstack((groupcsv, temp))
                    groupcsv = groupcsv[1:, :]
            groupcsvDF = pd.DataFrame(groupcsv)
            groupcsvDF.rename(index=str, columns={0: "Bundle Size Rank", 1: "Bundle ID", 2: "Steamlines ID",
                                                   3: "Point ID", 4: "length"})
            for i,ref in enumerate(references):
                groupcsvDF.rename(index=str, columns={5+i: ref})
            print('writing')
            groupcsvDF.to_csv(csv_summary, header=["Bundle Size Rank", "Bundle ID", "Streamlines ID", "Point ID", "Length"] + references)
            print(f'Writing bundle stats for {group} and {region_connection} to {csv_summary}')
        else:
            print(f'The file {csv_summary} already exists and no overwrite enabled: skipping')



    """
    # color by line-average fa
    renderer = window.Renderer()
    renderer.clear()
    renderer = window.Renderer()
    stream_actor3 = actor.line(group_clusters.clusters[bundle_id], np.array(bundle_fa), lookup_colormap=cmap)
    renderer.add(stream_actor3)
    bar = actor.scalar_bar(cmap)
    renderer.add(bar)
    # Uncomment the line below to show to display the window
    window.show(renderer, size=(600, 600), reset_camera=False)
    """
