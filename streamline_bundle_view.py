from dipy.io.streamline import load_trk, save_trk
from dipy.viz import window, actor
import os
import pickle
from tract_visualize import show_bundles, setup_view, view_test, setup_view_colortest
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
from dipy.denoise.enhancement_kernel import EnhancementKernel
#from dipy.tracking.fbcmeasures import FBCMeasures
import fury
#import pyximport
#pyximport.install()
import pandas as pd
from fbcmeasures import FBCMeasures

computer_name = socket.gethostname()

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

#project = 'AD_Decode'
project = 'AMD'

fixed = True
record = ''


inclusive = False
symmetric = True
write_txt = True
ratio = 1
top_percentile = 100
num_bundles = 10
num_bundles_toview = 4
distance = 3
num_points = 50

if project == 'AD_Decode':
    # ,(23,30)
    target_tuples = [(9, 1), (24, 1), (22, 1), (58, 57), (64, 57)]
    target_tuples = [(9, 1), (24, 1), (22, 1), (58, 57), (23, 24), (64, 57)]
    target_tuples = [(58, 57), (9, 1), (24, 1), (22, 1), (64, 57), (23, 24), (24, 30), (23, 30)]
    target_tuples = [(24, 30), (23, 24)]
    target_tuples = [(80, 58)]
    target_tuples = [(58, 57)]
    target_tuples = [(64,57)]
    #genotype_noninclusive
    target_tuples = [(9, 1), (24, 1), (58, 57), (64, 57), (22, 1)]
    #target_tuples = [(24, 1)]
    #genotype_noninclusive_volweighted_fa
    #target_tuples = [(9, 1), (57, 9), (61, 23), (84, 23), (80, 9)]

    #sex_noninclusive
    #target_tuples = [(64, 57), (58, 57), (9, 1), (64, 58), (80,58)]
    #sex_noninclusive_volweighted_fa
    #target_tuples = [(58, 24), (58, 30), (64, 30), (64, 24), (58,48)]

    groups = ['APOE4', 'APOE3']
    #groups = ['Male','Female']

    mainpath = os.path.join(mainpath, project, 'Analysis')
    anat_path = os.path.join(mainpath,'../../mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_b0.nii.gz')
    space_param = '_MDT'

if project == 'AMD':
    #target_tuples = [(9, 1),(24, 1), (76, 42), (76, 64), (77, 9), (43, 9)]
    groups_all = ['Paired 2-YR AMD','Paired 2-YR Control','Paired Initial Control','Paired Initial AMD',
                  'Initial AMD', 'Initial Control']
    groups_set = {'Initial':[2,3],'2Year':[0,1]}
    #groups = ['Paired Initial Control', 'Paired Initial AMD', 'Paired 2-YR Control', 'Paired 2-YR AMD']
    #groups = ['Paired Initial Control','Paired Initial AMD']
    #groups = ['Paired 2-YR Control', 'Paired 2-YR AMD']
    mainpath = os.path.join(mainpath, project)
    anat_path = os.path.join(mainpath,'../../mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz')

    space_param = '_affinerigid'

    #TN-PCA
    #target_tuples = [(62, 28), (56, 45), (77, 43), (58, 45), (79, 45), (56, 50), (28, 9), (62, 1), (28, 1), (62, 9), (22, 9), (56, 1),(77, 43), (76, 43), (61, 29), (63, 27), (73, 43), (53, 43)]
    #VBA
    #target_tuples = [(27,29), (61,63),(30, 16), (24, 16),(28, 31), (28, 22),(22, 31)]
    #TN-PCA / VBA combination
    #target_tuples = [(62, 28), (58, 45), (28, 9), (62, 1), (77, 43), (61, 29)]
    #target_tuples = [(62, 28), (58, 45),(77, 43), (61, 29)]

    group_select = 'Initial'
    groups = [groups_all[x] for x in groups_set[group_select]]

    target_tuples = [(62, 28), (58, 45),(77, 43), (61, 29)]
    target_tuples = [(61, 29)]
    plane = 'x'


write_stats = False

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



# target_tuple = (24,1)
# target_tuple = [(58, 57)]
# target_tuples = [(64, 57)]


ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio', '')
# target_tuple = (9,77)

_, _, index_to_struct, _ = atlas_converter(ROI_legends)

# figures_path = '/Volumes/Data/Badea/Lab/human/AMD/Figures_MDT_non_inclusive/'
# centroid_folder = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/'

figures_path = os.path.join(mainpath, f'Figures_distance10{space_param}{inclusive_str}{symmetric_str}{folder_ratio_str}')
centroid_folder = os.path.join(mainpath, f'Centroids{space_param}{inclusive_str}{symmetric_str}{folder_ratio_str}')
trk_folder = os.path.join(mainpath, f'Centroids{space_param}{inclusive_str}{symmetric_str}{folder_ratio_str}')
stats_folder = os.path.join(mainpath, f'Statistics{space_param}{inclusive_str}{symmetric_str}{folder_ratio_str}')

mkcdir([figures_path, centroid_folder, stats_folder])

# groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
#          'Paired Initial AMD']

# anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'


# superior frontal right to cerebellum right

#group cluster parameter

feature = ResampleFeature(nb_points=num_points)
metric = AveragePointwiseEuclideanMetric(feature=feature)


scene = None
selection = 'num_streams'

coloring_options = ['centroids_fa_mean_coloring','centroids_id_coloring','streams_fa_mean_coloring',
                    'streams_fa_points_coloring','streams_id_coloring','coherence_coloring_bundle','coherence_coloring_streams','coherence_coloring_points']
#coloring_options = ['streams_fa_points_coloring','streams_id_coloring']
#coloring_options = ['streams_fa_points_coloring','coherence_coloring_bundle','coherence_coloring_streams','coherence_coloring_points']
#coloring_options = ['coherence_coloring_bundle','coherence_coloring_streams','coherence_coloring_points']
coloring_options = ['streams_fa_points_coloring','streams_fa_mean_coloring','coherence_coloring_bundle','coherence_coloring_streams']
#coloring_options = ['streams_id_coloring']
fa_scale_range = (0.1, 0.3)
coherence_scale_range = (0.1, 0.6)

for coloring_option in coloring_options:
    if 'coherence' in coloring_option:
        D33 = 1
        D44 = 0.02
        t = 1
        k = EnhancementKernel(D33, D44, t)
        break
#coloring = 'bundles_fa_coloring'
#coloring = 'bundles_id_coloring'
#coloring = 'streams_fa_mean_coloring'
coloring = 'streams_fa_points_coloring'
coloring = 'streams_id_coloring'
#coloring = 'centroids_fa_mean_coloring'
#coloring = 'centroids_fa_point_coloring'

references = ['fa']




firsttest = True

interactive=True

for target_tuple in target_tuples:

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

    for group in groups:

        print(f'Setting up group {group}')
        group_str = group.replace(' ', '_')

        group_connection_str = group_str + space_param + ratio_str + '_' + region_connection
        if write_stats:
            stats_path = os.path.join(stats_folder, group_connection_str + '_bundle_stats.xlsx')
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

        centroid_file_path = os.path.join(centroid_folder,
                                          group_connection_str + '_centroid.py')
        fa_path = os.path.join(centroid_folder, group_connection_str + '_fa_lines.py')
        md_path = os.path.join(centroid_folder, group_connection_str + '_md_lines.py')
        trk_path = os.path.join(trk_folder, group_connection_str + '_streamlines.trk')
        fa_points_path = (os.path.join(centroid_folder, group_connection_str + '_' + 'fa' + '_points.py'))

        if os.path.exists(fa_path):
            with open(fa_path, 'rb') as f:
                fa_lines = pickle.load(f)
        if os.path.exists(md_path):
            with open(md_path, 'rb') as f:
                md_lines = pickle.load(f)
        # '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/Centroids_MDT_non_inclusive_symmetric_100/APOE4_MDT_ratio_100_ctx-lh-inferiorparietal_left_to_ctx-lh-inferiortemporal_left_streamlines.trk'
        if os.path.exists(trk_path):
            try:
                streamlines_data = load_trk(trk_path, 'same')
            except:
                streamlines_data = load_trk_spe(trk_path, 'same')


        #streamlines_2 = streamlines_data.remove_invalid_streamlines()

        streamlines = streamlines_data.streamlines

        if 'fa_lines' in locals():
            cutoff = np.percentile(fa_lines, 100 - top_percentile)
            select_streams = fa_lines >= cutoff
            fa_lines = list(compress(fa_lines, select_streams))
            streamlines = list(compress(streamlines, select_streams))
            streamlines = nib.streamlines.ArraySequence(streamlines)
            if np.shape(streamlines)[0] != np.shape(fa_lines)[0]:
                raise Exception('Inconsistency between streamlines and fa lines')
        else:
            txt = f'Cannot find {fa_path}, could not select streamlines based on fa'
            warnings.warn(txt)
            fa_lines = [None]

        group_qb = QuickBundles(threshold=distance, metric=metric)
        group_clusters = group_qb.cluster(streamlines)

        selected_bundles = []
        if selection =='num_streams':
            num_streamlines = [np.shape(cluster)[0] for cluster in group_clusters.clusters]
            top_bundles = sorted(range(len(num_streamlines)), key=lambda i: num_streamlines[i], reverse=True)[:num_bundles]
        for bundle in top_bundles:
            selected_bundles.append(group_clusters.clusters[bundle])
        bun_num = 0

        #colors_list = [window.colors.green, window.colors.yellow, window.colors.red, window.colors.brown,
        #          window.colors.orange, window.colors.blue, window.colors.pink, window.colors.violet,
        #          window.colors.cyan, window.colors.purple]

        for coloring in coloring_options:
            colorbar = True

            figures_coloring_path = os.path.join(mainpath, figures_path, coloring)
            #figures_coloring_path = os.path.join(mainpath, figures_path, 'test_zone')
            mkcdir(figures_coloring_path)

            if coloring == 'coherence_coloring_points':
                from dipy.io.image import load_nifti

                fbc_bundles = []
                bundles_clrs_points = []
                for bundle in selected_bundles:
                    fbc = FBCMeasures(streamlines[bundle.indices], k)
                    fbc_sl_orig, lfbc_orig, rfbc_orig = \
                        fbc.get_points_rfbc_thresholded(0, emphasis=0.01)
                    fbc_bundles.append(nib.streamlines.ArraySequence(fbc_sl_orig))
                    bundle_clrs_points = []
                    for stream_colors in lfbc_orig:
                        #for points_colors in stream_colors:
                        bundle_clrs_points.append(stream_colors)
                    bundles_clrs_points.append(bundle_clrs_points)

                coloring_vals = bundles_clrs_points
                trkobject = fbc_bundles
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=coherence_scale_range)

            elif coloring == 'coherence_coloring_streams':

                from dipy.tracking.fbcmeasures import FBCMeasures
                from dipy.io.image import load_nifti

                fbc_bundles = []
                bundles_clrs_streams = []
                for bundle in selected_bundles:
                    fbc = FBCMeasures(streamlines[bundle.indices], k)
                    fbc_sl_orig, lfbc_orig, rfbc_orig = \
                        fbc.get_points_rfbc_thresholded(0, emphasis=0.01)
                    fbc_bundles.append(nib.streamlines.ArraySequence(fbc_sl_orig))
                    bundles_clrs_streams.append(rfbc_orig)

                coloring_vals = bundles_clrs_streams
                trkobject = fbc_bundles
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=coherence_scale_range)

            elif coloring == 'coherence_coloring_bundle':

                from dipy.tracking.fbcmeasures import FBCMeasures
                from dipy.io.image import load_nifti

                fbc_bundles = []
                bundles_clrs_bundle = []
                for bundle in selected_bundles:
                    fbc = FBCMeasures(streamlines[bundle.indices], k)
                    fbc_sl_orig, lfbc_orig, rfbc_orig = \
                        fbc.get_points_rfbc_thresholded(0, emphasis=0.01)

                    fbc_bundles.append(nib.streamlines.ArraySequence(fbc_sl_orig))
                    templist=[]
                    templist.extend(np.mean(rfbc_orig) for i in range(len(fbc_sl_orig)))
                    bundles_clrs_bundle.append(templist)

                coloring_vals = bundles_clrs_bundle
                trkobject = fbc_bundles
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=coherence_scale_range)


            elif coloring == 'centroids_fa_point_coloring':
                if os.path.exists(fa_points_path):
                    with open(fa_points_path, 'rb') as f:
                        fa_points = pickle.load(f)
                if 'select_streams' in locals():
                    fa_points = list(compress(fa_points, select_streams))
                bundles_fa = []
                bundles_fa_mean = []
                bundle_fa_points_1 = []
                bundle_fa_points_2 = []
                for bundle in selected_bundles:
                    bundle_fa = []
                    for idx in bundle.indices:
                        bundle_fa.append(fa_points[idx])
                        for idx_point in range(len(fa_points[idx])):
                            bundle_fa_points_2.append(fa_points[idx][idx_point])

                    bundle_fa_points_1.append(np.array(bundle_fa))

                coloring_vals = bundle_fa_points_1
                trkobject = selected_bundles
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=fa_scale_range)

            elif coloring == 'centroids_fa_mean_coloring':
                bundles_fa = []
                bundles_fa_mean = []
                for bundle in selected_bundles:
                    bundle_fa = []
                    for idx in bundle.indices:
                        bundle_fa.append(fa_lines[idx])
                    bundles_fa.append(bundle_fa)
                    bundles_fa_mean.append(np.mean(bundle_fa))
                coloring_vals = bundles_fa_mean
                trkobject = selected_bundles
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=fa_scale_range)

            elif coloring == 'centroids_id_coloring':
                bundles_fa = []
                bundles_fa_mean = []
                for bundle in selected_bundles:
                    bundle_fa = []
                    for idx in bundle.indices:
                        bundle_fa.append(fa_lines[idx])
                    bundles_fa.append(bundle_fa)
                    bundles_fa_mean.append(np.mean(bundle_fa))
                trkobject = selected_bundles
                coloring_vals = fury.colormap.distinguishable_colormap(nb_colors=np.size(selected_bundles))

                if np.size(selected_bundles)>np.size(coloring_vals):
                    raise Exception('Not enough colors for number of bundles')
                else:
                    coloring_vals = coloring_vals[:np.size(selected_bundles)]
                lut_cmap = None

            elif coloring == 'streams_fa_mean_coloring':
                bundle_streamlines = []
                bundles_fa = []
                bundles_fa_mean = []
                for bundle in selected_bundles:
                    bundle_streamlines.append(streamlines[bundle.indices])
                    bundle_fa = []
                    for idx in bundle.indices:
                        bundle_fa.append(fa_lines[idx])
                    bundles_fa.append(bundle_fa)
                    bundles_fa_mean.append(np.mean(bundle_fa))
                coloring_vals = bundles_fa
                trkobject = bundle_streamlines
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=fa_scale_range)

            elif coloring == 'streams_fa_points_coloring':
                if os.path.exists(fa_points_path):
                    with open(fa_points_path, 'rb') as f:
                        fa_points = pickle.load(f)
                if 'select_streams' in locals():
                    fa_points = list(compress(fa_points, select_streams))
                bundle_fa_points = []
                bundle_streamlines = []
                for bundle in selected_bundles:
                    bundle_streamlines.append(streamlines[bundle.indices])
                    bundle_fa = []
                    for idx in bundle.indices:
                        bundle_fa.append(fa_points[idx])
                    bundle_fa_points.append(bundle_fa)
                coloring_vals = bundle_fa_points
                trkobject = bundle_streamlines
                lut_cmap = actor.colormap_lookup_table(
                    scale_range=fa_scale_range)

            elif coloring == 'streams_id_coloring':
                csv_bundleorder = os.path.join(stats_folder,
                                               group_select + '_' + region_connection + ratio_str + f'_bundle_order.csv')
                bundleorder = pd.read_csv(csv_bundleorder)
                neworder = bundleorder.to_dict()[group]

                bundle_streamlines = []
                neworder2 = np.zeros([num_bundles, 1])
                for i in np.arange(num_bundles):
                    neworder2[i] = neworder[i]

                selected_bundles = [selected_bundles[int(i)] for i in neworder2]

                for bundle in selected_bundles:
                    bundle_streamlines.append(streamlines[bundle.indices])
                coloring_vals = fury.colormap.distinguishable_colormap(nb_colors=np.size(selected_bundles))
                trkobject = bundle_streamlines[:num_bundles_toview]
                if np.size(selected_bundles)>len(coloring_vals):
                    raise Exception('Not enough colors for number of bundles')
                else:
                    coloring_vals = coloring_vals[:num_bundles_toview]
                lut_cmap = None
                colorbar=False

            else:
                trkobject = selected_bundles
                bundles_fa = None

            if write_stats:
                bun_num=0
                for bundle in selected_bundles:
                    l=0
                    worksheet.write(bun_num+1, l, bun_num+1)
                    l+=1
                    worksheet.write(bun_num + 1, l, np.shape(bundle)[0])
                    l+=1
                    for ref in references:
                        worksheet.write(bun_num+1, l+0, np.mean(bundles_fa[bun_num]))
                        worksheet.write(bun_num+1, l+1, np.min(bundles_fa[bun_num]))
                        worksheet.write(bun_num+1, l+2, np.max(bundles_fa[bun_num]))
                        worksheet.write(bun_num+1, l+3, np.std(bundles_fa[bun_num]))
                        l = l + 4
                    bun_num+=1
                workbook.close()

            record_path = os.path.join(figures_coloring_path, group_connection_str + f'_bundles_figure_distance_{str(distance)}.png')
            #scene = None
            #interactive = False
            #record_path = None
            scene = setup_view(trkobject, colors = lut_cmap,ref = anat_path, world_coords = True, objectvals = coloring_vals, colorbar=colorbar, record = record_path, scene = scene, plane = plane, interactive = interactive)
            #scene = setup_view_colortest(trkobject, colors=lut_cmap, ref=anat_path, world_coords=True, objectvals=coloring_vals,
            #                   colorbar=True, record=record_path, scene=scene, plane=plane, interactive=interactive)
            """
            test='record'
            if test is not None and firsttest:
                record_path = os.path.join(figures_coloring_path,'quick_test_2.png')
                view_test(scene,test,record_path=record_path)
                firsttest = False
            """
            interactive = False

        del(fa_lines,streamlines)

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