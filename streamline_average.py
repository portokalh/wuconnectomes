
import numpy as np
import copy
import pandas as pd
from nibabel import trackvis as tv
from dipy.tracking.streamline import Streamlines
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
from dipy.data import get_fnames
import vtk
from dipy.viz import window, actor
from dipy.tracking.streamline import Streamlines
from dipy.io.streamline import load_trk, save_trk
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric,mdf
from dipy.io.image import load_nifti
import warnings
from dipy.tracking import utils
from dipy.viz import window, actor
from time import sleep
from dipy.data import two_cingulum_bundles
from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines
import os, glob
from tract_save import unload_trk
import pickle
from dipy.tracking.utils import connectivity_matrix
from nifti_handler import getfa, getdiffdata_all, getdiffdata, getdiffpath, getgtab, getlabelmask, move_bvals, getmask, getb0s, getlabeltypemask
from file_tools import mkcdir
from tract_handler import ratio_to_str, gettrkpath
from convert_atlas_mask import convert_labelmask, atlas_converter


def show_bundles(bundles, colors=None, show=True, fname=None, fa=False, str_tube=False):
    ren = window.Renderer()
    ren.SetBackground(1., 1, 1)
    if str_tube:
        bundle_actor = actor.streamtube(bundles, colors, linewidth=0.5)
        ren.add(bundle_actor)
    else:
        for (i, bundle) in enumerate(bundles):
            color = colors[i]
            #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05

            lines_actor = actor.line(bundle, color, linewidth=2.5)
            # lines_actor.RotateX(-90)
            # lines_actor.RotateZ(90)
            ren.add(lines_actor)

    if fa:
        fa, affine_fa = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
        fa_actor = actor.slicer(fa, affine_fa)
        ren.add(fa_actor)

    if show:
        window.show(ren)
    if fname is not None:
        sleep(1)
        window.record(ren, n_frames=1, out_path=fname, size=(900, 900))

def get_grouping(grouping_xlsx):
    print('not done yet')

def get_diff_ref():
    print('not done yet')

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


samos = False
if samos:
    TRK_folder = '/mnt/paros_MRI/jacques/AMD/TRK_MDT_fixed'
    label_folder = '/mnt/paros_MRI/jacques/AMD/DWI'
    trkpaths = glob.glob(os.path.join(TRK_folder,'*trk'))
    figures_folder = '/mnt/paros_MRI/jacques/AMD/Figures_MDT'

santorini = True
if santorini:
    TRK_folder = '/Users/alex/jacques/AMD_testing_zone/TRK_MDT'
    label_folder = '/Users/alex/jacques/AMD_testing_zone/DWI'
    trkpaths = glob.glob(os.path.join(TRK_folder,'*trk'))
    figures_folder = '/Users/alex/jacques/AMD_testing_zone/Figures_MDT'
    pickle_folder = '/Users/alex/jacques/AMD_testing_zone/Pickle_MDT'
    centroid_folder = '/Users/alex/jacques/AMD_testing_zone/Centroids_MDT'
    mkcdir(figures_folder)
    mkcdir(pickle_folder)
    mkcdir(centroid_folder)
    if not os.path.exists(TRK_folder):
        raise Exception(f'cannot find TRK folder at {TRK_folder}')

ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"

str_identifier = ''
#reference_img refers to statistical values that we want to compare to the streamlines, say fa, rd, etc
references = []
stream_point = {}
stream = {}
verbose = True

groupstreamlines={}
groupLines = {}
groupPoints = {}
#locals()['groupstreamlines' + str(k + 3)].extend(target_str_aligned)
#locals()['Nativegroupstreamlines' + str(k + 3)].extend(native_target_stream_aligned)
#locals()['groupLinesFA' + str(k + 3)].extend(stream_fa)
#locals()['groupPointsFA' + str(k + 3)].extend(stream_point_fa)

groups_streamlines = {}

groups_subjects = {}
groups_subjects['Initial AMD'] = ['H27640', 'H27778', 'H29020', 'H26637', 'H27680', 'H26765', 'H27017', 'H26880', 'H28308', 'H28433', 'H28338', 'H26660', 'H28809', 'H27610', 'H26745', 'H27111', 'H26974', 'H27391', 'H28748', 'H29025', 'H29013', 'H27381', 'H26958', 'H28662', 'H26578', 'H28698', 'H27495', 'H28861', 'H28115', 'H28437', 'H26850', 'H28532', 'H28377', 'H28463', 'H26890', 'H28373', 'H28857', 'H27164', 'H27982']
groups_subjects['Paired 2-YR AMD'] = ['H22825', 'H21850', 'H29225', 'H29304', 'H29060', 'H23210', 'H21836', 'H29618', 'H22644', 'H22574', 'H22369', 'H29627', 'H29056', 'H22536', 'H23143', 'H22320', 'H22898', 'H22864', 'H29264', 'H22683']
groups_subjects['Initial Control'] = ['H26949', 'H27852', 'H28029', 'H26966', 'H27126', 'H28068', 'H29161', 'H28955', 'H26862', 'H28262', 'H28856', 'H27842', 'H27246', 'H27869', 'H27999', 'H29127', 'H28325', 'H26841', 'H29044', 'H27719', 'H27100', 'H29254', 'H27682', 'H29002', 'H29089', 'H29242', 'H27488', 'H27841', 'H28820', 'H27163', 'H28869', 'H28208', 'H27686']
groups_subjects['Paired 2-YR Control'] = ['H29403', 'H22102', 'H29502', 'H22276', 'H29878', 'H29410', 'H22331', 'H22368', 'H21729', 'H29556', 'H21956', 'H22140', 'H23309', 'H22101', 'H23157', 'H21593', 'H21990', 'H22228', 'H23028', 'H21915']
groups_subjects['Paired Initial Control'] = ['H27852', 'H28029', 'H26966', 'H27126', 'H29161', 'H28955', 'H26862', 'H27842', 'H27999', 'H28325', 'H26841', 'H27719', 'H27100', 'H27682', 'H29002', 'H27488', 'H27841', 'H28820', 'H28208', 'H27686']
groups_subjects['Paired Initial AMD'] = ['H29020', 'H26637', 'H27111', 'H26765', 'H28308', 'H28433', 'H26660', 'H28182', 'H27111', 'H27391', 'H28748', 'H28662', 'H26578', 'H28698', 'H27495', 'H28861', 'H28115', 'H28377', 'H26890', 'H28373', 'H27164']
groups_subjects['testing'] = ['H22825']
groups = ['Paired 2-YR AMD','Initial AMD','Initial Control','Paired 2-YR Control','Paired Initial Control','Paired Initial AMD']
#group_toview = 'Paired 2-YR AMD'
#group_toview = 'testing'
#roups = [group_toview]
"""
1 Cerebellum-Cortex_Right---Cerebellum-Cortex_Left 9 1 with weight of 3981.9602
2 rostralmiddlefrontal_Right---rostralmiddlefrontal_Left 76 42 with weight of 3781.9696
3 rostralmiddlefrontal_Right---middletemporal_Right 76 64 with weight of 3629.5854
4 superiorfrontal_Right---Cerebellum-Cortex_Right 77 9 with weight of 2834.8341
"""

ratio = 100
ratio_str = ratio_to_str(ratio)
str_identifier = '_MDT'+ratio_str
labeltype = 'lrordered'

verbose=True
picklesave=True


#1 lingual_Left---Cerebellum-Cortex_Right 28 9 with weight of 1194.3023
#2 lingual_Right---Cerebellum-Cortex_Left 62 1 with weight of 1070.304
#3 lingual_Left---Cerebellum-Cortex_Left 28 1 with weight of 986.7037
#4 lingual_Right---Cerebellum-Cortex_Right 62 9 with weight of 980.1034
#5 fusiform_Left---Cerebellum-Cortex_Right 22 9 with weight of 978.3528
target_tuple = (28, 9)
#target_tuple = (30, 50)

group_qb = {}
group_clusters = {}

for group in groups:
    groupstreamlines[group]=[]
    for ref in references:
        groupLines[group, ref]=[]
        groupPoints[group, ref]=[]
for group in groups:

    group_str = group.replace(' ', '_')
    _, _, index_to_struct, _ = atlas_converter(ROI_legends)
    centroid_file_path = os.path.join(centroid_folder, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]] + '.py')

    if not os.path.exists(centroid_file_path):

        subjects = groups_subjects[group]
        for subject in subjects:
            trkpath, exists = gettrkpath(TRK_folder, subject, str_identifier, pruned=False, verbose=verbose)
            if not exists:
                txt = f'Could not find subject {subject} at {TRK_folder} with {str_identifier}'
                warnings.warn(txt)
                continue
            #streamlines, header, _ = unload_trk(trkpath)
            trkdata = load_trk(trkpath, 'same')
            header = trkdata.space_attributes
            picklepath_grouping = os.path.join(pickle_folder, subject + str_identifier + '_grouping.p')
            grouping_xlsxpath = os.path.join(pickle_folder, subject + str_identifier + "_grouping.xlsx")
            if os.path.exists(picklepath_grouping):
                with open(picklepath_grouping, 'rb') as f:
                    grouping = pickle.load(f)
            elif os.path.exists(grouping_xlsxpath):
                grouping = get_grouping('grouping.xlsx')
            else:
                #affine_streams = np.eye(4)
                labelmask, labelaffine, labeloutpath, index_to_struct = getlabeltypemask(label_folder, 'MDT', ROI_legends, labeltype = labeltype, verbose = verbose)
                streamlines_2 = transform_streamlines(trkdata.streamlines, np.linalg.inv(trkdata.space_attributes[0]))
                M, grouping = connectivity_matrix(trkdata.streamlines, trkdata.space_attributes[0], labelmask, inclusive=True, symmetric=True,
                                    return_mapping=True,
                                    mapping_as_streamlines=False)
                if picklesave:
                    pickle.dump(grouping, open(picklepath_grouping, "wb"))
                    if verbose:
                        txt = ("The connectomes were saved at " + picklepath_grouping)
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
            target_qb = QuickBundles(threshold=distance1, metric=metric1)
            target_clusters = target_qb.cluster(target_streamlines_set)
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
        pickle.dump(group_clusters[group], open(centroid_file_path, "wb"))

    else:
        with open(centroid_file_path, 'rb') as f:
            group_clusters[group] = pickle.load(f)
        # registration
        #srr = StreamlineLinearRegistration()
        ##srm = srr.optimize(static=target_clusters_control.centroids, moving=target_clusters.centroids)
        #target_str_aligned = srm.transform(target_streamlines)
        #native_target_stream_aligned = transform_streamlines(target_str_aligned, np.linalg.inv(affine_fa))

from tract_save import save_trk_heavy_duty
"""
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

# groups = ['Paired 2-YR AMD','Initial AMD','Initial Control','Paired 2-YR Control','Paired Initial Control','Paired Initial AMD']
group_toview = 'Paired 2-YR AMD'
group_toview = 'Initial AMD'
viz_top_bundle = True
ref = None

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

    show_bundles(bundle_list, colors, fa=1)