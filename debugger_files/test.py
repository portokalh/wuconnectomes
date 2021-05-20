# %%
import numpy as np
import pandas as pd
import pickle
import copy
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
from dipy.tracking import utils
from dipy.viz import window, actor
from time import sleep
from dipy.data import two_cingulum_bundles
from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines

# %%
def show_bundles(bundles, colors=None, show=True, fname=None,fa = False, str_tube = False):

    ren = window.Renderer()
    ren.SetBackground(1., 1, 1)
    if str_tube:
        bundle_actor = actor.streamtube(bundles, colors, linewidth=0.5)
        ren.add(bundle_actor)
    else:
        for (i, bundle) in enumerate(bundles):
            color = colors[i]
    #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05

            lines_actor = actor.line(bundle, color,linewidth=2.5)
            #lines_actor.RotateX(-90)
            #lines_actor.RotateZ(90)
            ren.add(lines_actor)
        
    if fa:
        fa, affine_fa= load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfaN54900.nii.gz')
        fa_actor = actor.slicer(fa, affine_fa)
        ren.add(fa_actor)
    
    if show:
        window.show(ren)
    if fname is not None:
        sleep(1)
        window.record(ren, n_frames=1, out_path=fname, size=(900, 900))

# %%
#indicate the ROIs interested, note target_l<target_r
target_l = 51
target_r = 257-166+1000

# %%
#general information
l = ['N54717','N54718','N54719','N54720','N54722','N54759','N54760','N54761','N54762','N54763','N54764','N54765','N54766','N54770','N54771','N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806','N54807','N54818','N54824','N54825','N54826','N54837','N54838','N54843','N54844','N54856','N54857','N54858','N54859','N54860','N54861','N54873','N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892','N54893','N54897','N54898','N54899','N54900','N54915','N54916','N54917']
# gen4idx = [1,2,3,4,7,8,9,10,12,13,52,53,54]
# gen3idx = [14,15,16,17,18,19,30,21,22,23]
# gen4 = [l[i] for i in gen4idx]
# gen3 = [l[j] for j in gen3idx]
# gen = [gen3, gen4]

#exclude N54900
gen0idx = [5,6,11,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]
gen0 = [l[k] for k in gen0idx]
# oldidx = [5,6,11,24,25,26,27,28,29,30,31,32,33,34,45,46,47,48,49,50]
# old = [l[o] for o in oldidx]
oldselectidx = [6, 25, 27, 28, 29, 32, 33, 45, 46, 50]
old = [l[o] for o in oldselectidx]
youngidx = [35,36,37,38,39,40,41,42,43,44]
young = [l[y] for y in youngidx]
age = [young,old]


# %%
# import random
# random.seed(123)
# oldidx = random.sample(oldidx,10)
# old = [l[o] for o in oldidx]
# age = [young,old]

# %%
#set path
mypath = '/Users/alex/code/Wenlin/data'
outpath = '/Users/alex/code/Wenlin/Tracts_Registration/results'

# %%
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

# %%
#load the control animal
streams_control,hdr_control = load_trk(mypath+'/wenlin_results/N54900_bmCSA_detr_small.trk')
labels_control, affine_labels_control = load_nifti(mypath+'/wenlin_data/labels/fa_labels_warp_N54900_RAS.nii.gz') 
fa_control, affine_fa_control= load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfaN54900.nii.gz')

# %%
#uncomment this cell to generate both left and right
# labels_pair_control = copy.copy(labels_control)
# nonz_control = np.nonzero(labels_pair_control)
# for i in range(len(nonz_control[0])):
#     if labels_pair_control[nonz_control[0][i], nonz_control[1][i], nonz_control[2][i]]>=1000:
#         labels_pair_control[nonz_control[0][i], nonz_control[1][i], nonz_control[2][i]] -= 1000
# print('pair labels generated')

# #pair labels target control animals
# streams_fix_control = lambda : (sl for sl in streams_control if len(sl)>1)
# streamlines_control = Streamlines(streams_fix_control())
# M_control, grouping_control = utils.connectivity_matrix(streamlines_control, labels_pair_control, 
#                                                         affine=affine_labels_control, return_mapping=True,
#                                                         mapping_as_streamlines=True)

# target_streamlines_control = grouping_control[target_l, target_r]

# %%
#target control animals
streams_fix_control = lambda : (sl for sl in streams_control if len(sl)>1)
streamlines_control = Streamlines(streams_fix_control())
M_control, grouping_control = utils.connectivity_matrix(streamlines_control, labels_control, 
                                                        affine=affine_labels_control, return_mapping=True,
                                                        mapping_as_streamlines=True)

target_streamlines_control = grouping_control[target_l, target_r]

# %%
#cluster control animals
target_qb_control = QuickBundles(threshold=distance1,metric=metric1)
target_clusters_control = target_qb_control.cluster(target_streamlines_control)
print("Control Nb. clusters:", len(target_clusters_control))

# %%
#group calculation
for k in range(2):
    agegroup = age[k]
    createVar = locals()
    createVar['groupstreamlines'+str(k+1)] = [] #create empty list, 2-old, 1-young
    createVar['groupLinesFA'+str(k+1)] = [] #create empty list, 2-old, 1-young
    createVar['groupPointsFA'+str(k+1)] = [] #create empty list, 2-old, 1-young
    createVar['Nativegroupstreamlines'+str(k+1)] = [] #create empty list, 2-old, 1-young
    animallist = []
    print('Group'+str(k+1)+' started')
    for j in range(len(agegroup)):
        runno = agegroup[j]
        animallist.append(runno)
        streams,hdr = load_trk(mypath+'/wenlin_results/'+runno+'_bmCSA_detr_small.trk')
        labels, affine_labels = load_nifti(mypath+'/wenlin_data/labels/fa_labels_warp_'+runno+'_RAS.nii.gz') 
        fa, affine_fa= load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfa'+runno+'.nii.gz')
#         labels_pair = copy.copy(labels)
#         nonz = np.nonzero(labels_pair)
#         for i in range(len(nonz[0])):
#             if labels_pair[nonz[0][i], nonz[1][i], nonz[2][i]]>=1000:
#                 labels_pair[nonz[0][i], nonz[1][i], nonz[2][i]] -= 1000
#         #print('pair labels generated')
        
        #target moving animals
        streams_fix = lambda : (sl for sl in streams if len(sl)>1)
        streamlines = Streamlines(streams_fix())
        M, grouping = utils.connectivity_matrix(streamlines, labels, affine=affine_labels, 
                                                return_mapping=True,mapping_as_streamlines=True)

        target_streamlines_ = grouping[target_l, target_r]

        target_streamlines = set_number_of_points(target_streamlines_, nb_points=num_points2)
        
        #first clustering for transform matrix
        target_qb = QuickBundles(threshold=distance1,metric=metric1)
        target_clusters = target_qb.cluster(target_streamlines)
#         print('NO.'+str(j+1)+' '+runno+" Nb. clusters:", len(target_clusters))
        
        #attach fa information
        native_target_streamlines = transform_streamlines(target_streamlines, np.linalg.inv(affine_fa))
        stream_fa = []
        stream_point_fa = []
        for s in range(len(native_target_streamlines)):
            point_fa = [fa[int(k[0]),int(k[1]),int(k[2])] for k in native_target_streamlines[s]]
            stream_point_fa.append(point_fa)
            stream_fa.append(np.mean(point_fa))
        
        #registration
        srr = StreamlineLinearRegistration()
        srm = srr.optimize(static=target_clusters_control.centroids, moving=target_clusters.centroids)
        target_str_aligned = srm.transform(target_streamlines)
        native_target_stream_aligned = transform_streamlines(target_str_aligned, np.linalg.inv(affine_fa))
        
        
        locals()['groupstreamlines'+str(k+1)].extend(target_str_aligned)
        locals()['Nativegroupstreamlines'+str(k+1)].extend(native_target_stream_aligned)
        locals()['groupLinesFA'+str(k+1)].extend(stream_fa)
        locals()['groupPointsFA'+str(k+1)].extend(stream_point_fa)
        
        print('NO.'+str(j+1)+' '+runno+" Nb. streamlines:", len(target_str_aligned))
        
    print('agetype-'+str(k+1)+' finished. '+'total number of clusters for group'+ str(k+1) 
          + ': {}'.format(len(locals()['groupstreamlines'+str(k+1)])))
    print('animal list: ', animallist)
    print(' ')
        
        

# %%
group1_qb = QuickBundles(threshold=distance2,metric=metric2)
group1_clusters = group1_qb.cluster(groupstreamlines1)
group2_qb = QuickBundles(threshold=distance2,metric=metric2)
group2_clusters = group2_qb.cluster(groupstreamlines2)
print("Young Group Nb. clusters:", len(group1_clusters))
print("Old Group Nb. clusters:", len(group2_clusters))

# %%
#color by line-average fa

group = 2

if group == 1:
    Nativegroupstreamlines = Nativegroupstreamlines1
    groupLinesFA = groupLinesFA1
    name = 'Group_Young'
else:
    Nativegroupstreamlines = Nativegroupstreamlines2
    groupLinesFA = groupLinesFA2
    name = 'Group_Old'
    
cmap = actor.colormap_lookup_table(
    scale_range=(np.min(groupLinesFA), np.max(groupLinesFA)))

renderer = window.Renderer()
stream_actor = actor.line(Nativegroupstreamlines,np.array(groupLinesFA),lookup_colormap=cmap)
fa_actor = actor.slicer(fa_control, np.eye(4))
renderer.add(stream_actor)
renderer.add(fa_actor)
bar = actor.scalar_bar(cmap)
renderer.add(bar)
# Uncomment the line below to show to display the window
window.show(renderer, size=(600, 600), reset_camera=False)
window.record(renderer,size=(600, 600),
              out_path = outpath+'/'+str(target_l)+'--'+str(target_r)+name+' lineFA Viz.png')

# %%
##### color by points
group = 1

if group == 1:
    Nativegroupstreamlines = Nativegroupstreamlines1
    groupPointsFA = groupPointsFA1
    name = 'Group_Young'
else:
    Nativegroupstreamlines = Nativegroupstreamlines2
    groupPointsFA = groupPointsFA2
    name = 'Group_Old'

cmap = actor.colormap_lookup_table(
scale_range=(np.min(groupPointsFA), np.max(groupPointsFA)))

renderer = window.Renderer()
renderer.clear()
renderer = window.Renderer()
colors = []
for s in range(len(Nativegroupstreamlines)):
    stream = Nativegroupstreamlines[s]
    for idx in range(len(stream)):
        colors.append(groupPointsFA[s][idx])

stream_actor2 = actor.line(Nativegroupstreamlines, colors, linewidth=0.2,lookup_colormap=cmap)

renderer.add(stream_actor2)
fa_actor2 = actor.slicer(fa_control, np.eye(4))
renderer.add(fa_actor2)
bar = actor.scalar_bar(cmap)
renderer.add(bar)

# Uncomment the line below to show to display the window
window.show(renderer, size=(600, 600), reset_camera=False)
window.record(renderer,size=(600, 600),
              out_path = outpath+'/'+str(target_l)+'--'+str(target_r)+name+' PointFA Viz.png')

# %%
#viz a specific bundle with line FA 

group = 1
bundle_id = 40

if group == 1:
    group_clusters = group1_clusters
    groupLinesFA = groupLinesFA1
    name = 'Group_Young-Bundle '+str(bundle_id)
else:
    group_clusters = group2_clusters
    groupLinesFA = groupLinesFA2
    name = 'Group_Old-Bundle '+str(bundle_id)

k = group_clusters.clusters[bundle_id]
bundle_fa = []
for idx in k.indices:
    bundle_fa.append(groupLinesFA[idx])

#         cmap = actor.colormap_lookup_table(
#         scale_range=(np.min(bundle_fa), np.max(bundle_fa)))
cmap = actor.colormap_lookup_table(
scale_range=(0.1, 0.5))

#color by line-average fa
renderer = window.Renderer()
renderer.clear()
renderer = window.Renderer()
stream_actor3 = actor.line(group_clusters.clusters[bundle_id],np.array(bundle_fa),lookup_colormap=cmap)
renderer.add(stream_actor3)
bar = actor.scalar_bar(cmap)
renderer.add(bar)
# Uncomment the line below to show to display the window
window.show(renderer, size=(600, 600), reset_camera=False)
#window.record(renderer,size=(600, 600),
              #out_path = outpath+'/'+str(target_l)+'--'+str(target_r)+name+' lineFA Viz.png')

# %%
#viz top bundle
#swtich fa=0 to not show fa slice, switch to 1 to show fa slice
np.random.seed(123)
group = 1
num_of_bundles =6

if group == 1:
    group_clusters = group1_clusters
    groupLinesFA = groupLinesFA1
    name = 'Group_Young-Bundle top '+str(num_of_bundles)
else:
    group_clusters = group2_clusters
    groupLinesFA = groupLinesFA2
    name = 'Group_Old-Bundle '+str(num_of_bundles)
    
top_idx = sorted(range(len(group_clusters.clusters_sizes())), key=lambda i: group_clusters.clusters_sizes()[i],
             reverse=True)[:num_of_bundles]

bundle_list = [group_clusters.clusters[idx] for idx in top_idx]
color_list = [(np.random.randint(0, 255),np.random.randint(0, 255),np.random.randint(0, 255)) 
              for n in range(num_of_bundles)]
color_list_dis_all = [window.colors.green,window.colors.yellow,
                  window.colors.red,window.colors.brown,
                  window.colors.orange,window.colors.blue]
color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]
if num_of_bundles <= 6:
    colors = color_list_dis
else:
    colors = color_list
    
show_bundles(bundle_list,colors,fa = 0)

# %%
#viz top bundle with centroid
np.random.seed(123)
group = 2
num_of_bundles =4

if group == 1:
    group_clusters = group1_clusters
    groupLinesFA = groupLinesFA1
    name = 'Group_Young-Bundle top '+str(num_of_bundles)
else:
    group_clusters = group2_clusters
    groupLinesFA = groupLinesFA2
    name = 'Group_Old-Bundle '+str(num_of_bundles)
    
top_idx = sorted(range(len(group_clusters.clusters_sizes())), key=lambda i: group_clusters.clusters_sizes()[i],
             reverse=True)[:num_of_bundles]

bundle_list = [group_clusters.centroids[idx] for idx in top_idx]
color_list = [(np.random.randint(0, 255),np.random.randint(0, 255),np.random.randint(0, 255)) 
              for n in range(num_of_bundles)]
color_list_dis_all = [window.colors.green,window.colors.yellow,
                  window.colors.red,window.colors.brown,
                  window.colors.orange,window.colors.blue]
color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]

if num_of_bundles <= 6:
    colors = color_list_dis
else:
    colors = color_list
    
show_bundles(bundle_list,colors,str_tube=True,fa = 0)

# %%
#viz top bundle with centroid and bundles
np.random.seed(123)
group = 2
num_of_bundles =6

if group == 1:
    group_clusters = group1_clusters
    groupLinesFA = groupLinesFA1
    name = 'Group_Young-Bundle top '+str(num_of_bundles)
else:
    group_clusters = group2_clusters
    groupLinesFA = groupLinesFA2
    name = 'Group_Old-Bundle '+str(num_of_bundles)

top_idx = sorted(range(len(group_clusters.clusters_sizes())), key=lambda i: group_clusters.clusters_sizes()[i],
             reverse=True)[:num_of_bundles]

color_list = [(np.random.randint(0, 255),np.random.randint(0, 255),np.random.randint(0, 255)) 
              for n in range(num_of_bundles)]
color_list_dis_all = [window.colors.green,window.colors.yellow,
                  window.colors.red,window.colors.brown,
                  window.colors.orange,window.colors.blue]
color_list_dis = [color_list_dis_all[i] for i in range(num_of_bundles)]

if num_of_bundles <= 6:
    colors = color_list_dis
else:
    colors = color_list

fa, affine_fa= load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfaN54900.nii.gz')
fa_actor = actor.slicer(fa, affine_fa)

        
bundle_list_center = [group_clusters.centroids[idx] for idx in top_idx]
bundle_list = [group_clusters.clusters[idx] for idx in top_idx]
ren = window.Renderer()
ren.SetBackground(1., 1, 1)

bundle_actor = actor.streamtube(bundle_list_center, colors, linewidth=0.1)


for (i, bundle) in enumerate(bundle_list):
    color = colors[i]
    lines_actor = actor.line(bundle, color,opacity=0.9)
    ren.add(lines_actor)
    
ren.add(bundle_actor) 
#uncomment this to show fa
#ren.add(fa_actor)
window.show(ren)        

# %%
num_of_bundles = 20
top_idx_group1 = sorted(range(len(group1_clusters.clusters_sizes())),
                        key=lambda i: group1_clusters.clusters_sizes()[i],reverse=True)[:num_of_bundles]
top_idx_group2 = sorted(range(len(group2_clusters.clusters_sizes())),
                        key=lambda i: group2_clusters.clusters_sizes()[i],reverse=True)[:num_of_bundles]

# %%
bundle_group1 = [group1_clusters.centroids[idx] for idx in top_idx_group1]
bundle_group2 = [group2_clusters.centroids[idx] for idx in top_idx_group2]

# %%
from dipy.segment.metric import ResampleFeature, AveragePointwiseEuclideanMetric,mdf
dist_all = np.zeros((num_of_bundles,num_of_bundles))
for g1 in range(len(bundle_group1)):
    for g2 in range(len(bundle_group2)):
        id1 = top_idx_group1[g1]
        id2 = top_idx_group2[g2]
        dist_all[g1,g2] = (mdf(group1_clusters.centroids[id1],group2_clusters.centroids[id2]))

# %%
import copy
dist_all_fix = copy.copy(dist_all)
dist_all_idx = []
for i in range(len(bundle_group1)):
        idx = np.argmin(dist_all_fix[i,:])
        dist_all_idx.append([i,idx])
        dist_all_fix[:,idx] = 100000
#dist_all_idx

dist_group1_idx = [dist_all_idx[iii][0] for iii in range(num_of_bundles)]#size id
dist_group2_idx = [dist_all_idx[iii][1] for iii in range(num_of_bundles)]#size id

correspond_bundle_id = [] #bundle id
for i in range(6):
    correspond_bundle_id.append([top_idx_group1[dist_all_idx[i][0]],top_idx_group2[dist_all_idx[i][1]]])
    print(str(top_idx_group1[dist_all_idx[i][0]])+'--'+str(top_idx_group2[dist_all_idx[i][1]]))

# %%
group1List = [top_idx_group1[dist_all_idx[i][0]] for i in range(6)]
group2List = [top_idx_group2[dist_all_idx[i][1]] for i in range(6)]
print(group1List)
print(group2List)
print(correspond_bundle_id)

# %%
#viz bundle for loop
num_of_top_bundle = 6
group1List = [top_idx_group1[dist_all_idx[i][0]] for i in range(6)]
group2List = [top_idx_group2[dist_all_idx[i][1]] for i in range(6)]
groupList = [[],group1List,group2List]
for groupidx in [1,2]:
    for num in groupList[groupidx]:   #[0,2,6,3,8,1] for group2
        cor_id = groupList[groupidx].index(num)
        group = groupidx
        bundle_id = num

        if group == 1:
            group_clusters = group1_clusters
            groupLinesFA = groupLinesFA1
            groupPointsFA = groupPointsFA1
            name = 'Group_Young-Bundle_'+str(bundle_id)+'_corID'+str(cor_id)
        if group == 2:
            group_clusters = group2_clusters
            groupLinesFA = groupLinesFA2
            groupPointsFA = groupPointsFA2
            name = 'Group_Old-Bundle_'+str(bundle_id)+'_corID'+str(cor_id)

        bundle_point_fa = []
        bundle_fa = []
        k = group_clusters.clusters[bundle_id]
        for idx in k.indices:
            bundle_fa.append(groupLinesFA[idx])
            for idx_point in range(len(groupPointsFA[idx])):
                bundle_point_fa.append(groupPointsFA[idx][idx_point])


#         cmap = actor.colormap_lookup_table(
#         scale_range=(np.min(bundle_fa), np.max(bundle_fa)))
        cmap = actor.colormap_lookup_table(
        scale_range=(0.1, 0.6))

        #color by line-average fa
        renderer = window.Renderer()
        renderer.clear()
        renderer = window.Renderer()
#         stream_actor3 = actor.line(group_clusters.clusters[bundle_id],np.array(bundle_fa),lookup_colormap=cmap)
        stream_actor3 = actor.line(group_clusters.clusters[bundle_id],
                                   np.array(bundle_point_fa),lookup_colormap=cmap)
        renderer.add(stream_actor3)
        bar = actor.scalar_bar(cmap)
        renderer.add(bar)
        # Uncomment the line below to show to display the window
        #window.show(renderer, size=(600, 600), reset_camera=False)
        window.record(renderer,size=(600, 600),
                      out_path = outpath+'/'+str(target_l)+'--'+str(target_r)+name+' PointFA Viz.png')

# %%
group1csv = np.zeros((1,6))
for i in range(6):
    idsize = dist_group1_idx[i]
    idbundle = group1List[i]
    fa = []
    for s in group1_clusters.clusters[idbundle].indices:
            temp = np.hstack((idsize*np.ones((num_points2,1)),
                              idbundle*np.ones((num_points2,1)),
                              s*np.ones((num_points2,1)),
                              np.array(range(num_points2)).reshape(num_points2,1),
                              np.array(groupPointsFA1[s]).reshape(num_points2,1),
                             list(utils.length([groupstreamlines1[s]]))*np.ones((num_points2,1))))
            group1csv = np.vstack((group1csv,temp))
group1csv = group1csv[1:,:]
group1csvDF = pd.DataFrame(group1csv)
group1csvDF.rename(index=str, columns={0:"Bundle Size Rank",1:"Bundle ID",2:"Steamlines ID",
                                       3:"Point ID", 4:"FA", 5:"length"})
group1csvDF.to_csv(outpath+'/'+str(target_l)+'--'+str(target_r)+'group1FA.csv'
                   ,header = ["Bundle Size Rank","Bundle ID","Streamlines ID",
                                                     "Point ID", "FA","Length"])

# %%
group2csv = np.zeros((1,6))
for i in range(6):
    idsize = dist_group2_idx[i]
    idbundle = group2List[i]
    fa = []
    for s in group2_clusters.clusters[idbundle].indices:
            temp = np.hstack((idsize*np.ones((num_points2,1)),
                              idbundle*np.ones((num_points2,1)),
                              s*np.ones((num_points2,1)),
                              np.array(range(num_points2)).reshape(num_points2,1),
                              np.array(groupPointsFA2[s]).reshape(num_points2,1),
                             list(utils.length([groupstreamlines2[s]]))*np.ones((num_points2,1))))
            group2csv = np.vstack((group2csv,temp))
group2csv = group2csv[1:,:]  
group2csvDF = pd.DataFrame(group2csv)
group2csvDF.rename(index=str, columns={0:"Bundle Size Rank",1:"Bundle ID",2:"Steamlines ID",
                                       3:"Point ID", 4:"FA", 5:"length"})
group2csvDF.to_csv(outpath+'/'+str(target_l)+'--'+str(target_r)+'group2FA.csv'
                   ,header = ["Bundle Size Rank","Bundle ID","Streamlines ID",
                                                     "Point ID", "FA","Length"])

# %%
#before alignment - centroids for moving and reference
show_bundles([target_clusters_control.centroids,target_clusters.centroids],
             colors = [window.colors.orange,window.colors.green])

# %%
#after alignment - centroids for moving and reference
show_bundles([target_clusters_control.centroids,target_str_aligned],
             colors = [window.colors.orange,window.colors.green])

# %%
