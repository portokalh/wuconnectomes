from dipy.io.streamline import load_trk, save_trk
from dipy.viz import window, actor
import os
import pickle
from tract_visualize import show_bundles, setup_view
from convert_atlas_mask import convert_labelmask, atlas_converter
from tract_handler import ratio_to_str, gettrkpath
from itertools import compress
import numpy as np
import nibabel as nib


figures_path = '/Volumes/Data/Badea/Lab/human/AMD/Figures_MDT_non_inclusive/'

record = ''

target_tuple = (9,1)
#target_tuple = (9,77)
ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
_, _, index_to_struct, _ = atlas_converter(ROI_legends)

centroid_folder = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/'
groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
          'Paired Initial AMD']

anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'

ratio = 1
ratio_str = ratio_to_str(ratio)

text_path = os.path.join(figures_path,index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]]+'_stats.txt')

testfile = open(text_path, "w")
testfile.write("Parameters for groups\n")
testfile.close()

top_percentile = 2

#superior frontal right to cerebellum right

scene = None
print(index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]])

for group in groups:

    print(f'Setting up group {group}')
    group_str = group.replace(' ', '_')

    centroid_file_path = os.path.join(centroid_folder,
                                      group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                      index_to_struct[target_tuple[1]] + '_centroid.py')
    fa_path = os.path.join(centroid_folder,
                                      group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                      index_to_struct[target_tuple[1]] + '_fa_lines.py')
    md_path = os.path.join(centroid_folder,
                                      group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                      index_to_struct[target_tuple[1]] + '_md_lines.py')
    trk_path = os.path.join(centroid_folder,
                                      group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                      index_to_struct[target_tuple[1]] + '_streamlines.trk')

    if os.path.exists(fa_path):
        with open(fa_path, 'rb') as f:
            fa_lines = pickle.load(f)
    if os.path.exists(md_path):
        with open(md_path, 'rb') as f:
            md_lines = pickle.load(f)

    streamlines_data = load_trk(trk_path,'same')
    streamlines = streamlines_data.streamlines

    testfile = open(text_path, "a")
    testfile.write(f"Mean, Median, max and std FA for streamlines in group {group}: \n{np.mean(fa_lines)}, "
                   f"{np.median(fa_lines)}, {np.max(fa_lines)}, {np.std(fa_lines)}\n")
    testfile.write(f"Mean, Median, max and std MD for streamlines in group {group}: \n{np.mean(md_lines)}, "
                   f"{np.median(md_lines)}, {np.max(md_lines)}, {np.std(md_lines)}\n")
    testfile.write(f"Number of streamlines in group {group}: \n{np.shape(streamlines)[0]}\n")
    testfile.close()

    cutoff = np.percentile(fa_lines,100 - top_percentile)
    select_streams = fa_lines>cutoff
    fa_lines = list(compress(fa_lines, select_streams))
    streamlines = list(compress(streamlines, select_streams))
    streamlines = nib.streamlines.ArraySequence(streamlines)

    """
    scene = window.Scene()
    
    scene.clear()
    stream_actor4 = actor.line(streamlines, (1., 0.5, 0), linewidth=0.1)
    
    scene.add(stream_actor4)
    """
    hue = (0.5, 0.5)  # blue only, probably should change soonish
    saturation = (0.0, 1.0)  # black to white
    lut_cmap = actor.colormap_lookup_table(
        scale_range=(0, 1),
        hue_range=hue,
        saturation_range=saturation)

    lut_cmap = actor.colormap_lookup_table(
        scale_range=(0.45, 0.55))

    record_path = os.path.join(figures_path, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                      index_to_struct[target_tuple[1]] + '_figure.png')

    #scene = None
    #scene = setup_view(streamlines[:], colors = lut_cmap,ref = anat_path, world_coords = True, objectvals = fa_lines[:], colorbar=True, record = record_path, scene = scene)
    #add something to help make the camera static over multiple iterations? Would be VERY nice.


    ##write text file here that summarizes fa, md for each group, would be very helpful


"""
centroids_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive_100/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_centroid.py'
fa_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive_100/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_fa_lines.py'
md_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive_100/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_fa_lines.py'
trk_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive_100/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_streamlines.trk'


centroids_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_centroid.py'
fa_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_fa_lines.py'
md_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_fa_lines.py'
trk_path = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/Initial_AMD_MDT_ratio_100_right-cerebellum-cortex_right_to_left-cerebellum-cortex_left_streamlines.trk'
"""