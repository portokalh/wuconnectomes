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
import warnings

project = 'AD_Decode'

fixed = True
record = ''

computer_name = socket.gethostname()

inclusive = False
symmetric = True
write_txt = True
ratio = 1
top_percentile = 2
groups = ['APOE4', 'APOE3']
#groups = ['Male','Female']
#,(23,30)
target_tuples = [(9, 1), (24,1), (22, 1), (58, 57), (64, 57)]
target_tuples = [(9, 1), (24,1), (22, 1), (58, 57),  (23,24), (64, 57)]
target_tuples = [(58, 57), (9, 1), (24,1), (22, 1), (64, 57),(23,24),(24,30),(23,30)]
target_tuples = [(24,30),(23,24),(24,30)]
target_tuples = [(24,30)]


changewindow_eachtarget = False

if inclusive:
    inclusive_str = '_inclusive'
else:
    inclusive_str = '_non_inclusive'

if symmetric:
    symmetric_str = '_symmetric'
else:
    symmetric_str = '_non_symmetric'

#if fixed:
#    fixed_str = '_fixed'
#else:
#    fixed_str = ''

samos = False
if 'samos' in computer_name:
    mainpath = '/mnt/paros_MRI/jacques/'
    ROI_legends = "/mnt/paros_MRI/jacques/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'santorini' in computer_name:
    #mainpath = '/Users/alex/jacques/'
    mainpath = '/Volumes/Data/Badea/Lab/human/'
    ROI_legends = "/Volumes/Data/Badea/ADdecode.01/Analysis/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
elif 'blade' in computer_name:
    mainpath = '/mnt/munin6/Badea/Lab/human/'
    ROI_legends = "/mnt/munin6/Badea/Lab/atlases/IITmean_RPI/IITmean_RPI_index.xlsx"
else:
    raise Exception('No other computer name yet')

#target_tuple = (24,1)
#target_tuple = [(58, 57)]
#target_tuples = [(64, 57)]


ratio_str = ratio_to_str(ratio)
print(ratio_str)
if ratio_str == '_all':
    folder_ratio_str = ''
else:
    folder_ratio_str = ratio_str.replace('_ratio', '')
#target_tuple = (9,77)

_, _, index_to_struct, _ = atlas_converter(ROI_legends)

if project == 'AMD':
    mainpath = os.path.join(mainpath,project)
    groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
              'Paired Initial AMD']
    anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'

if project == 'AD_Decode':
    mainpath = os.path.join(mainpath, project, 'Analysis')
    anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_21ADDecode03_IITmean_RPI_fullrun-work/dwi/SyN_0p5_3_0p5_fa/faMDT_NoNameYet_n37_i6/median_images/MDT_b0.nii.gz'


#figures_path = '/Volumes/Data/Badea/Lab/human/AMD/Figures_MDT_non_inclusive/'
#centroid_folder = '/Volumes/Data/Badea/Lab/human/AMD/Centroids_MDT_non_inclusive/'
figures_path = os.path.join(mainpath, f'Figures_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
centroid_folder = os.path.join(mainpath, f'Centroids_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
trk_folder = os.path.join(mainpath, f'Centroids_MDT{inclusive_str}{symmetric_str}{folder_ratio_str}')
mkcdir([figures_path, centroid_folder])

#groups = ['Initial AMD', 'Paired 2-YR AMD', 'Initial Control', 'Paired 2-YR Control', 'Paired Initial Control',
#          'Paired Initial AMD']

#anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'



#superior frontal right to cerebellum right

scene = None

interactive = True

for target_tuple in target_tuples:


    print(target_tuple[0], target_tuple[1])
    print(index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[target_tuple[1]])

    if write_txt:
        text_path = os.path.join(figures_path, index_to_struct[target_tuple[0]] + '_to_' + index_to_struct[
            target_tuple[1]] + '_stats.txt')
        testfile = open(text_path, "w")
        testfile.write("Parameters for groups\n")
        testfile.close()

    if changewindow_eachtarget:
        firstrun = True

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
        trk_path = os.path.join(trk_folder,
                                          group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                          index_to_struct[target_tuple[1]] + '_streamlines.trk')

        if os.path.exists(fa_path):
            with open(fa_path, 'rb') as f:
                fa_lines = pickle.load(f)
        if os.path.exists(md_path):
            with open(md_path, 'rb') as f:
                md_lines = pickle.load(f)
        if os.path.exists(trk_path):
            try:
                streamlines_data = load_trk(trk_path,'same')
            except:
                streamlines_data = load_trk_spe(trk_path, 'same')
        streamlines = streamlines_data.streamlines

        if write_txt:
            testfile = open(text_path, "a")
            testfile.write(f"Mean, Median, max and std FA for streamlines in group {group}: \n{np.mean(fa_lines)}, "
                           f"{np.median(fa_lines)}, {np.max(fa_lines)}, {np.std(fa_lines)}\n")
            testfile.write(f"Mean, Median, max and std MD for streamlines in group {group}: \n{np.mean(md_lines)}, "
                           f"{np.median(md_lines)}, {np.max(md_lines)}, {np.std(md_lines)}\n")
            testfile.write(f"Number of streamlines in group {group}: \n{np.shape(streamlines)[0]}\n")
            testfile.close()

        if 'fa_lines' in locals():
            cutoff = np.percentile(fa_lines,100 - top_percentile)
            select_streams = fa_lines>cutoff
            fa_lines_new = list(compress(fa_lines, select_streams))
            streamlines_new = list(compress(streamlines, select_streams))
            streamlines_new = nib.streamlines.ArraySequence(streamlines_new)
            if np.shape(streamlines)[0] != np.shape(fa_lines)[0]:
                raise Exception('Inconsistency between streamlines and fa lines')
        else:
            txt = f'Cannot find {fa_path}, could not select streamlines based on fa'
            warnings.warn(txt)
            fa_lines = [None]

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

        #lut_cmap = actor.colormap_lookup_table(
         #   scale_range=(0.01, 0.55))
        lut_cmap = actor.colormap_lookup_table(
            scale_range=(0.1, 0.25))

        record_path = os.path.join(figures_path, group_str + '_MDT' + ratio_str + '_' + index_to_struct[target_tuple[0]] + '_to_' +
                                          index_to_struct[target_tuple[1]] + '_figure.png')

        #scene = None
        #interactive = True
        #record_path = None
        scene = setup_view(streamlines_new[:], colors = lut_cmap,ref = anat_path, world_coords = True, objectvals = fa_lines_new[:], colorbar=True, record = record_path, scene = scene, interactive = interactive)
        del(fa_lines,fa_lines_new,streamlines,streamlines_new)
        interactive = False