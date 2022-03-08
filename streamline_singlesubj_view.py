import os
from dipy.viz import window, actor
from streamline_nocheck import load_trk as load_trk_spe
from dipy.io.streamline import load_trk, save_trk
from tract_visualize import show_bundles, setup_view

trk_folder = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_MDT/'
trk_name = 'S03350_stepsize_2_all_wholebrain_pruned.trk'
figures_path = '/Users/alex/jacques/AD_decode_abstract/figures/'
anat_path = '/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_Control_n72_i6/median_images/MDT_dwi.nii.gz'

trk_file_path = os.path.join(trk_folder,trk_name)

if os.path.exists(trk_file_path):
    try:
        streamlines_data = load_trk(trk_file_path, 'same')
    except:
        streamlines_data = load_trk_spe(trk_file_path, 'same')
else:
    raise Exception('cannot find file')

streamlines = streamlines_data.streamlines

hue = (0.5, 0.5)  # blue only, probably should change soonish
saturation = (0.0, 1.0)  # black to white
lut_cmap = actor.colormap_lookup_table(
    scale_range=(0, 1),
    hue_range=hue,
    saturation_range=saturation)

# lut_cmap = actor.colormap_lookup_table(
#   scale_range=(0.01, 0.55))
lut_cmap = actor.colormap_lookup_table(
    scale_range=(0.1, 0.25))

record_path = os.path.join(figures_path, trk_name[0:6]+'_figure.png')

scene = None
interactive = True
# record_path = None
scene = setup_view(streamlines[:], colors = lut_cmap ,ref = anat_path, world_coords = True, objectvals = [None], colorbar=True, record = record_path, scene = scene, interactive = interactive)
interactive = False