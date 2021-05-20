from dipy.io.image import load_nifti, save_nifti
from dipy.viz import window, actor, ui
import numpy as np

labels_volume = "/Volumes/Data/Badea/Lab/mouse/VBM_19BrainChAMD01_IITmean_RPI_with_2yr-results/connectomics//H29056/H29056_IITmean_RPI_labels.nii.gz"
labels = load_nifti(labels_volume)

if np.size(np.shape(labels)) == 1:
    labels = labels[0]
if np.size(np.shape(labels)) == 4:
    labels = labels[:, :, :, 0]
print("Mask shape is " + str(np.shape(labels)))

affine = np.eye(4)
value_range = (0, 1)
scene = window.Scene()
scene.background((0.5, 0.5, 0.5))
slice_actor = actor.slicer(labels, affine, value_range)
scene.add(slice_actor)
slice_actor2 = slice_actor.copy()

slice_actor2.display(slice_actor2.shape[0]//2, None, None)
scene.add(slice_actor2)

scene.reset_camera()
scene.zoom(1.4)
window.record(scene, out_path='/Users/alex/jacques/slices4.png', size=(600, 600),
              reset_camera=False)
window.show(scene, size=(600, 600), reset_camera=False)
