import numpy as np
from scipy.ndimage.morphology import binary_dilation

from dipy.data import read_stanford_labels, read_stanford_t1
from dipy.direction import peaks
from dipy.reconst import shm
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
from dipy.tracking.streamline import Streamlines

outpath = BIGGUS_DISKUS + "/C57_JS/Testzone/"

hardi_img, gtab, labels_img = read_stanford_labels()
data = hardi_img.get_data()
labels = labels_img.get_data()

t1 = read_stanford_t1()
t1_data = t1.get_data()

white_matter = binary_dilation((labels == 1) | (labels == 2))
csamodel = shm.CsaOdfModel(gtab, 6)
csapeaks = peaks.peaks_from_model(model=csamodel,
                                  data=data,
                                  sphere=peaks.default_sphere,
                                  relative_peak_threshold=.8,
                                  min_separation_angle=45,
                                  mask=white_matter)

affine = np.eye(4)
seeds = utils.seeds_from_mask(white_matter, affine, density=1)
stopping_criterion = BinaryStoppingCriterion(white_matter)

streamline_generator = LocalTracking(csapeaks, stopping_criterion, seeds,
                                     affine=affine, step_size=0.5)
streamlines = Streamlines(streamline_generator)

cc_slice = labels == 2
cc_streamlines = utils.target(streamlines, affine, cc_slice)
cc_streamlines = Streamlines(cc_streamlines)

other_streamlines = utils.target(streamlines, affine, cc_slice,
                                 include=False)
other_streamlines = Streamlines(other_streamlines)
assert len(other_streamlines) + len(cc_streamlines) == len(streamlines)

from dipy.viz import window, actor, colormap as cmap

# Enables/disables interactive visualization
interactive = False

# Make display objects
color = cmap.line_colors(cc_streamlines)
cc_streamlines_actor = actor.line(cc_streamlines,
                                  cmap.line_colors(cc_streamlines))
cc_ROI_actor = actor.contour_from_roi(cc_slice, color=(1., 1., 0.),
                                      opacity=0.5)

vol_actor = actor.slicer(t1_data)

vol_actor.display(x=40)
vol_actor2 = vol_actor.copy()
vol_actor2.display(z=35)

# Add display objects to canvas
r = window.Renderer()
r.add(vol_actor)
r.add(vol_actor2)
r.add(cc_streamlines_actor)
r.add(cc_ROI_actor)

# Save figures
window.record(r, n_frames=1, out_path=outpath+'corpuscallosum_axial.png',
              size=(800, 800))
if interactive:
    window.show(r)
r.set_camera(position=[-1, 0, 0], focal_point=[0, 0, 0], view_up=[0, 0, 1])
window.record(r, n_frames=1, out_path=outpath+'corpuscallosum_sagittal.png',
              size=(800, 800))
if interactive:
    window.show(r)

M, grouping = utils.connectivity_matrix(cc_streamlines, affine, labels,
                                        return_mapping=True,
                                        mapping_as_streamlines=True)
M[:3, :] = 0
M[:, :3] = 0

import numpy as np
import matplotlib.pyplot as plt
plt.imshow(np.log1p(M), interpolation='nearest')
plt.savefig(outpath+"connectivity.png")

lr_superiorfrontal_track = grouping[11, 54]
shape = labels.shape
dm = utils.density_map(lr_superiorfrontal_track, affine, shape)

import nibabel as nib
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import save_trk

# Save density map
dm_img = nib.Nifti1Image(dm.astype("int16"), hardi_img.affine)
dm_img.to_filename(outpath+"lr-superiorfrontal-dm.nii.gz")

lr_sf_trk = Streamlines(lr_superiorfrontal_track)

# Save streamlines
sft = StatefulTractogram(lr_sf_trk, dm_img, Space.VOX)
save_trk(sft, "lr-superiorfrontal.trk")