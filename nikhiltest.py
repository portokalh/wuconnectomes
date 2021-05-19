#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:54:49 2020

@author: alex
"""

from dipy.viz import window, actor
from time import sleep
from dipy.data import two_cingulum_bundles

cb_subj1, cb_subj2 = two_cingulum_bundles()

from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
from dipy.io.streamline import load_trk

trkprunedata = load_trk("/Users/alex/jacques/Nikhiltest/N57442_small_stepsize_2_pruned.trk", "same")
affine = trkprunedata._affine
trkprunedata.to_vox()
pruned_streamlines_SL = trkprunedata.streamlines

streamlines_test = list(pruned_streamlines_SL)
#endpoints = [sl[0::len(sl) - 1] for sl in streamlines_test]
#lin_T, offset = _mapping_to_voxel(affine)
#endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)
#i, j, k = endpoints.T

cb_subj1 = set_number_of_points(cb_subj1, 20)
cb_subj2 = set_number_of_points(cb_subj2, 20)

srr = StreamlineLinearRegistration()

srm = srr.optimize(static=cb_subj1, moving=cb_subj2)

cb_subj2_aligned = srm.transform(cb_subj2)


def show_both_bundles(bundles, colors=None, show=True, fname=None):

    scene = window.Scene()
    scene.SetBackground(1., 1, 1)
    for (i, bundle) in enumerate(bundles):
        color = colors[i]
        lines_actor = actor.streamtube(bundle, color, linewidth=0.3)
        lines_actor.RotateX(-90)
        lines_actor.RotateZ(90)
        scene.add(lines_actor)
    if show:
        window.show(scene)
    if fname is not None:
        sleep(1)
        window.record(scene, n_frames=1, out_path=fname, size=(900, 900))


show_both_bundles([cb_subj1, cb_subj2],
                  colors=[window.colors.orange, window.colors.red],
                  show=True,
                  fname=None)

show_both_bundles([cb_subj1, cb_subj2_aligned],
                  colors=[window.colors.orange, window.colors.red],
                  show=True,
                  fname=None)

show_both_bundles([streamlines_test],
                  colors=[window.colors.orange, window.colors.red],
                  show=True,
                  fname=None)

print('hi')