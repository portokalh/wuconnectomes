
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
from dipy.tracking import utils
from dipy.viz import window, actor
from time import sleep
from dipy.data import two_cingulum_bundles
from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
from dipy.tracking.streamline import transform_streamlines


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
        fa, affine_fa = load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfaN54900.nii.gz')
        fa_actor = actor.slicer(fa, affine_fa)
        ren.add(fa_actor)

    if show:
        window.show(ren)
    if fname is not None:
        sleep(1)
        window.record(ren, n_frames=1, out_path=fname, size=(900, 900))

