#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Small tools for sorting through tracts
"""


import nibabel as nib
import numpy as np
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from dipy.io.streamline import load_trk
from os.path import splitext
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
import pickle

from types import ModuleType, FunctionType
from gc import get_referents

import smtplib

import os, re, sys, io, struct, socket, datetime
from email.mime.text import MIMEText
import glob

from dipy.tracking.utils import unique_rows


from time import time
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.local_tracking import LocalTracking
from dipy.direction import peaks
from nibabel.streamlines import detect_format
from dipy.io.utils import (create_tractogram_header,
                           get_reference_info)
from dipy.viz import window, actor

from dipy.segment.mask import segment_from_cfa
from dipy.segment.mask import bounding_box

import multiprocessing
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
import multiprocessing.pool



from scipy.ndimage.morphology import binary_dilation
from dipy.tracking import utils
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
from dipy.tracking.streamline import Streamlines
import matplotlib.pyplot as plt

#from dipy.denoise.localpca import mppca
from denoise_processes import mppca
from dipy.denoise.gibbs import gibbs_removal

from random import randint

from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib
import matplotlib.pyplot as plt
#import dipy.tracking.life as life
import JSdipy.tracking.life as life
from dipy.viz import window, actor, colormap as cmap
import dipy.core.optimize as opt
from functools import wraps

from bvec_handler import fix_bvals_bvecs#, extractbvec_fromheader
from figures_handler import denoise_fig, show_bundles, window_show_test, LifEcreate_fig
from tract_eval import bundle_coherence, LiFEvaluation
from dif_to_trk import make_tensorfit, QCSA_tractmake
from BIAC_tools import MyPool

def longstring(string,margin=0):

    teststring = list(string)
    val1 = teststring.index(next(filter(lambda i: i != 0, teststring)))
    teststring.reverse()
    val2 = np.size(string) - teststring.index(next(filter(lambda i: i != 0, teststring)))
    if (val1 - margin) < 0:
        val1 = 0
    else:
        val1 = val1 - margin
    if (val2 + margin) > np.size(string):
        val2 = np.size(string)
    else:
        val2 = val2 + margin
    string[val1:val2] = True
    return string


def _with_initialize(generator):
    """Allows one to write a generator with initialization code.

    All code up to the first yield is run as soon as the generator function is
    called and the first yield value is ignored.
    """
    @wraps(generator)
    def helper(*args, **kwargs):
        gen = generator(*args, **kwargs)
        next(gen)
        return gen

    return helper

@_with_initialize
def target(streamlines, affine, target_mask, include=True, strict=False):
    """Filters streamlines based on whether or not they pass through an ROI.

    Parameters
    ----------
    streamlines : iterable
        A sequence of streamlines. Each streamline should be a (N, 3) array,
        where N is the length of the streamline.
    affine : array (4, 4)
        The mapping between voxel indices and the point space for seeds.
        The voxel_to_rasmm matrix, typically from a NIFTI file.
    target_mask : array-like
        A mask used as a target. Non-zero values are considered to be within
        the target region.
    include : bool, default True
        If True, streamlines passing through `target_mask` are kept. If False,
        the streamlines not passing through `target_mask` are kept.

    Returns
    -------
    streamlines : generator
        A sequence of streamlines that pass through `target_mask`.

    Raises
    ------
    ValueError
        When the points of the streamlines lie outside of the `target_mask`.

    See Also
    --------
    density_map
    """
    target_mask = np.array(target_mask, dtype=bool, copy=True)
    lin_T, offset = _mapping_to_voxel(affine)
    yield
    # End of initialization

    for sl in streamlines:
        try:
            ind = _to_voxel_coordinates(sl, lin_T, offset)
            i, j, k = ind.T
            state = target_mask[i, j, k]
        except IndexError:
            raise ValueError("streamlines points are outside of target_mask")
        if state.any() == include:
            if strict == 'strict':
                yield sl[state == include]
            elif strict == 'longstring':
                longsl=longstring(state == include,margin=2)
                yield sl[longsl]
            else:
                yield sl


def prune_streamlines(streamline, mask, cutoff=2, harshcut=None, verbose=None):
    """
    gets rid of extraneous and spurious streamlines by getting rid of voxels outside of mask
    streamline: the list of streamlines
    mask: the mask usually coressponding to a brain extracted diffusion (b0) image
    cutoff: the minimum number of voxels necessary for any one streamline to be included. Must be at least 2 for
    any evaluation of the streamlines via Life
    harshcut: if set at True, harshcut will delete any voxels of a streamline outside of a mask as well as any following
    voxel, if set to None or False, harshcut will only delete voxels outside of mask and treat subsequent voxels as normal

    returns the pruned streamline
    """
    delstream=[]
    if verbose:
        print("Starting the pruning")
    duration = time()
    if mask is not None:
        for idx,s in enumerate(streamline): #iterate through all streams
            s_vox=unique_rows(np.round(s).astype(np.intp))
            cutlist=[]
            for vox in range(np.shape(s_vox)[0]):
                if not mask[tuple(s_vox[vox])]:
                    cutlist.append(vox)         #if local mask is 0, add voxel to list of voxels to cut
            if harshcut:                        #if harshcut, destroy all voxels folloxing the out of mask voxel in streamline
                startcut = np.min(cutlist)
                cutlist = range(startcut, len(np.shape(s_vox)[0]))
                s = np.delete(s, cutlist)
            else:
                s = np.delete(s, cutlist, axis=0)   #else, cut all voxels of streamlines not in mask
            streamline[idx] = s                       #replace old streamline with new spurious streamline
            streamshape = np.shape(np.asarray(s))
            if streamshape[0] < cutoff: #minimum number of voxels required in streamline, default and minimum at 2
                delstream.append(idx)       #if number of voxels in streamline too low, cut it
                if verbose == "oververbose":
                    print("Skipped stream " + str(idx) + " out of " + str(len(streamline)) + " streamlines for being too small")
    else:
        for idx, s in enumerate(streamline):
            streamshape = np.shape(np.asarray(s))
            if streamshape[0] < cutoff:     #minimum number of voxels required in streamline, default and minimum at 2
                delstream.append(idx)
                if verbose:
                    print("Skipped stream " + str(idx) + " out of " + str(len(streamline)) + " streamlines for being too small")
    if verbose:
        print("Obtaining fiber signal process done in " + str(time() - duration) + "s")
        if len(delstream) > 0:
            print("Skipped " + str(len(delstream)) + " out of " + str(
                len(streamline)) + " due to size constraints (tiny streamlines)")
    for idx in reversed(delstream):
        streamline.pop(idx)
    return streamline

def save_trk_heavy_duty(fname, streamlines, affine, vox_size=None, shape=None, header=None):
    """ Saves tractogram files (*.trk)

    Parameters
    ----------
    fname : str
        output trk filename
    streamlines : list of 2D arrays, generator or ArraySequence
        Each 2D array represents a sequence of 3D points (points, 3).
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline points.
    vox_size : array_like (3,), optional
        The sizes of the voxels in the reference image (default: None)
    shape : array, shape (dim,), optional
        The shape of the reference image (default: None)
    header : dict, optional
        Metadata associated to the tractogram file(*.trk). (default: None)
    """
    if vox_size is not None and shape is not None:
        if not isinstance(header, dict):
            header = {}
        header[Field.VOXEL_TO_RASMM] = affine.copy()
        header[Field.VOXEL_SIZES] = vox_size
        header[Field.DIMENSIONS] = shape
        header[Field.VOXEL_ORDER] = "".join(aff2axcodes(affine))

    tractogram = nib.streamlines.LazyTractogram(streamlines)
    tractogram.affine_to_rasmm = affine
    trk_file = nib.streamlines.TrkFile(tractogram, header=header)
    nib.streamlines.save(trk_file, fname)


def unload_trk(tractogram_path, reference='same'):
    """ Similar functionality as the older version of load_trk, as it directly
    extracts the streams and header instead of returning a Tractogram object

    Parameters
    ----------
    tractogram_path: the file path of the tractogram data ( path/tract.trk )
    reference: the file used for the header information. if 'same', use the hdr from tractogram file
    """

    if reference.lower() == "same":
        print("Reference taken directly from file")
    tract_obj = load_trk(tractogram_path, reference)
    streams_control = tract_obj.streamlines
    hdr_control = tract_obj.space_attribute
    return streams_control, hdr_control, tract_obj