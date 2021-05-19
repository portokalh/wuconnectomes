#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:34:51 2020

@author: alex
Now part of DTC pipeline. Used to save and unload heavy trk files
"""


import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from dipy.io.streamline import load_trk

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