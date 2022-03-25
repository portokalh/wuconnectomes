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
from dipy.io.utils import create_tractogram_header
import time, os
import numpy as np

def make_tractogram_object(fname, streamlines, affine, vox_size=None, shape=None, header=None):
    """ Saves tractogram object for future use

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
    return tractogram, trk_file

def save_trk_header(filepath, streamlines, header, affine=np.eye(4), fix_streamlines = False, verbose=False):

    myheader = create_tractogram_header(filepath, *header)
    trk_sl = lambda: (s for s in streamlines)
    if verbose:
        print(f'Saving streamlines to {filepath}')
        time1 = time.perf_counter()
    save_trk_heavy_duty(filepath, streamlines=trk_sl,
                        affine=affine, header=myheader, fix_streamlines = fix_streamlines, return_tractogram=False)
    if verbose:
        time2 = time.perf_counter()
        print(f'Saved in {time2 - time1:0.4f} seconds')

def save_trk_heavy_duty(fname, streamlines, affine, vox_size=None, shape=None, header=None, fix_streamlines = False,
                        return_tractogram=False,sftp=None):
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
    if fix_streamlines:
        tractogram.remove_invalid_streamlines()
    trk_file = nib.streamlines.TrkFile(tractogram, header=header)
    if sftp is not None:
        nib.streamlines.save(trk_file, fname)
    else:
        temp_path = f'{os.path.join(os.path.expanduser("~"), os.path.basename(fname))}'
        nib.streamlines.save(trk_file, temp_path)
        try:
            sftp.put(temp_path, fname)
            os.remove(temp_path)
        except Exception as e:
            os.remove(temp_path)
            raise Exception(e)
    if return_tractogram:
        return tractogram, trk_file

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
    try:
        hdr_control = tract_obj.space_attribute
    except:
        hdr_control = tract_obj.space_attributes
    return streams_control, hdr_control, tract_obj