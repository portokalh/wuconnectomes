#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Small tools for sorting through tracts
"""

import numpy as np
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
#from dipy.tracking.utils import unique_rows
from time import time
from functools import wraps
import pandas as pd
from BIAC_tools import isempty
from dipy.tracking.streamline import Streamlines
from dipy.io.utils import create_tractogram_header
from tract_save import save_trk_heavy_duty
from dipy.io.streamline import load_trk
from dipy.tracking.utils import length
import os

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


def get_trk_params(streamlines, verbose = False):
    #trkdata = load_trk(trkpath, "same")
    if verbose:
        print("loaded ")
    # trkdata.to_vox()
    #header = trkdata.space_attribute
    #affine = trkdata._affine
    lengths = length(streamlines)
    del streamlines
    # lengths = list(length(trkstreamlines))
    lengths = list(lengths)
    numtracts = np.size(lengths)
    minlength = np.min(lengths)
    maxlength = np.max(lengths)
    meanlength = np.mean(lengths)
    stdlength = np.std(lengths)
    return numtracts, minlength, maxlength, meanlength, stdlength


def get_tract_params(mypath, subject, str_identifier, verbose = False):

    trkpath = gettrkpath(mypath, subject, str_identifier, verbose)
    trkdata = load_trk(trkpath, "same")
    verbose = True
    if verbose:
        print("loaded ")
    # trkdata.to_vox()
    header = trkdata.space_attribute
    affine = trkdata._affine
    lengths = length(trkdata.streamlines)
    #del trkdata
    # lengths = list(length(trkstreamlines))
    lengths = list(lengths)
    numtracts = np.size(lengths)
    minlength = np.min(lengths)
    maxlength = np.max(lengths)
    meanlength = np.mean(lengths)
    stdlength = np.std(lengths)
    if verbose:
        print("For subject " + subject + " the number of tracts is " + numtracts + ", the minimum length is " + minlength + ", the maximum length is " + maxlength + ", the mean length is " + meanlength + ", the std is " + stdlength)
    return subject, numtracts, minlength, maxlength, meanlength, stdlength, header, affine, trkdata


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
    voxel_counter = 0
    outmask_counter = 0
    if mask is not None:
        num_voxel=0
        for idx,s in enumerate(streamline): #iterate through all streams
            j = 0
            s_vox = np.round(s).astype(np.intp)
            voxel_counter += len(s_vox)
            cutlist=[]
            for vox in range(np.shape(s_vox)[0]):
                #if np.array_equal(s_vox[vox], [38, 149, 167]):
                #    print(mask[tuple(s_vox[vox])])
                try:
                    if not mask[tuple(s_vox[vox])]:
                        cutlist.append(vox)  # if local mask is 0, add voxel to list of voxels to cut
                        j += 1
                        num_voxel += 1
                except IndexError:
                    cutlist.append(vox)         #if local mask is 0, add voxel to list of voxels to cut
                    j += 1
                    num_voxel += 1
                    outmask_counter += 1
                    print("Out of bounds streamline")
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
                    print("Skipped stream " + str(idx) + " out of " + str(len(streamline)) + " streamlines for being too small after cutting off "+ str(j) +" voxels")

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
        if num_voxel > 0:
            print("Dropped " + str(num_voxel) + " voxels for being out of mask")
    if outmask_counter > 0:
        print("Found " + str(outmask_counter) + " voxels out of mask")
    for idx in reversed(delstream):
        streamline.pop(idx)
    return streamline


def save_roisubset(streamlines, roislist, roisexcel, labelmask, stringstep, ratios, trkpath, subject, affine, header):
    
    #atlas_legends = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"
    
    df = pd.read_excel(roisexcel, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()    
    
    for rois in roislist:
        if len(rois)==1:
            roiname = "_" + rois[0] + "_"
        elif len(rois)>1:
            roiname="_"
            for roi in rois:
                roiname = roiname + roi[0:4]
            roiname = roiname + "_"    
            
        labelslist=[]#fimbria

        for roi in rois:
            rslt_df = df.loc[df['Structure'] == roi.lower()]
            if roi.lower() == "wholebrain" or roi.lower() == "brain":
                labelslist=None
            else:
                labelslist=np.concatenate((labelslist,np.array(rslt_df.index2)))

        if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
            txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
            print(txt)

#bvec_orient=[-2,1,3]    
    
        if isempty(labelslist):
            roimask = np.where(labelmask == 0, False, True)
        else:
            if labelmask is None:
                raise ("Bad label data, could not define ROI for streams")
            roimask = np.zeros(np.shape(labelmask),dtype=int)
            for label in labelslist:
                roimask = roimask + (labelmask == label)
        
        if not isempty(labelslist):
            trkroipath = trkpath + '/' + subject + roiname + "_stepsize_" + stringstep + '.trk'
            if not os.path.exists(trkroipath):
                affinetemp=np.eye(4)
                trkstreamlines = target(streamlines, affinetemp, roimask, include=True, strict="longstring")
                trkstreamlines = Streamlines(trkstreamlines)
                myheader = create_tractogram_header(trkroipath, *header)
                roi_sl = lambda: (s for s in trkstreamlines)
                save_trk_heavy_duty(trkroipath, streamlines=roi_sl,
                            affine=affine, header=myheader)
            else:
                trkdata = load_trk(trkroipath, 'same')
                trkdata.to_vox()
                if hasattr(trkdata, 'space_attribute'):
                    header = trkdata.space_attribute
                elif hasattr(trkdata, 'space_attributes'):
                    header = trkdata.space_attributes
                trkstreamlines = trkdata.streamlines
                
        for ratio in ratios:
            if ratio != 1:
                trkroiminipath = trkpath + '/' + subject + '_ratio_' + ratios + roiname + "_stepsize_" + stringstep + '.trk'
                if not os.path.exists(trkroiminipath):
                    ministream = []
                    for idx, stream in enumerate(trkstreamlines):
                        if (idx % ratio) == 0:
                            ministream.append(stream)
                    trkstreamlines = ministream
                    myheader = create_tractogram_header(trkminipath, *header)
                    ratioed_roi_sl_gen = lambda: (s for s in trkstreamlines)
                    if allsave:
                        tract_save.save_trk_heavy_duty(trkroiminipath, streamlines=ratioed_roi_sl_gen,
                                            affine=affine, header=myheader)
                else:
                    trkdata = load_trk(trkminipath, 'same')
                    trkdata.to_vox()
                    if hasattr(trkdata, 'space_attribute'):
                        header = trkdata.space_attribute
                    elif hasattr(trkdata, 'space_attributes'):
                        header = trkdata.space_attributes
                    trkstreamlines = trkdata.streamlines