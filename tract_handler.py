#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Small tools for sorting through tracts, part of DTC pipeline
get_streamvals looks at tracts and catches the values in the same space on a different image (FA, MD, etc)
target catches the streamlines going through the ROI
gettrkpath catches the right file,
gettrkparams gets the parameters from a set of streamlines
get_tract_params is essentially the same but also catches the right file beforehand with gettrkpath
get_connectome_attributes is a WIP that will be used for TNPCA, obtaining more parameters based on fa, md, rd, etc
pruned_streamlines will... prune streamlines
save_roisubset saves streamlines dependent on right ROI
"""

import numpy as np
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
#from dipy.tracking.utils import unique_rows
from time import time
from functools import wraps
import pandas as pd
from BIAC_tools import isempty, send_mail
from dipy.tracking.streamline import Streamlines
from dipy.io.utils import create_tractogram_header
from tract_save import save_trk_heavy_duty
from dipy.io.streamline import load_trk
from dipy.tracking.utils import length
from dipy.viz import window, actor
import glob
import argparse
import os, shutil
from tract_save import save_trk_header
from streamline_nocheck import load_trk as load_trk_spe
import nibabel as nib
from computer_nav import glob_remote, load_trk_remote
import fnmatch


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

def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a

def catch_unique(list_in):
   # intilize an empty list
   unq_list = []
   for x in list_in:
       x = totuple(x)
       if x not in unq_list:
         unq_list.append(x)
   return(unq_list)

def array_to_tuplelist(array):
    list = []
    for x in array:
        x = totuple(x)
        list.append(x)
    return list

def get_streamvals(streamlines, affine, value_array, voxel_array = None):
    """With a set of streamlines and a value matrix,
    iterates through the set of streamlines in order to determine what is the average and maximum values found in that array
    expected to work with FA, MD, etc (must be properly registered and in same voxel space)"""

    #vals = np.ndarray(shape=(3, 0), dtype=int)
    #locations = np.ndarray(shape=(3, 0), dtype=int)
    if voxel_array is None:
        stream_voxels = []
        lin_T, offset = _mapping_to_voxel(affine)
        for sl, _ in enumerate(streamlines):
            # Convert streamline to voxel coordinates
            entire = _to_voxel_coordinates(streamlines[sl], lin_T, offset)
            stream_voxels_temp = array_to_tuplelist(entire)
            stream_voxels = stream_voxels + stream_voxels_temp

        stream_voxels_uniq = catch_unique(stream_voxels) #to ensure that voxels arent counted multiple times. Remove if voxel weight is considered a positive
        voxel_array = np.array(stream_voxels_uniq)

    i,j,k = voxel_array.T
    if np.size(value_array) == 2:
        value_array = value_array[0]
    vals = value_array[i,j,k]
    meanvals = np.mean(vals)
    maxvals = np.max(vals)

    return meanvals, maxvals, stream_voxels_uniq

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


def ratio_to_str(ratio):
    if ratio == 1:
        saved_streamlines = "_all"
    elif ratio > 1:
        saved_streamlines = "_ratio_" + str(ratio)
    else:
        raise Exception(f'invalid ratio {str(ratio)}')
    return saved_streamlines

def get_ratio(filepath):
    filename = os.path.basename(filepath)
    if '_all' in filename:
        ratio = 1
    elif '_ratio' in filename:
        ratio = filename.split('_ratio')[1].split('_')[0]
    else:
        ratio = None
    return(ratio)


def gettrkpath(trkpath, subject, str_identifier, pruned=False, verbose=False, sftp=None):

    if sftp is None:
        #Finds the right trk file based on certain established parameters (folder, subject, extra identifiers)
        if os.path.isfile(trkpath):
            if os.path.splitext(trkpath)[1] == ".trk":
                return trkpath, True
            else:
                trkpath = os.path.abspath(trkpath)
        if pruned:
            filepath = os.path.join(trkpath,subject + str_identifier + '_pruned.trk')
        else:
            filepath = os.path.join(trkpath, subject + str_identifier + '.trk')
        trkpaths = glob.glob(filepath)
        if trkpaths:
            trkfile = trkpaths[0]
            if verbose:
                print("Subject " + subject + " was found at " + trkfile)
            return trkfile, True
        else:
            print("Could not find "+filepath)
            return filepath, False
    else:
        if np.size(glob_remote(trkpath,sftp))>0:
            if os.path.splitext(trkpath)[1] == ".trk":
                return trkpath, True
            else:
                trkpath = os.path.abspath(trkpath)

        if pruned:
            filepath = os.path.join(trkpath,subject + str_identifier + '_pruned.trk')
        else:
            filepath = os.path.join(trkpath, subject + str_identifier + '.trk')

        trkpaths = glob_remote(filepath,sftp)
        if trkpaths:
            trkfile = trkpaths[0]
            if verbose:
                print("Subject " + subject + " was found at " + trkfile)
            return trkfile, True
        else:
            print("Could not find " + filepath)
            return filepath, False

def gettrkpath_testsftp(trkpath, subject, str_identifier, sftp = None, pruned=False, verbose=False):
    #Finds the right trk file based on certain established parameters (folder, subject, extra identifiers)
    if os.path.isfile(trkpath):
        if os.path.splitext(trkpath)[1] == ".trk":
            return trkpath, True
        else:
            trkpath = os.path.abspath(trkpath)
    if pruned:
        filepath = os.path.join(trkpath,subject + str_identifier + '_pruned.trk')
    else:
        filepath = os.path.join(trkpath, subject + str_identifier + '.trk')

    try:
        sftp.stat(filepath)
        if verbose:
            print("Subject " + subject + " was found at " + filepath)
        return filepath, True
    except:
        print("Could not find "+filepath)
        return filepath, False


def get_trk_params(streamlines, verbose = False):
    #Gets extra parameters from a set of streamlines and returns them (num tracts, average/min/max length, std length)
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


def viewclusters(clusters,streamlines, outpath=None, interactive=False):
    #Linked to viewing clusters. If outpath given, will save info to right location, if interactive, will show window
    colormap = actor.create_colormap(np.ravel(clusters.centroids))
    colormap_full = np.ones((len(streamlines), 3))
    for cluster, color in zip(clusters, colormap):
        colormap_full[cluster.indices] = color

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(streamlines, colormap_full))
    window.record(scene, out_path=outpath, size=(600, 600))

    # Enables/disables interactive visualization
    if interactive:
        window.show(scene)


def get_connectome_attributes(streamlines, affine, fa, md, verbose):
    #WIP, gets more parameters linked to streamlines and connectomes for TNPCA
    numtracts, minlength, maxlength, meanlength, stdlength = get_trk_params(streamlines, verbose)
    if fa is None:
        print("Fa not found")
        meanfa, maxfa = 0, 0
    else:
        meanfa, maxfa, voxel_array = get_streamvals(streamlines, affine, fa)
    if md is None:
        print("Fa not found")
        meanmd, maxmd = 0, 0
    else:
        meanmd, maxmd = get_streamvals(streamlines, affine, md, voxel_array = voxel_array)

    return numtracts, minlength, maxlength, meanlength, stdlength, meanfa, maxfa, meanmd, maxmd

def get_tract_params(mypath, subject, str_identifier='', pruned = False, verbose = False, sftp=None):

    trkpath, exists = gettrkpath(mypath, subject, str_identifier, pruned, verbose,sftp=sftp)
    if trkpath is not None:
        trkdata = load_trk_remote(trkpath, "same",sftp=sftp)
        verbose = True
        if verbose:
            print("loaded ")
        # trkdata.to_vox()
        if hasattr(trkdata, 'space_attribute'):
            header = trkdata.space_attribute
        elif hasattr(trkdata, 'space_attributes'):
            header = trkdata.space_attributes
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
            print("For subject " + subject + " the number of tracts is " + str(numtracts) + ", the minimum length is " +
                  str(minlength) + ", the maximum length is " + str(maxlength) + ", the mean length is " + str(meanlength)
                  + ", the std is " + str(stdlength))
        return subject, numtracts, minlength, maxlength, meanlength, stdlength, header, affine, trkdata
    else:
        print("Error, trkfile not found")

def streamline_checker(streamline, voxdim, verbose = False):

    from itertools import product

    mymin = np.min(streamline, axis=0)
    mymax = np.max(streamline, axis=0)
    np.asarray(list(product(*zip(mymin, mymax))))
    bbox_min = np.min(streamline, axis=0)
    bbox_max = np.max(streamline, axis=0)
    bbox_corners = np.asarray(list(product(*zip(bbox_min, bbox_max))))
    if np.any(bbox_corners[:, 0] > voxdim[0]) or \
            np.any(bbox_corners[:, 1] > voxdim[1]) or \
            np.any(bbox_corners[:, 2] > voxdim[2]):
        if verbose:
            print('deleted a streamline')
        return False
    else:
        return True

def trk_fixer(trkpath, trk_newpath, verbose=False):

    if verbose:
        t1 = time()
        txt = f'Beginning to load {trkpath}'
    trkdata = load_trk_spe(trkpath, 'same')
    if verbose:
        t2 = time()
        txt = f'It took {str(t2-t1)} seconds to load {trkpath}'
        print(txt)
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes
    remove, keep = trkdata.remove_invalid_streamlines()
    if verbose:
        t3 = time()
        duration = t3 - t2
        txt = f"it took {duration} to do the test on {trkpath}, and {np.size(remove)} streamlines were removed. saving it to {trk_newpath}"
        print(txt)
    if np.size(remove)==0:
        shutil.copy(trkpath, trk_newpath)
    else:
        new_streamlines = trkdata.streamlines
        save_trk_header(trk_newpath, new_streamlines, header, verbose=verbose)
    if verbose:
        t4 = time()
        txt = f'It took {str(t4-t3)} seconds to save {trk_newpath}'
        print(txt)

def trk_fixer_old(trkpath, trk_newpath, verbose = False):

    if verbose:
        t1 = time()
        txt = f'Beginning to load {trkpath}'
    trkdata = load_trk_spe(trkpath, 'same')
    trk_streamlines = trkdata.streamlines
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes

    voxdim = trkdata.dimensions

    if verbose:
        t2 = time()
        txt = f'It took {str(t2-t1)} seconds to load {trkpath}'
        print(txt)
    vverbose=True

    orig_size = np.shape(trk_streamlines)[0]

    trk_streamlines[:] = (x for x in trk_streamlines if streamline_checker(x, voxdim,vverbose))

    """
    for i, streamline in enumerate(trk_streamlines):
        mymin = np.min(trk_streamlines, axis=0)
        mymax = np.max(trk_streamlines, axis=0)
        np.asarray(list(product(*zip(mymin, mymax))))
        bbox_min = np.min(trk_streamlines, axis=0)
        bbox_max = np.max(trk_streamlines, axis=0)
        bbox_corners = np.asarray(list(product(*zip(bbox_min, bbox_max))))
        if np.any(bbox_corners[:, 0] > voxdim[0]) or \
        np.any(bbox_corners[:, 1] > voxdim[1]) or \
        np.any(bbox_corners[:, 2] > voxdim[2]):
            np.pop()
    """

    new_size = np.shape(trk_streamlines)[0]
    cut_streamlines = orig_size - new_size

    if verbose:
        t3 = time()
        duration = t3 - t2
        txt = f"it took {duration} to do the test on {trkpath}, and {str(cut_streamlines)} streamlines were removed. saving it to {trk_newpath}"
        print(txt)

    #np.asarray(list(product(*zip(bbox_min, bbox_max))))

    if cut_streamlines!=0:
        save_trk_header(trk_newpath, trk_streamlines, header, verbose = verbose)
    else:
        shutil.copy(trkpath, trk_newpath)

    if verbose:
        t4 = time()
        txt = f'It took {str(t4-t3)} seconds to save {trk_newpath}'
        print(txt)

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
                #if s_vox[vox][2] > 135:
                #    print("hi")
                try:
                    if not mask[tuple(s_vox[vox])]:
                        cutlist.append(vox)  # if local mask is 0, add voxel to list of voxels to cut
                        j += 1
                        num_voxel += 1
                    if np.any(vox < 0):
                        cutlist.append(vox)  # if local mask is 0, add voxel to list of voxels to cut
                        j += 1
                        num_voxel += 1
                except:
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
            j = 0
            s_vox = np.round(s).astype(np.intp)
            voxel_counter += len(s_vox)
            cutlist=[]
            num_voxel = 0
            for vox in range(np.shape(s_vox)[0]):
                if np.any(vox < 0):
                    cutlist.append(vox)  # if local mask is 0, add voxel to list of voxels to cut
                    j += 1
                    num_voxel += 1
            if num_voxel>1:
                print("cut "+num_voxel+" for being out of bounds (bounding box issue)")
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
                        save_trk_heavy_duty(trkroiminipath, streamlines=ratioed_roi_sl_gen,
                                            affine=affine, header=myheader)
                else:
                    trkdata = load_trk(trkminipath, 'same')
                    trkdata.to_vox()
                    if hasattr(trkdata, 'space_attribute'):
                        header = trkdata.space_attribute
                    elif hasattr(trkdata, 'space_attributes'):
                        header = trkdata.space_attributes
                    trkstreamlines = trkdata.streamlines


def reducetractnumber(oldtrkfile, newtrkfilepath, getdata=False, ratio=10, return_affine= False, verbose=False):

    if verbose:
        print("Beginning to read " + oldtrkfile)
    trkdata = load_trk(oldtrkfile, "same")
    if verbose:
        print("loaded the file " + oldtrkfile)
    trkdata.to_vox()
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes
    affine = trkdata._affine
    trkstreamlines = trkdata.streamlines

    ministream=[]
    for idx, stream in enumerate(trkstreamlines):
        if (idx % ratio) == 0:
            ministream.append(stream)
    del trkstreamlines
    myheader = create_tractogram_header(newtrkfilepath, *header)
    ratioed_sl = lambda: (s for s in ministream)
    save_trk_heavy_duty(newtrkfilepath, streamlines=ratioed_sl,
                                   affine=affine, header=myheader)
    if verbose:
        print("The file " + oldtrkfile + " was reduced to one "+str(ratio)+"th of its size and saved to "+newtrkfilepath)

    if getdata:
        if return_affine:
            return(ministream,affine)
        else:
            return(ministream)
    else:
        if return_affine:
            return(affine)
        else:
            return


"""
###TRK to TCK converter, courtesy of Marc Cote, https://gist.github.com/MarcCote/ea6842cc4c3950f7596fc3c8a0be0154
def build_argparser():
    DESCRIPTION = "Convert tractograms (TRK -> TCK)."
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('tractograms', metavar='bundle', nargs="+", help='list of tractograms.')
    p.add_argument('-f', '--force', action="store_true", help='overwrite existing output files.')
    return p
"""

def trktotck(trk_path, overwrite=False):

    import warnings
    try:
        tractogram = load_trk(trk_path, 'same')
    except:
        tractogram = load_trk_spe(trk_path, 'same')

    if nib.streamlines.detect_format(tractogram) is not nib.streamlines.TrkFile:
        warnings.warn("Skipping non TRK file: '{}'".format(tractogram))

    output_filename = tractogram[:-4] + '.tck'
    if os.path.isfile(output_filename) and not overwrite:
        warnings.warn("Skipping existing file: '{}'. Set overwrite to true".format(output_filename))

    trk = nib.streamlines.load(tractogram)
    nib.streamlines.save(trk.tractogram, output_filename)


def reducetractnumber_all(trkpath, str_identifier1, str_identifier2,  subject, ratio, verbose):

        trkoldpath = gettrkpath(trkpath, subject, str_identifier1 + "_pruned", verbose)
        trknewpath = gettrkpath(trkpath, subject, str_identifier2 + "_pruned", verbose)
        if trknewpath is None:
            trknewpath = (trkpath + '/' + subject + "" + str_identifier2 + '.trk')
            reducetractnumber(trkoldpath,trknewpath, getdata=False, ratio=ratio, return_affine=False, verbose=True)
        else:
            if verbose:
                txt = ("Subject "+ subject +" is already done")
                send_mail(txt, subject="reduce code in process")
                print(txt)
