#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:48:38 2020

@author: Jacques Stout

Part of the DTC pipeline
This is the big one. Almost all the DTC launchers call to it, and it essentially contains most if not all of the
main 'steps', such as 'get diffusion preprocessing, launch the tractography, start the connectomes' without having those
functions themselves. Essentially a 'transition file' between the launcher and the processing functions.
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import nibabel as nib
import numpy as np
#from dipy.io.streamline import load_trk
from streamline_nocheck import load_trk
from os.path import splitext
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
import pickle
import pandas as pd
import smtplib
import pathlib
from dipy.denoise.localpca import mppca
from dif_to_trk import QCSA_tractmake
from file_tools import mkcdir
import os, re, sys, io, struct, socket, datetime
from email.mime.text import MIMEText
import glob
from bvec_handler import extractbvals
from dipy.tracking.utils import unique_rows
from numpy import ravel_multi_index

from collections import defaultdict, OrderedDict

from time import time
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from tract_save import save_trk_heavy_duty
from dipy.tracking.local_tracking import LocalTracking
from dipy.direction import peaks
from nibabel.streamlines import detect_format
from dipy.io.utils import (create_tractogram_header)
from dipy.viz import window, actor
from dipy.tracking.utils import connectivity_matrix
from dipy.segment.mask import segment_from_cfa
from dipy.segment.mask import bounding_box

import multiprocessing
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
import multiprocessing.pool


from dipy.reconst import shm
from scipy.ndimage.morphology import binary_dilation
from dipy.tracking import utils
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion, ThresholdStoppingCriterion

from dipy.tracking.streamline import Streamlines
import matplotlib.pyplot as plt

#from dipy.denoise.localpca import mppca


from random import randint

from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib
import matplotlib.pyplot as plt
import xlsxwriter

from figures_handler import shore_scalarmaps
from tract_eval import bundle_coherence, LiFEvaluation
from BIAC_tools import send_mail, isempty
from tract_handler import target, prune_streamlines, get_trk_params, get_tract_params, gettrkpath, reducetractnumber, reducetractnumber_all
from nifti_handler import getfa, getdiffdata_all, getdiffdata, getdiffpath, getgtab, getlabelmask, move_bvals, getmask, getb0s, getrefdata
from connectome_handler import connectivity_matrix_custom
from diff_preprocessing import dwi_to_mask, denoise_pick, make_tensorfit
import tract_save
import xlrd
import warnings
import shutil
#import dipy.reconst.dki as dki
import dipy.reconst.msdki as msdki

from multiprocessing import Pool
from convert_atlas_mask import convert_labelmask, atlas_converter
#from connectivity_own import connectivity_matrix_special
from excel_management import connectomes_to_excel, grouping_to_excel
from computer_nav import load_trk_remote, checkfile_exists_remote, glob_remote, load_nifti_remote, pickledump_remote, \
    remove_remote

def strfile(string):
    # Converts strings into more usable 'file strings (mostly takes out . and turns it into _
    if string == 'any':
        return ''
        # if the string is any, that means that it is unspecified in filepath, therefore null
    else:
        try:
            string = str(string)
            string = string.replace(".", "_")
            return(string)
        except AttributeError:
            raise AttributeError("strfile error: not a usable number or string ")

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_str_identifier(stepsize, ratio, trkroi):

    if ratio == 1:
        saved_streamlines = "_all"
    else:
        saved_streamlines = "_ratio_" + str(ratio)

    if len(trkroi) == 1:
        roistring = "_" + trkroi[0]  # + "_"
    elif len(trkroi) > 1:
        roistring = "_"
        for roi in trkroi:
            roistring = roistring + roi[0:4]
        roistring = roistring  # + "_"
    str_identifier = '_stepsize_' + str(stepsize).replace(".", "_") + saved_streamlines + roistring
    return str_identifier

def almicedf_fix(df, verbose=None):
    # masterFile='/Users/alex/AlexBadea_MyPapers/FIGURES/mwm/mwm_master_organized.csv'

    df = df.replace({'runno': {"N54716/N54915": "N54915", "N54714/N54916": "N54916",
                               "N54715/N54917": "N54917", "N54718//N54721": "N54718",
                               "N54760/N54754/N54755": "N54760", "N54757/N54759": "N54759",
                               "N54805/N54800": "N54805", "N54891/N54900LRspecific": "N54891"}})
    df = df.replace({'Day': {"ProbeTrial1": "ProbeTrial", "ProbeTrial2": "ProbeTrial"}})
    df = df.dropna(subset=['runno'])

    alldays = df.Day.unique()
    df = df.dropna()

    if verbose:
        for day in alldays:
            df_day = df.loc[df['Day'] == day]
            sf_day = df_day.groupby('runno')['Acq Number'].nunique()
            print('For the Day: ' + day + ' those animals had less than the standard number of tests:')
            print(sf_day.where(sf_day < np.max(sf_day)).dropna())

    return (df)
    # na_columns=df.columns[df.isna().any()]
    # df_withna=df[df.isnull().any(axis=1)][:].head()
    # df = mice_day2.groupby('runno')['Acq Number'].nunique()

def testsnr():

    corpus_mask = np.where(labels == 121, 1, 0) + np.where(labels == 1121, 1, 0)
    # meed to change threshold, possibly ROI, that better corresponds to mice (pick area with high FA)
    # probably should investigate area with
    threshold = (0.6, 1, 0, 0.1, 0, 0.1)
    mask_cc_part, cfa = segment_from_cfa(tensor_fit,corpus_mask,threshold,return_cfa = True)
    cfa_img = nib.Nifti1Image((cfa * 255).astype(np.uint8), affine)
    mask_cc_part_img = nib.Nifti1Image(corpus_mask.astype(np.uint8), affine)
    nib.save(mask_cc_part_img, '/Users/alex/jacques/mask_CC_part.nii.gz')

    region = 30
    fig = plt.figure('Corpus callosum segmentation')

    plt.subplot(1, 3, 1)
    plt.title("Corpus callosum (CC)")
    plt.axis('off')
    red = cfa[..., 0]
    plt.imshow(np.rot90(corpus_mask[region, ...]))

    plt.subplot(1, 3, 2)
    plt.title("Corpus callosum (CC)")
    plt.axis('off')
    red = cfa[..., 0]
    plt.imshow(np.rot90(red[region, ...]))

    plt.subplot(1, 3, 3)
    plt.title("CC mask used for SNR computation")
    plt.axis('off')
    plt.imshow(np.rot90(mask_cc_part[region, ...]))
    fig.savefig("CC_segmentation.png", bbox_inches='tight')

    mean_signal = np.mean(data[mask_cc_part], axis=0)
    from scipy.ndimage.morphology import binary_dilation
    mask_noise = binary_dilation(mask, iterations=10)
    mask_noise[..., :mask_noise.shape[-1] // 2] = 1
    mask_noise = ~mask_noise
    mask_noise_img = nib.Nifti1Image(mask_noise.astype(np.uint8), affine)
    nib.save(mask_noise_img, 'mask_noise.nii.gz')

    noise_std = np.std(data[mask_noise, :])
    print('Noise standard deviation sigma= ', noise_std)

    # Exclude null bvecs from the search
    idx = np.sum(gtab.bvecs, axis=-1) == 0
    gtab.bvecs[idx] = np.inf
    axis_X = np.argmin(np.sum((gtab.bvecs - np.array([1, 0, 0])) ** 2, axis=-1))
    axis_Y = np.argmin(np.sum((gtab.bvecs - np.array([0, 1, 0])) ** 2, axis=-1))
    axis_Z = np.argmin(np.sum((gtab.bvecs - np.array([0, 0, 1])) ** 2, axis=-1))

    for direction in [0, axis_X, axis_Y, axis_Z]:
        SNR = mean_signal[direction] / noise_std
        if direction == 0:
            print("SNR for the b=0 image is :", SNR)
        else:
            print("SNR for direction", direction, " ",
                  gtab.bvecs[direction], "is :", SNR)


def analysis_diffusion_figures(dwipath,outpath,subject, bvec_orient=[1,2,3], verbose=None):

    if verbose:
        print('Running the ' + subject + ' file')

    diff_data, affine, gtab, labelmask, vox_size, fdwipath, _, ref_info = getdiffdata_all(dwipath, subject, bvec_orient)
    outpath_fig = outpath + subject + "SHORE_maps.png"
    shore_scalarmaps(diff_data, gtab, outpath_fig, verbose = verbose)

def dwiconnectome_analysis(dwipath,outpath,subject, whitematter_labels, targetrois, labelslist, bvec_orient=[1,2,3], verbose=None):

    diff_data, affine, gtab, labelmask, vox_size, fdwipath, header, ref_info = getdiffdata_all(dwipath, subject, bvec_orient, verbose)
    #from dipy.data import read_stanford_labels, read_stanford_t1
    #hardi_img, gtab_hardi, labels_hardi = read_stanford_labels()
    #data = hardi_img.get_data()
    #labels = labels_img.get_data()

    #t1 = read_stanford_t1()
    #t1_data = t1.get_data()

    import numpy as np
    whitemask = np.zeros(np.shape(labelmask),dtype=int)
    for label in whitematter_labels:
        whitemask = whitemask + (labelmask == label)

    roimask = np.zeros(np.shape(labelmask),dtype=int)
    for label in labelslist:
        roimask = roimask + (labelmask == label)

    white_matter = binary_dilation(whitemask)
    csamodel = shm.CsaOdfModel(gtab, 6)
    csapeaks = peaks.peaks_from_model(model=csamodel,
                                      data=diff_data,
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

    cc_slice = roimask == 1
    #cc_slice = roilabels
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

    #vol_actor = actor.slicer(t1_data)
    vol_actor = actor.slicer(diff_data[:,:,:,0])


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
    targetrois = targetrois[0]
    window.record(r, n_frames=1, out_path=outpath + subject + targetrois + '_axial.png',
                  size=(800, 800))
    if interactive:
        window.show(r)
    r.set_camera(position=[-1, 0, 0], focal_point=[0, 0, 0], view_up=[0, 0, 1])
    window.record(r, n_frames=1, out_path=outpath + subject + targetrois + '_sagittal.png',
                  size=(800, 800))
    if interactive:
        window.show(r)

    atlas_legends = "/Users/alex/jacques/atlases/CHASSSYMM3AtlasLegends.xlsx"

    M, grouping = utils.connectivity_matrix(cc_streamlines, affine, labelmask,
                                            return_mapping=True,
                                            mapping_as_streamlines=True)
    M2, grouping2 = utils.connectivity_matrix(cc_streamlines, affine, labelmask,
                                            return_mapping=True,
                                            mapping_as_streamlines=True)

    M[:3, :] = 0
    M[:, :3] = 0

    import numpy as np
    import matplotlib.pyplot as plt
    plt.imshow(np.log1p(M), interpolation='nearest')
    plt.savefig(outpath + subject + "connectivity.png")

    lr_superiorfrontal_track = grouping[11, 54]
    shape = labelmask.shape
    dm = utils.density_map(lr_superiorfrontal_track, affine, shape)

    import nibabel as nib
    from dipy.io.stateful_tractogram import Space, StatefulTractogram
    from dipy.io.streamline import save_trk

    # Save density map
    dm_img = nib.Nifti1Image(dm.astype("int16"), affine)
    dm_img.to_filename(outpath + subject + "lr-superiorfrontal-dm.nii.gz")

    lr_sf_trk = Streamlines(lr_superiorfrontal_track)

    # Save streamlines
    sft = StatefulTractogram(lr_sf_trk, dm_img, Space.VOX)
    save_trk(sft, outpath + subject + "lr-superiorfrontal.trk")


def deprecation(message):
    warnings.warn(message, DeprecationWarning, stacklevel=2)

def makemask_fromdiff(dwipath, subject, bvec_orient):

    data, affine, _, mask, vox_size, fdwipath, header, ref_info = getdiffdata_all(dwipath, subject, bvec_orient)
    outpath = os.path.join(dwipath, subject)
    mask = dwi_to_mask(data, subject, affine, outpath, makefig = False)
    return(mask)


"""
def prunestreamline(trkorigpath, trkprunepath, cutoff = 4, forcestart = False):

    #trkprunepath = trkpath + '/' + subject + str_identifier + '_pruned.trk'
    #trkpaths = glob.glob(trkpath + '/' + subject + '_' + tractsize + strproperty + 'stepsize_' + str(stepsize) + '.trk')

    if not os.path.isfile(trkprunepath) or forcestart:

        trkdata = load_trk(trkorigpath, "same")
        affine = trkdata._affine
        trkdata.to_vox()
        trkstreamlines = trkdata.streamlines
        if hasattr(trkdata,'space_attribute'):
            header = trkdata.space_attribute
        elif hasattr(trkdata,'space_attributes'):
            header = trkdata.space_attributes

        cutoff=4
        pruned_streamlines = prune_streamlines(list(trkstreamlines), labelmask, cutoff=cutoff, verbose=verbose)
        del(trkdata)
        pruned_streamlines_SL = Streamlines(pruned_streamlines)
        myheader = create_tractogram_header(trkprunepath, *header)
        prune_sl = lambda: (s for s in pruned_streamlines)
        if prunesave:
            tract_save.save_trk_heavy_duty(trkprunepath, streamlines=prune_sl, affine=affine, header=myheader)
        return(pruned_streamlines)
"""

def excel_extract(roi_path):

    #from __future__ import print_function
    from os.path import join, dirname, abspath
    xl_workbook = xlrd.open_workbook(roi_path)

    xl_sheet = xl_workbook.sheet_by_index(0)
    print('Sheet name: %s' % xl_sheet.name)
    data = np.zeros(xl_sheet.nrows-1,xl_sheet.nrows-1)
    num_cols = xl_sheet.ncols  # Number of columns
    for row_idx in range(1, xl_sheet.nrows):  # Iterate through rows
        for col_idx in range(1, num_cols):  # Iterate through columns
            cell_obj = xl_sheet.cell(row_idx, col_idx)  # Get cell object by row, col
            data[row_idx-1,col_idx-1]=cell_obj
    return data

def copylabels(folder1, folder2, subject, verbose = False):
    _, _, labelspath = getlabelmask(folder1, subject, verbose)
    filename = os.path.basename(labelspath)
    newlabelspath = os.path.join(folder2, filename)
    if os.path.isfile(newlabelspath):
        return
    else:
        shutil.copy2(labelspath, newlabelspath)

def makedir(dir):
    if os.path.isdir(dir):
        return
    else:
        os.mkdir(dir)

def tract_connectome_analysis(diffpath, trkpath, str_identifier, outpath, subject, ROI_excel, bvec_orient, masktype = "T1",
                              inclusive = False, function_processes = 1, overwrite = False, picklesave = True, labeltype='orig',
                              symmetric = True, reference_weighting_type = None, volume_weighting=False, verbose = None, sftp = None):

    picklepath_connect = os.path.join(outpath, subject + str_identifier + '_connectomes.p')
    connectome_xlsxpath = os.path.join(outpath, subject + str_identifier + "_connectomes.xlsx")
    picklepath_grouping = os.path.join(outpath, subject + str_identifier + '_grouping.p')
    grouping_xlsxpath = os.path.join(outpath, subject + str_identifier + "_grouping.xlsx")
    picklepath_index = os.path.join(outpath, subject + str_identifier + '_index_struct.p')

    if reference_weighting_type is not None:
        reference_weight, ref_affine, _, ref_fpath, _, _ = getrefdata(diffpath, subject, reference_weighting_type, verbose, sftp=sftp)
        if reference_weight is None:
            txt = f'Subject is {subject} did not have reference {reference_weighting_type} at {diffpath}'
            warnings.warn(txt)
            return None
    else:
        reference_weight = None

    if volume_weighting:
        picklepath_connect_vol = os.path.join(outpath, subject + str_identifier + '_volweighted' + '_connectomes.p')
        connectome_xlsxpath_vol = os.path.join(outpath, subject + str_identifier + '_volweighted' + "_connectomes.xlsx")
    if reference_weighting_type is not None:
        picklepath_connect_ref = os.path.join(outpath, subject + str_identifier + '_' + reference_weighting_type + '_connectomes.p')
        connectome_xlsxpath_ref = os.path.join(outpath, subject + str_identifier + '_' + reference_weighting_type + "_connectomes.xlsx")
    if volume_weighting and reference_weighting_type is not None:
        picklepath_connect_volref = os.path.join(outpath, subject + str_identifier + '_volweighted' + '_' + reference_weighting_type + '_connectomes.p')
        connectome_xlsxpath_volref = os.path.join(outpath, subject + str_identifier + '_volweighted' + '_' + reference_weighting_type + "_connectomes.xlsx")

    mkcdir(outpath, sftp=sftp)
    print(f'Volume weighting is {volume_weighting}, reference weighting is {reference_weighting_type}\ninclusive is {inclusive}, symmetric is {symmetric}, saving connectomes at {outpath}')

    if sftp is None and os.path.exists(picklepath_connect) and os.path.exists(connectome_xlsxpath) and os.path.exists(grouping_xlsxpath) and not overwrite:
        print(f"The writing of pickle and excel of {str(subject)} at {connectome_xlsxpath} is already done")
        return
    if sftp is not None and checkfile_exists_remote(picklepath_connect,sftp) and checkfile_exists_remote(connectome_xlsxpath,sftp) and checkfile_exists_remote(grouping_xlsxpath,sftp) and not overwrite:
        print(f"The writing of pickle and excel of {str(subject)} at {connectome_xlsxpath} is already done")
        return

    mask, affinemask = getmask(diffpath,subject,masktype,verbose,sftp=sftp)
    trkfilepath, trkexists = gettrkpath(trkpath, subject, str_identifier, pruned = False, verbose = verbose, sftp = sftp)
    trkprunepath, trkpruneexists = gettrkpath(trkpath, subject, str_identifier, pruned = True, verbose = verbose, sftp = sftp)
    labelmask, labelaffine, labelpath = getlabelmask(diffpath, subject, verbose, sftp=sftp)
    if mask is None:         # Build Brain Mask
        if masktype == 'dwi':
            if verbose:
                print("Beginning to read the dwifile of subject " + subject + " at " + diffpath)
            diff_data, diffaffine, _, diff_fpath, _, _ = getdiffdata(diffpath, subject, verbose,sftp=sftp)
            if verbose:
                print("loaded the file " + diff_fpath)
            outpathmask = str(pathlib.Path(diff_fpath).parent.absolute())
            if verbose:
                print("Creating mask for subject " + subject + " at " + outpathmask)
            vol_idx = getb0s(diffpath, subject,sftp=sftp)
            mask, _ = dwi_to_mask(diff_data, subject, diffaffine, outpathmask, makefig=False, vol_idx=vol_idx, median_radius=5,
                                  numpass=6, dilate=2)
    mypath = diffpath

    import numpy as np
    prunesave = True
    if np.size(np.shape(labelmask)) == 1:
        labelmask = labelmask[0]
    if np.size(np.shape(labelmask)) == 4:
        labelmask = labelmask[:, :, :, 0]
    print("Mask shape is " + str(np.shape(labelmask)))
    cutoff = 2

    converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(ROI_excel)


    print('what is the labeltype: ')
    print(labeltype)
    if labeltype == 'combined':
        labeloutpath = labelpath.replace('.nii.gz','_comb.nii.gz')
        if not checkfile_exists_remote(labeloutpath, sftp):
            labelmask = convert_labelmask(labelmask, converter_comb, atlas_outpath=labeloutpath,
                                          affine_labels=labelaffine, sftp=sftp)
        else:
            labelmask, labelaffine = load_nifti(labeloutpath)
        index_to_struct = index_to_struct_comb
    elif labeltype == 'lrordered':
        labeloutpath = labelpath.replace('.nii.gz','_lr_ordered.nii.gz')
        if not checkfile_exists_remote(labeloutpath,sftp):
            labelmask = convert_labelmask(labelmask, converter_lr, atlas_outpath=labeloutpath,
                                      affine_labels=labelaffine, sftp=sftp)
        else:
            labelmask, labelaffine, _, _, _ = load_nifti_remote(labeloutpath,sftp)
        index_to_struct = index_to_struct_lr
    else:
        raise TypeError("Cannot recognize label type (this error raise is a THEORETICALLY a temp patch")

    #pickle.dump(index_to_struct, open(picklepath_index, "wb"))

    if (trkfilepath is not None and trkpruneexists is False and prunesave):

        if 'dwi_data' not in locals():
            if verbose:
                print("Beginning to read the dwifile of subject " + subject + " at "+diffpath)
            dwi_data, _, _, diffpath, _, _ = getdiffdata(diffpath, subject, verbose,sftp=sftp)
            if verbose:
                print("loaded the file " + diffpath)

        if verbose:
            print("Beginning to read " + trkfilepath)
        if sftp is None:
            trkdata = load_trk(trkfilepath, "same")
        else:
            trkdata = load_trk_remote(trkfilepath, "same",sftp=sftp)
        if verbose:
            print("loaded the file " + trkfilepath)
        print(trkfilepath)
        affine = trkdata._affine
        trkdata.to_vox()
        trkstreamlines = trkdata.streamlines
        trkprunepath = os.path.dirname(trkfilepath) + '/' + subject + str_identifier + '_pruned.trk'

        if np.size(np.shape(dwi_data)) == 1:
            dwi_data = dwi_data[0]
        if np.size(np.shape(dwi_data)) == 4:
            dwi_data = dwi_data[:, :, :, 0]
        print("Mask shape is " + str(np.shape(dwi_data)))

        pruned_streamlines = prune_streamlines(list(trkstreamlines), mask, cutoff=cutoff, verbose=verbose)
        pruned_streamlines_SL = Streamlines(pruned_streamlines)
        if hasattr(trkdata,'space_attribute'):
            header = trkdata.space_attribute
        elif hasattr(trkdata,'space_attributes'):
            header = trkdata.space_attributes
        myheader = create_tractogram_header(trkprunepath, *header)
        prune_sl = lambda: (s for s in pruned_streamlines)
        if prunesave:
            tract_save.save_trk_heavy_duty(trkprunepath, streamlines=prune_sl, affine=affine, header=myheader,sftp=sftp)
        del(prune_sl,pruned_streamlines,trkdata)
    elif trkpruneexists:
        if verbose:
            print("Beginning to read " + trkprunepath)
        #trkprunepath = '/Volumes/Data/Badea/Lab/human/AD_Decode/Analysis/TRK_MPCA_100/S02715_ratio100_badrun.trk'
        #trkprunedata = load_trk(trkprunepath, "same")
        if sftp is None:
            trkprunedata = load_trk(trkprunepath, "same")
        else:
            trkprunedata = load_trk_remote(trkprunepath, "same",sftp=sftp)
        if verbose:
            print("loaded the file " + trkprunepath)
        affine = trkprunedata._affine
        trkprunedata.to_vox()
        pruned_streamlines_SL = trkprunedata.streamlines

        streamlines_test = list(pruned_streamlines_SL)
        from figures_handler import viewstreamlines_anat
        #viewstreamlines_anat(streamlines_test, mask, affinemask, ratio=100, threshold=10., verbose=False)
        """
        endpoints = [sl[0::len(sl) - 1] for sl in streamlines_test]
        lin_T, offset = _mapping_to_voxel(affine)
        endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)
        i, j, k = endpoints.T
        try:
            mask[i, j, k]
        except:
            pruned_streamlines = prune_streamlines(list(pruned_streamlines_SL), mask, cutoff=cutoff,
                                                   verbose=False)
            pruned_streamlines_SL = Streamlines(pruned_streamlines)
        """
        del(trkprunedata)

    cutoff = 4
    #pruned_streamlines = prune_streamlines(list(pruned_streamlines_SL), labelmask, cutoff=cutoff, verbose=verbose)
    #pruned_streamlines_SL = Streamlines(pruned_streamlines)

    affine_streams = np.eye(4)
    t = time()

    if function_processes > 1:
        listcut = []

        n = function_processes
        size_SL = np.size(pruned_streamlines_SL)
        listcut.append(0)
        for i in np.arange(n - 1):
            listcut.append(np.int(((i + 1) * size_SL) / n))
            print(size_SL, i + 1, n)
        listcut.append(size_SL)
        print(listcut)
        pruned_cut = []
        for i in np.arange(n):
            pruned_cut.append(pruned_streamlines_SL[listcut[i]:listcut[i+1]])
        pool = Pool()
        return_mapping = True
        mapping_as_streamlines = False

        if verbose:
            print("The streamline is split into "+str(function_processes)+" of size "+str(np.int(size_SL / n)))

        print(f'This run has inclusive: {inclusive} and symmetric: {symmetric}')
        connectomic_results = pool.starmap_async(connectivity_matrix_custom, [(Streamlines(pruned_streamlines_SL[listcut[i]:listcut[i+1]]), affine_streams, labelmask,
                                                                              inclusive, symmetric, return_mapping,
                                                                              mapping_as_streamlines,reference_weight,volume_weighting) for i in np.arange(n)]).get()

        matrix =  np.zeros(np.shape(connectomic_results[0][0]))
        matrix_vol = np.zeros(np.shape(connectomic_results[0][0]))
        matrix_refweighted = np.zeros(np.shape(connectomic_results[0][0]))
        matrix_vol_refweighted = np.zeros(np.shape(connectomic_results[0][0]))
        grouping = {}
        i=0
        for connectome_results in connectomic_results:
            matrix += connectome_results[0]
            if volume_weighting:
                matrix_vol += connectome_results[1]
            if reference_weight is not None:
                matrix_refweighted += connectome_results[2]
            if volume_weighting and reference_weight is not None:
                matrix_vol_refweighted += connectome_results[3]

            for key, val in connectome_results[4].items():
                if key in grouping:
                    grouping[key].extend([j+listcut[i] for j in val])
                else:
                    grouping[key] = val
            i = i + 1
    else:
        matrix, matrix_vol, matrix_refweighted, matrix_vol_refweighted, grouping = connectivity_matrix_custom(pruned_streamlines_SL, affine_streams, labelmask, inclusive=inclusive, symmetric=symmetric,
                                            return_mapping=True,
                                            mapping_as_streamlines=False, reference_weighting = reference_weight, volume_weighting=volume_weighting)

        n = 1

    if verbose:
        print("Time taken for the accelerated calculation with " + str(n) + " processes " + str(- t + time()))

    if picklesave:
        pickledump_remote(matrix, picklepath_connect, sftp=sftp)
        pickledump_remote(grouping, picklepath_grouping, sftp=sftp)
        if volume_weighting:
            pickledump_remote(matrix_vol, picklepath_connect_vol, sftp=sftp)
        if reference_weighting_type is not None:
            pickledump_remote(matrix_refweighted, picklepath_connect_ref, sftp=sftp)
        if volume_weighting and reference_weighting_type is not None:
            pickledump_remote(matrix_vol_refweighted, picklepath_connect_volref, sftp=sftp)

        if verbose:
            txt = ("The connectomes were saved at "+picklepath_connect)
            send_mail(txt, subject="Pickle save")
            print(txt)

    matrix_sl = np.empty(np.shape(matrix), dtype=object)
    for i in np.arange(np.shape(matrix_sl)[0]):
        for j in np.arange(np.shape(matrix_sl)[1]):
            matrix_sl[i, j] = []
    for key in grouping.keys():
        matrix_sl[key] = grouping[key]
        matrix_sl[tuple(np.flip(key))] = grouping[key]

    matrix = np.delete(matrix, 0, 0)
    matrix = np.delete(matrix, 0, 1)
    matrix_vol = np.delete(matrix_vol, 0, 0)
    matrix_vol = np.delete(matrix_vol, 0, 1)
    matrix_refweighted = np.delete(matrix_refweighted, 0, 0)
    matrix_refweighted = np.delete(matrix_refweighted, 0, 1)
    matrix_vol_refweighted = np.delete(matrix_vol_refweighted, 0, 0)
    matrix_vol_refweighted = np.delete(matrix_vol_refweighted, 0, 1)
    matrix_sl = np.delete(matrix_sl, 0, 0)
    matrix_sl = np.delete(matrix_sl, 0, 1)

    connectomes_to_excel(matrix, index_to_struct, connectome_xlsxpath, overwrite=overwrite, verbose=verbose, sftp=sftp)
    if volume_weighting:
        connectomes_to_excel(matrix_vol, index_to_struct, connectome_xlsxpath_vol, overwrite=overwrite, verbose=verbose, sftp=sftp)
    if reference_weighting_type is not None:
        connectomes_to_excel(matrix_refweighted, index_to_struct, connectome_xlsxpath_ref, overwrite=overwrite, verbose=verbose, sftp=sftp)
    if volume_weighting and reference_weighting_type is not None:
        connectomes_to_excel(matrix_vol_refweighted, index_to_struct, connectome_xlsxpath_volref, overwrite=overwrite, verbose=verbose, sftp=sftp)

    grouping_to_excel(matrix_sl, index_to_struct, grouping_xlsxpath, overwrite=overwrite, verbose=verbose, sftp=sftp)
    if verbose:
        txt = ("The excelfile was saved at "+grouping_xlsxpath)
        send_mail(txt, subject="Excel save")
        print(txt)
    
    del(matrix, grouping,pruned_streamlines_SL)


def tract_connectome_analysis_pickle(dwipath, trkpath, str_identifier, outpath, subject, ROI_excel, bvec_orient, masktype="T1",
                              inclusive=False, function_processes=1, forcestart=False, picklesave=True,
                              verbose=None):

    picklepath_connect = os.path.join(outpath, subject + str_identifier + '_connectomes.p')
    connectome_xlsxpath = os.path.join(outpath, subject + str_identifier + "_connectomes.xlsx")
    grouping_xlsxpath = os.path.join(outpath, subject + str_identifier + "_grouping.xlsx")
    picklepath_grouping = os.path.join(outpath, subject + str_identifier + '_grouping.p')
    picklepath_connect_test = os.path.join(outpath, subject + str_identifier + '_connectomes_test.p')
    picklepath_temp = os.path.join(outpath, subject + str_identifier + '_index_struct.p')
    mkcdir(outpath)

    with open(picklepath_connect,'rb') as f:
        M = pickle.load(f)
    with open(picklepath_grouping,'rb') as f:
        grouping = pickle.load(f)
    with open(picklepath_connect_test,'rb') as f:
        M2 = pickle.load(f)
    with open(picklepath_temp,'rb') as f:
        index_to_struct = pickle.load(f)

    matrix_sl = np.empty(np.shape(M), dtype=object)
    for i in np.arange(np.shape(matrix_sl)[0]):
        for j in np.arange(np.shape(matrix_sl)[1]):
            matrix_sl[i, j] = []
    for key in grouping.keys():
        matrix_sl[key] = grouping[key]
        matrix_sl[tuple(np.flip(key))] = grouping[key]

    M = np.delete(M, 0, 0)
    M = np.delete(M, 0, 1)

    matrix_sl = np.delete(matrix_sl, 0, 0)
    matrix_sl = np.delete(matrix_sl, 0, 1)

    if os.path.exists(picklepath_connect) and os.path.exists(connectome_xlsxpath) and os.path.exists(grouping_xlsxpath) and forcestart:
        os.remove(picklepath_connect)
        os.remove(connectome_xlsxpath)
        os.remove(grouping_xlsxpath)

    connectomes_to_excel(M, index_to_struct, connectome_xlsxpath)
    grouping_to_excel(matrix_sl, index_to_struct, grouping_xlsxpath)
    if verbose:
        txt = ("The excelfile was saved at "+grouping_xlsxpath)
        send_mail(txt, subject="Excel save")
        print(txt)

def tract_connectome_analysis_old(dwipath, trkpath, str_identifier, outpath, subject, whitematter_labels, targetrois, labelslist, ROI_excel, bvec_orient=[1,2,3], verbose=None):

    #str_identifier = roistring + tractsize + '_stepsize_' + str(stepsize)
    trkfile = gettrkpath(trkpath, subject, str_identifier, verbose)
    fa_data, _, gtab, vox_size, hdr, header = getfa(dwipath, subject, bvec_orient, verbose)
    labelmask, _ = getlabelmask(dwipath, subject, bvec_orient, verbose)

    import numpy as np
    whitemask = np.zeros(np.shape(labelmask),dtype=int)
    for label in whitematter_labels:
        whitemask = whitemask + (labelmask == label)

    roimask = np.zeros(np.shape(labelmask),dtype=int)
    for label in labelslist:
        roimask = roimask + (labelmask == label)


    """
    white_matter = binary_dilation(whitemask)
    
    csamodel = shm.CsaOdfModel(gtab, 6)
    csapeaks = peaks.peaks_from_model(model=csamodel,
                                      data=fdwi_data,
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
    """

    """
    trkdata = load_trk(trkfile, "same")
    affine = trkdata._affine
    trkdata.to_vox()
    trkstreamlines = trkdata.streamlines

    fullstreamlines = Streamlines(trkstreamlines)
    """

    trkprunepath = trkpath + '/' + subject + str_identifier + '_pruned.trk'
    #trkpaths = glob.glob(trkpath + '/' + subject + '_' + tractsize + strproperty + 'stepsize_' + str(stepsize) + '.trk')

    prunesave = True
    pruneforcestart = False
    if not os.path.isfile(trkprunepath) or pruneforcestart:

        trkdata = load_trk(trkfile, "same")
        print(trkfile)
        affine = trkdata._affine
        trkdata.to_vox()
        trkstreamlines = trkdata.streamlines

        cutoff=4
        pruned_streamlines = prune_streamlines(list(trkstreamlines), labelmask, cutoff=cutoff, verbose=verbose)
        pruned_streamlines_SL = Streamlines(pruned_streamlines)
        if hasattr(trkdata,'space_attribute'):
            header = trkdata.space_attribute
        elif hasattr(trkdata,'space_attributes'):
            header = trkdata.space_attributes
        myheader = create_tractogram_header(trkprunepath, *header)
        prune_sl = lambda: (s for s in pruned_streamlines)
        if prunesave:
            tract_save.save_trk_heavy_duty(trkprunepath, streamlines=prune_sl, affine=affine, header=myheader)
        del(prune_sl,pruned_streamlines,trkdata)
    else:
        trkprunedata = load_trk(trkprunepath, "same")
        trkprunedata.to_vox()
        pruned_streamlines_SL = trkprunedata.streamlines
        del(trkprunedata)

    affine_streams = np.eye(4)
    """
    #cc_slice = roimask == 1
    #cc_slice = roilabels
    #cc_streamlines = utils.target(trkstreamlines, affine_streams, cc_slice)
    #cc_streamlines = Streamlines(cc_streamlines)
    """

    M, grouping = utils.connectivity_matrix(pruned_streamlines_SL, affine_streams, labelmask,
                                            return_mapping=True,
                                            mapping_as_streamlines=True)

    print("The nunmber of tracts associated with the label 0 for subject " + subject+" is " + str(np.sum(M[0,:])))
    M = np.delete(M, 0, 0)
    M = np.delete(M, 0, 1)

    picklepath_connect = outpath + subject + "_" + str_identifier + '_connectomes.p'
    picklepath_grouping = outpath + subject + str_identifier + '_grouping.p'
    pickle.dump(M, open(picklepath_connect,"wb"))
    #pickle.dump(grouping, open(picklepath_grouping,"wb"))

    if verbose:
        txt= ("The connectomes were saved at "+picklepath_connect)
        send_mail(txt, subject="Pickle save")
        print(txt)

    excel_path = outpath + subject + "_" + str_identifier + "_connectomes.xlsx"
    connectomes_to_excel(M, ROI_excel, excel_path)
    if verbose:
        txt= ("The excelfile was saved at "+excel_path)
        send_mail(txt, subject="Excel save")
        print(txt)

    #whitem_slice = whitemask == 1
    #white_streamlines = utils.target(trkstreamlines, affine, whitem_slice)
    #white_streamlines = Streamlines(white_streamlines)
    """
    cc_slice = roimask == 1
    #cc_slice = roilabels
    cc_streamlines = utils.target(trkstreamlines, affine_streams, cc_slice)
    cc_streamlines = Streamlines(cc_streamlines)

    other_streamlines = utils.target(trkstreamlines, affine_streams, cc_slice,
                                     include=False)
    other_streamlines = Streamlines(other_streamlines)
    assert len(other_streamlines) + len(cc_streamlines) == len(trkstreamlines)

    from dipy.viz import window, actor, colormap as cmap

    # Enables/disables interactive visualization
    interactive = False

    # Make display objects
    color = cmap.line_colors(cc_streamlines)
    cc_streamlines_actor = actor.line(cc_streamlines,
                                      cmap.line_colors(cc_streamlines))
    cc_ROI_actor = actor.contour_from_roi(cc_slice, color=(1., 1., 0.),
                                          opacity=0.5)

    #vol_actor = actor.slicer(t1_data)
    vol_actor = actor.slicer(fa_data)

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
    targetrois = targetrois[0]
    window.record(r, n_frames=1, out_path=outpath + subject + targetrois + '_axial.png',
                  size=(800, 800))
    if interactive:
        window.show(r)
    r.set_camera(position=[-1, 0, 0], focal_point=[0, 0, 0], view_up=[0, 0, 1])
    window.record(r, n_frames=1, out_path=outpath + subject + targetrois + '_sagittal.png',
                  size=(800, 800))
    if interactive:
        window.show(r)

    atlas_legends = "/Users/alex/jacques/atlases/CHASSSYMM3AtlasLegends.xlsx"

    M, grouping = utils.connectivity_matrix(cc_streamlines, affine_streams, labelmask,
                                            return_mapping=True,
                                            mapping_as_streamlines=True)

    picklepath_connect = outpath + subject + "_" + tractsize + '_connectomes.p'
    picklepath_grouping = outpath + subject + tractsize + '_grouping.p'
    pickle.dump(M, open(picklepath_connect,"wb"))
    pickle.dump(grouping, open(picklepath_grouping,"wb"))

    M[:3, :] = 0
    M[:, :3] = 0

    import numpy as np
    import matplotlib.pyplot as plt
    plt.imshow(np.log1p(M), interpolation='nearest')
    plt.savefig(outpath + subject + "connectivity.png")

    lr_superiorfrontal_track = grouping[11, 54]
    shape = labelmask.shape
    dm = utils.density_map(lr_superiorfrontal_track, affine_streams, shape)

    import nibabel as nib
    from dipy.io.stateful_tractogram import Space, StatefulTractogram
    from dipy.io.streamline import save_trk

    # Save density map
    dm_img = nib.Nifti1Image(dm.astype("int16"), affine_streams)
    dm_img.to_filename(outpath + subject + "lr-superiorfrontal-dm.nii.gz")

    lr_sf_trk = Streamlines(lr_superiorfrontal_track)

    # Save streamlines
    sft = StatefulTractogram(lr_sf_trk, dm_img, Space.VOX)
    save_trk(sft, outpath + subject + "lr-superiorfrontal.trk")
    """


def ROI_labels_mask(diff_data, labelsmask, labelslist):

    mask = np.zeros(np.shape(labelsmask), dtype=int)

    for label in labelslist:
        mask = mask + (labelsmask == label)
    diff_data_masked = diff_data * np.repeat(mask[:, :, :, None], np.shape(diff_data)[3], axis=3)

    return(diff_data_masked, mask)

def get_diffusionattributes(diffpath, outpath, subject, str_identifier, vol_b0, ratio, bvec_orient,
                                createmask, overwrite, verbose):

    diff_data, affine, gtab, vox_size, fdiffpath, header, ref_info = getdiffdata_all(diffpath, subject, bvec_orient, verbose)

    if createmask:         # Build Brain Mask
        #changes to the outpath reading, to check next time
        mask, _ = diff_to_mask(diff_data, subject, affine, outpath, makefig=False, vol_idx=vol_b0, median_radius=5,
                              numpass=6, dilate=2)

    averages = msdki.mean_signal_bvalue(diff_data, gtab, bmag=None)
    b0 = averages[0]
    diff = averages[1]
    b0path = os.path.join(outpath, subject + "_b0.nii.gz")
    diffpath = os.path.join(outpath, subject + "_diff.nii.gz")
    if not os.path.exists(b0path):
        save_nifti(b0path, b0, affine)
    if not os.path.exists(diffpath):
        save_nifti(diffpath, diff, affine)

    """
    msdki_model = msdki.MeanDiffusionKurtosisModel(gtab)
    msdki_fit = msdki_model.fit(dwi_data, mask=mask)

    MSD = msdki_fit.msd
    MSK = msdki_fit.msk

    dki_model = dki.DiffusionKurtosisModel(gtab)
    dki_fit = dki_model.fit(dwi_data, mask=mask)

    MD = dki_fit.md
    MK = dki_fit.mk(0, 3)
    """

def diff_preprocessing(diffpath,outpath,subject, bvec_orient, denoise="none",savefa="yes",processes=1, createmask = True, vol_b0 = None, verbose = False):

    diff_fpath = getdiffpath(diffpath, subject, verbose)

    """
    subjfolder = glob.glob(os.path.join(datapath, "*" + identifier + "*"))[0]
    subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
    for bxhfile in subjbxh:
        bxhtype = checkbxh(bxhfile, False)
        if bxhtype == "dwi":
            dwipath = bxhfile.replace(".bxh", ".nii.gz")
            break
    """

    extractbvals(diff_fpath, subject)
    fbvals, fbvecs = move_bvals(diffpath, subject, outpath)
    gtab = getgtab(outpath, subject, bvec_orient)

    diff_data, affine, vox_size, diff_fpath, header, ref_info = getdiffdata(diff_fpath, subject, verbose)

    if verbose:
        print('Running the preprocessing for subject ' + subject)

    if createmask:         # Build Brain Mask
        mask, _ = diff_to_mask(diff_data, subject, affine, outpath, makefig=False, vol_idx=vol_b0, median_radius=5, numpass=6, dilate=2)
    else:
        mask, _ = getlabelmask(diffpath, subject, verbose=True)


    outpathdenoise = os.path.join(outpath, subject)
    diff_data, outpath_denoise = denoise_pick(diff_data, affine, header, outpathdenoise, mask, denoise, processes=processes, verbose=verbose) #accepts mpca, gibbs, all, none

    print(savefa)
    if savefa == "yes" or savefa == "y" or savefa == 1 or savefa is True or savefa == "all":
        outpathbmfa = make_tensorfit(diff_data,mask,gtab,affine,subject,outpath=outpath, verbose=verbose)
    else:
        print('FA was not calculated')
        outpathbmfa = None

    return (outpath_denoise)

import re


def check_dif_ratio(trkpath, subject, strproperty, ratio,sftp = None):

    string_list = ["_ratio_100", "_ratio_10", "_all"]
    strremoved = strproperty
    for word in string_list:
        strremoved = strremoved.replace(word, "_alt")

    for word in string_list:
        strnew = strremoved.replace("_alt", word)
        filepath = (os.path.join(trkpath, subject + "*" + strnew + '*.trk'))
        if sftp is None:
            trkpaths = glob.glob(filepath)
        else:
            trkpaths = glob_remote(filepath, sftp)
        if trkpaths:
            trkfile = trkpaths[0]
            trkname = os.path.basename(trkfile)
            trknewpath = os.path.join(trkpath, subject + strproperty + '.trk')
            ratiostart = trkname.find("ratio")
            if os.path.exists(trknewpath):
                return
            elif trkname.find("all") != -1:
                if ratio==1:
                    return
                elif ratio>1:
                    reducetractnumber(trkfile,trknewpath, getdata=False, ratio=ratio, return_affine=False, verbose=True)
                    return
                else:
                    raise Warning("Error with the input ratio")
            elif ratiostart != -1:
                temp = re.findall(r'\d+', trkname[ratiostart:])
                res = list(map(int, temp))
                newratio = ratio/res[0]
                if newratio>1:
                    reducetractnumber(trkfile,trknewpath, getdata=False, ratio=newratio, return_affine=False, verbose=True)
                return
    #trkfilepath = gettrkpath(trkpath, subject, strproperty, verbose)

def create_tracts(diffpath, outpath, subject, figspath, step_size, peak_processes, strproperty = "", ratio = 1, masktype="binary",
                      classifier="FA", labelslist=None, bvec_orient=[1,2,3], doprune=False, overwrite=False,
                  get_params=False, denoise="", verbose=None, sftp=None):

    check_dif_ratio(outpath, subject, strproperty, ratio)
    outpathtrk, trkexists = gettrkpath(outpath, subject, strproperty, pruned=doprune, verbose=False,sftp=sftp)

    if trkexists and overwrite is False:
        print("The tract creation of subject " + subject + " is already done")
        if get_params:
            subject, numtracts, minlength, maxlength, meanlength, stdlength, _, _, _ = \
                get_tract_params(outpathtrk, subject, strproperty, verbose = verbose)
            params = [numtracts, minlength, maxlength, meanlength, stdlength]
            return outpathtrk, None, params
        else:
            return outpathtrk, None, None

    if verbose:
        print('Running the ' + subject + ' file')

    mask, _ = getmask(diffpath, subject, masktype, verbose,sftp=sftp)

    diff_data, affine, gtab, vox_size, fdiffpath, header, ref_info = getdiffdata_all(diffpath, subject, bvec_orient, denoise=denoise, verbose=verbose,sftp=sftp)
    #fdiffpath = getdiffpath(diffpath, subject, denoise=denoise, verbose=verbose)

    if masktype == "dwi" and mask is None:
        warnings.warn(f'Did not find mask, assuming that the diffusion path {fdiffpath} is already masked')
        mask,_ = dwi_to_mask(diff_data, subject, affine, diffpath, masking='extract', makefig=False, header=header, verbose=True, sftp=sftp)

    if np.size(np.shape(mask)) == 1:
        mask = mask[0]
    if np.size(np.shape(mask)) == 4:
        mask = mask[:, :, :, 0]
    print("Mask shape is " + str(np.shape(mask)))

    if np.mean(diff_data) == 0:
        print("The subject " + subject + "could not be found at " + diffpath)
        return

    if classifier == "FA":
        outpathbmfa, mask = make_tensorfit(diff_data,mask,gtab,affine,subject,outpath=diffpath,verbose=verbose)

    print(verbose)
    if verbose:
        txt = ("The QCSA Tractmake is ready to launch for subject " + subject)
        print(txt)
        send_mail(txt,subject="QCSA main function start")
        print("email sent")

    if doprune:
        outpathtrk_noprune = os.path.join(outpath, subject + strproperty + '.trk')
        if os.path.exists(outpathtrk_noprune) and not os.path.exists(outpathtrk):
            cutoff = 2
            trkdata = load_trk_remote(outpathtrk_noprune, "same",sftp=sftp)
            affine = trkdata._affine
            trkdata.to_vox()
            trkstreamlines = trkdata.streamlines
            if np.size(np.shape(diff_data)) == 1:
                diff_data = diff_data[0]
            if np.size(np.shape(diff_data)) == 4:
                diff_data = diff_data[:, :, :, 0]
            print("Mask shape is " + str(np.shape(diff_data)))
            pruned_streamlines = prune_streamlines(list(trkstreamlines), mask, cutoff=cutoff, verbose=verbose)
            pruned_streamlines_SL = Streamlines(pruned_streamlines)
            if hasattr(trkdata, 'space_attribute'):
                ref_info = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                ref_info = trkdata.space_attributes
            myref = create_tractogram_header(outpathtrk, *ref_info)
            prune_sl = lambda: (s for s in pruned_streamlines)
            tract_save.save_trk_heavy_duty(outpathtrk, streamlines=prune_sl, affine=affine, header=myref)
            del (prune_sl, pruned_streamlines, trkdata)
            if get_params:
                numtracts, minlength, maxlength, meanlength, stdlength = get_tract_params(outpathtrk, subject,
                                                                                          strproperty,
                                                                                          verbose)
                params = [numtracts, minlength, maxlength, meanlength, stdlength]
                return subject, outpathtrk, params
            else:
                params = None
                return subject, outpathtrk, params

    outpathtrk, trkstreamlines, params = QCSA_tractmake(diff_data, affine, vox_size, gtab, mask, masktype, ref_info,
                                                        step_size, peak_processes, outpathtrk, subject, ratio,
                                                        overwrite, get_params, doprune, figspath=figspath,
                                                        verbose=verbose,sftp=sftp)

    if labelslist:
        print('In process of implementing')
        """
        labelslist= [59, 1059, 62, 1062]
        labelmask=mask
        roiname = "_hyptsept_"
        strproperty = roiname

        ratios = [1]
        roislist = [['fimbria'], ['corpus_callosum'], ['hypothalamus', 'septum'], ['primary_motor_cortex']]
        print("reached this spot")
        #save_roisubset(streamlines_generator, roislist, roisexcel, labelmask, stringstep, ratios, trkpath, subject, affine, header)

        from dipy.tracking.utils import length
        lengths = length(sg)
        del trkdata
        # lengths = list(length(trkstreamlines))
        lengths = list(lengths)
        numtracts = np.size(lengths)
        minlength = np.min(lengths)
        maxlength = np.max(lengths)
        meanlength = np.mean(lengths)
        stdlength = np.std(lengths)
        print("Numtracts is "+ numtracts + ", minlength is "+minlength+", maxlength is "+maxlength+", meanlength is "+meanlength+", stdlength is"+stdlength)
        """

    return subject, outpathtrk, params


def evaluate_coherence(dwipath,trkpath,subject,stepsize, tractsize, labelslist=None, outpathpickle=None,
                       outpathfig=None, processes=1, allsave=False, display=True, strproperty="", ratio=1,
                       verbose=None):

    _, affine, gtab, labelmask, vox_size, fdwipath, _, _ = getdiffdata_all(dwipath, subject)

    if isempty(labelslist):
        if labelmask is None:
            roimask = (diff_data[:, :, :, 0] > 0)
        else:
            roimask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at "+dwipath+ " for subject " + subject)
        roimask = np.zeros(np.shape(labelmask),dtype=int)
        for label in labelslist:
            roimask = roimask + (labelmask == label)
            
    if roimask.dtype == "bool":
        roimask = roimask.astype(int)
    """     
    if outpathpickle is None:
        outpathpickle=outpathtrk
    if outpathfig is None:
        outpathfig=outpathtrk
    #picklepath="/mnt/BIAC/munin3.dhe.duke.edu/Badea/Lab/mouse/C57_JS/Figures/N57433_roi.p"
    print(outpathpickle)
    print(subject)
    #outpathpickle=outpathpickle+"/"+subject
    #pickle.dump(roimask, open(outpicklepath+"_roimask.p", "wb"))
    #roimask = labels > 1

    #Changing the outpathfig to include the subject prefix so I dont carry that particular value around needlessly
    """
    trkpaths = glob.glob(trkpath + '/' + subject + '_' + tractsize + strproperty + 'stepsize_' + str(stepsize) + '.trk')
    trkfile = trkpaths[0]
    outpathtrk = trkpath + '/'    
    
    outpathfig = outpathfig+'/'
    print("Beginning Tract Evaluation of " + subject)
    stepsize = strfile(stepsize)
    doprune = True
    """
    if doprune:
        trkprunefile = trkpath + '/' + subject + '_' + tractsize + '_stepsize_' + stepsize + '_pruned.trk'
        if not os.path.exists(trkprunefile):

            trkdata = load_trk(trkfile, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkstreamlines = trkdata.streamlines
            trkstreamlines=prune_streamlines(list(trkstreamlines), fdwi_data[:, :, :, 0], cutoff=cutoff, verbose=verbose)
            myheader = create_tractogram_header(trkprunefile, *header)
            prune_sl = lambda: (s for s in trkstreamlines)
            tract_save.save_trk_heavy_duty(trkprunefile, streamlines=prune_sl,
                                           affine=affine, header=myheader)
    """

    trkdata = load_trk(trkfile, "same")
    trkdata.to_vox()
    trkstreamlines = trkdata.streamlines

    if verbose:
        txt="Beginning the evaluation of subject " + subject + "\n We will evaluate for one in " + str(ratio) \
            + " streamlines"
        print(txt)
        send_mail(txt,subject="LifE start msg ")

    kpath = (outpathpickle+"/Enhancekernel.p")

    from dipy.denoise.enhancement_kernel import EnhancementKernel
    D33 = 1.0
    D44 = 0.02
    t = 1
    k = EnhancementKernel(D33, D44, t)

    k = bundle_coherence(trkstreamlines,affine, k, t1_data=None,interactive=False)


def tract_getroi(trkstreamlines, affine, myheader, labelslist, labelmask, trkroipath, verbose = False):

    roimask = np.zeros(np.shape(labelmask), dtype=int)
    for label in labelslist:
        roimask = roimask + (labelmask == label)

    affinetemp=np.eye(4)
    trkroistreamlines = target(trkstreamlines, affinetemp, roimask, include=True, strict="longstring")
    trkroistreamlines = Streamlines(trkroistreamlines)
    #trkroipath = trkpath + '/' + subject + '_' + tractsize + strproperty + stepsize + '.trk'
    #myheader = create_tractogram_header(trkroipath, *header)
    roi_sl = lambda: (s for s in trkroistreamlines)
    tract_save.save_trk_heavy_duty(trkroipath, streamlines=roi_sl,
                affine=affine, header=myheader)
    if verbose:
        txt = "Successfully saved trk roi at "+trkroipath
    return trkroistreamlines, roi_sl

    
def evaluate_tracts(diffpath,trkpath,subject,stepsize, tractsize, labelslist=None, outpathpickle=None, outpathfig=None, processes=1, allsave=False, display=True, strproperty = "", verbose=None):

    """
    fdwi = dwipath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    dwidata, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    ffalabels = dwipath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'
    labels, affine_labels = load_nifti(ffalabels)

    roimask = (labels == 163) + (labels == 1163) # + (labels == 120) + (labels == 1120)
    
    fbvals = dwipath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'
    fbvecs = dwipath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'
    try:
        bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)
    except FileNotFoundError:
        raise Error("Missing bvals and/or bvecs")

    gtab = gradient_table(bvals, bvecs)

    """
    if tractsize == "small":
        ratio=100
    elif tractsize == "all":
        ratio=1
    elif tractsize == "tiny":
        ratio=1000

    diff_data, affine, gtab, labelmask, vox_size, fdiffpath, _, ref_info = getdiffdata_all(diffpath, subject)
    if isempty(labelslist):
        if labelmask is None:
            roimask = (diff_data[:, :, :, 0] > 0)
        else:
            roimask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at "+diffpath+ " for subject " + subject)
        roimask = np.zeros(np.shape(labelmask),dtype=int)
        for label in labelslist:
            roimask = roimask + (labelmask == label)

    if roimask.dtype == "bool":
        roimask=roimask.astype(int)

    save_nifti(diffpath+"/"+subject+"_roimask.nii.gz",roimask,affine)
    if outpathpickle is None:
        outpathpickle=outpathtrk
    if outpathfig is None:
        outpathfig=outpathtrk
    #picklepath="/mnt/BIAC/munin3.dhe.duke.edu/Badea/Lab/mouse/C57_JS/Figures/N57433_roi.p"
    print(outpathpickle)
    print(subject)
    #outpathpickle=outpathpickle+"/"+subject
    #pickle.dump(roimask, open(outpicklepath+"_roimask.p", "wb"))
    #roimask = labels > 1

    #Changing the outpathfig to include the subject prefix so I dont carry that particular value around needlessly
    outpathfig = outpathfig+'/'
    print("Beginning Tract Evaluation of " + subject)
    stepsize = strfile(stepsize)

    if verbose:
        txt="Beginning the evaluation of subject " + subject + "\n We will evaluate for one in " + str(ratio) \
            + " streamlines"
        print(txt)
        send_mail(txt,subject="LifE start msg ")

    trkpaths = glob.glob(trkpath + '/' + subject + '*_' + tractsize + '_*' + stepsize + '.trk')
    trkfile = trkpaths[0]
    outpathtrk = trkpath + '/'

    doprune=True
    cutoff = 2
    if doprune:
        trkprunefile = trkpath + '/' + subject + '_' + tractsize + '_stepsize_' + stepsize + '_pruned.trk'
        if not os.path.exists(trkprunefile):

            trkdata = load_trk(trkfile, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkstreamlines = trkdata.streamlines
            trkstreamlines=prune_streamlines(list(trkstreamlines), diff_data[:, :, :, 0], cutoff=cutoff, verbose=verbose)
            myheader = create_tractogram_header(trkprunefile, *header)
            prune_sl = lambda: (s for s in trkstreamlines)
            tract_save.save_trk_heavy_duty(trkprunefile, streamlines=prune_sl,
                                           affine=affine, header=myheader)

        trkfile = trkprunefile

    if ratio != 1:
        trkminipath = trkpath + '/' + subject + '_' + tractsize + '_brain_ratio_' + str(ratio) + '_stepsize_' + stepsize + '.trk'
        if not os.path.exists(trkminipath):
            #if os.path.exists(headerpath) and os.path.exists(streamlinespath):
            #    trkorigstreamlines=pickle.load(open(streamlinespath, "rb"))
            #    header=pickle.load(open(headerpath, "rb"))
            #else:
            trkdata = load_trk(trkfile, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkorigstreamlines = trkdata.streamlines
            ministream = []
            for idx, stream in enumerate(trkorigstreamlines):
                if (idx % ratio) == 0:
                    ministream.append(stream)
            myheader = create_tractogram_header(trkminipath, *header)
            ratioed_sl_gen = lambda: (s for s in ministream)
            if allsave:
                tract_save.save_trk_heavy_duty(trkminipath, streamlines=ratioed_sl_gen,
                                    affine=affine, header=myheader)
        elif isempty(labelslist):
            trkdata = load_trk(trkminipath, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkstreamlines = trkdata.streamlines

    if not isempty(labelslist):
        trkroipath = trkpath + '/' + subject + '_' + tractsize + strproperty + "_stepsize_" + stepsize + '.trk'
        if not os.path.exists(trkroipath):
            if not 'trkorigstreamlines' in locals():
                trkdata = load_trk(trkfile, 'same')
                trkdata.to_vox()
                if hasattr(trkdata, 'space_attribute'):
                    header = trkdata.space_attribute
                elif hasattr(trkdata, 'space_attributes'):
                    header = trkdata.space_attributes
                trkorigstreamlines = trkdata.streamlines

            affinetemp=np.eye(4)
            trkstreamlines = target(trkorigstreamlines, affinetemp, roimask, include=True, strict="longstring")
            del(trkorigstreamlines)
            trkstreamlines = Streamlines(trkstreamlines)
            trkroipath = trkpath + '/' + subject + '_' + tractsize + strproperty + stepsize + '.trk'
            myheader = create_tractogram_header(trkroipath, *header)
            roi_sl = lambda: (s for s in trkstreamlines)
            if allsave:
                tract_save.save_trk_heavy_duty(trkroipath, streamlines=roi_sl,
                            affine=affine, header=myheader)
        else:
            trkdata = load_trk(trkroipath, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkstreamlines = trkdata.streamlines

        if ratio != 1:
            trkroiminipath = trkpath + '/' + subject + '_' + tractsize + strproperty + "ratio_" + str(ratio) + '_stepsize_' + stepsize + '.trk'
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
    
    else:
        if not 'trkstreamlines' in locals():
            trkdata = load_trk(trkfile, 'same')
            trkdata.to_vox()
            if hasattr(trkdata, 'space_attribute'):
                header = trkdata.space_attribute
            elif hasattr(trkdata, 'space_attributes'):
                header = trkdata.space_attributes
            trkstreamlines = trkdata.streamlines
    
    if display:
        saveconnectivitymatrix=True

        from dipy.viz import window, actor, colormap as cmap

        # Enables/disables interactive visualization
        interactive = False

        # Make display objects
        if outpathfig is not None:
            color = cmap.line_colors(trkstreamlines)
            streamlines_actor = actor.line(trkstreamlines,
                                              cmap.line_colors(trkstreamlines))
            ROI_actor = actor.contour_from_roi(roimask, color=(1., 1., 0.),
                                                  opacity=0.5)
            labels_actor = actor.contour_from_roi(labelmask, color=(3., 3., 30.),
                                                  opacity=0.2)
            vol_actor = actor.slicer(labelmask)

            vol_actor.display(x=40)
            vol_actor2 = vol_actor.copy()
            vol_actor2.display(z=35)

            # Add display objects to canvas
            r = window.Renderer()
            r.add(vol_actor)
            r.add(vol_actor2)
            r.add(streamlines_actor)
            r.add(ROI_actor)
            r.add(labels_actor)

            # Save figures
            window.record(r, n_frames=1, out_path=outpathfig+ subject+'_fimbria_axial.png',
                          size=(800, 800))
            if verbose:
                print("saving figure to: "+outpathfig+ subject+'_fimbria_axial.png')
            interactive=False
            if interactive:
                window.show(r)
            r.set_camera(position=[-1, -1, -1], focal_point=[0, 0, 0], view_up=[0, 0, 1])
            window.record(r, n_frames=1, out_path=outpathfig+ subject + '_fimbria_sagittal.png',
                          size=(800, 800))
        #if interactive:
        #    window.show(r)
        """

        if saveconnectivitymatrix:
            M, grouping = utils.connectivity_matrix(trkstreamlines, affine,
                                                    roimask.astype(np.uint8),
                                                    return_mapping=True,
                                                    mapping_as_streamlines=True)
            M[:3, :] = 0
            M[:, :3] = 0
            plt.imshow(np.log1p(M), interpolation='nearest')
            if outpathfig is not None:
                plt.savefig(outpathfig+"connectivity.png")
        """

    #if display:
    #    window_show_test(trkstreamlines, mask_roi, anat, interactive= True, outpath=None)


    #trk_streamlines = [s[0] for s in nib.trackvis.read(trkfile, points_space='voxel')[0]]
    #len_sl=len(trk_streamlines)
    #% memit fiber_fit = fiber_model.fit(data, trk_streamlines[2 ** 12], affine=np.eye(4))
    display = False
    duration=time()
    strproperty = tractsize + strproperty + "ratio_" + str(ratio) + "_stepsize_" +stepsize
    model_error, mean_error = LiFEvaluation(diff_data, trkstreamlines, gtab, subject=subject, header=header,
                                            roimask=roimask, affine=affine,display=display, outpathpickle=outpathpickle,
                                            outpathtrk=outpathtrk, processes=processes, outpathfig=outpathfig,
                                            strproperty=strproperty, verbose=verbose)
    #picklepath = '/Users/alex/jacques/test_pickle_subj'+subject+'.p'
    #results=[outpathtrk,model_error,mean_error]
    #if subject == "N54859":
    #    pickle.dump(results, open(picklepath, "wb"))

    if verbose:
        txt=("Finished life evaluation of subject " + (subject)+ ", whole process took " + str(time()-duration) + " s")
        print(txt)
        send_mail(txt,subject="LifE save msg ")
    #picklepath = trkpath+subject+'lifevals.p'
    #pickle.dump(tracteval_results, open(picklepath,"wb"))
    
    return [outpathfig, model_error, mean_error]
