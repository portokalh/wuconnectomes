#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:48:38 2020

@author: Jacques Stout
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import nibabel as nib
import numpy as np
from dipy.io.streamline import load_trk
from os.path import splitext
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
import pickle
import pandas as pd
import smtplib

import os, re, sys, io, struct, socket, datetime
from email.mime.text import MIMEText
import glob

from dipy.tracking.utils import unique_rows


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

from dipy.segment.mask import segment_from_cfa
from dipy.segment.mask import bounding_box

import multiprocessing
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
import multiprocessing.pool


from dipy.reconst import shm

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
import xlsxwriter

from figures_handler import denoise_fig, shore_scalarmaps
from tract_eval import bundle_coherence, LiFEvaluation
from dif_to_trk import make_tensorfit, QCSA_tractmake
from BIAC_tools import send_mail, isempty
from tract_handler import target, prune_streamlines, get_trk_params, get_tract_params
from nifti_handler import getfa, getdwidata
import tract_save
import xlrd
import warnings


def string_inclusion(string_option,allowed_strings,option_name):
    "checks if option string is part of the allowed list of strings for that option"
    try:
        string_option=string_option.lower()
    except AttributeError:
        if string_option is None:
            #raise Warning(option_name + " stated as None value, option will not be implemented")
            print(option_name + " stated as None value, option will not be implemented")
        else:
            raise AttributeError('Unrecognized value for ' + option_name)

    if not any(string_option == x for x in allowed_strings):
        raise ValueError(string_option + " is an unrecognized string, please check your input for " + option_name)
    if string_option == "none":
        print(option_name + " stated as None value, option will not be implemented")

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


def reducetractnumber(oldtrkfile, newtrkfilepath, getdata=True, ratio=10, verbose=False):

    trkdata = load_trk(oldtrkfile, "same")
    if verbose:
        print("laoaded ")
    trkdata.to_vox()
    header = trkdata.space_attribute
    affine = trkdata._affine
    trkstreamlines = trkdata.streamlines

    ministream=[]
    for idx, stream in enumerate(trkstreamlines):
        if (idx % ratio) == 0:
            ministream.append(stream)
    del trkstreamlines
    myheader = create_tractogram_header(newtrkfilepath, *header)
    ratioed_sl = lambda: (s for s in ministream)
    tract_save.save_trk_heavy_duty(newtrkfilepath, streamlines=ratioed_sl,
                                   affine=affine, header=myheader)
    if verbose:
        print("The file " + oldtrkfile + "was reduced to one "+str(ratio)+"th of its size and saved to "+newtrkfilepath)

    if getdata:
        return(ministream)
    else:
        return


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

    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, _, header = getdwidata(dwipath, subject, bvec_orient)
    outpath_fig = outpath + subject + "SHORE_maps.png"
    shore_scalarmaps(fdwi_data, gtab, outpath_fig, verbose = verbose)

def dwiconnectome_analysis(dwipath,outpath,subject, whitematter_labels, targetrois, labelslist, bvec_orient=[1,2,3], verbose=None):

    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, header, hdr = getdwidata(dwipath, subject, bvec_orient, verbose)
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
    vol_actor = actor.slicer(fdwi_data[:,:,:,0])


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

#def gettrkpath(trkpath, subject, tractsize, strproperty, stepsize, verbose=False):
def gettrkpath(trkpath, subject, str_identifier, verbose=False):
    filepath=(trkpath + '/' + subject + "*" + str_identifier + '.trk')
    trkpaths = glob.glob(filepath)
    if trkpaths:
        trkfile = trkpaths[0]
        if verbose:
            print("Subject " + subject + " was found at " + trkfile)
    else:
        print("Could not find "+filepath)
        return
    return trkfile

def deprecation(message):
    warnings.warn(message, DeprecationWarning, stacklevel=2)

def getlabelmask(mypath, subject,verbose=None):


    labelsoption = glob.glob(mypath + '/' + subject + '/' + subject + '*labels.nii.gz')
    if np.size(labelsoption)>0:
        labelspath = labelsoption[0]
    elif os.path.exists(mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'):
        labelspath = mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'):
        labelspath = mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'):
        labelspath = mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labelspath = mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labelspath = mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labelspath = mypath + '/mask.nii.gz'
    elif os.path.exists(mypath + '/mask.nii'):
        labelspath = mypath + '/mask.nii'

    if 'labelspath' in locals():
        labels, affine_labels = load_nifti(labelspath)
        if verbose:
            print("Label mask taken from " + labelspath)
    else:
        print('mask not found')
        txt = ("Label mask taken from " + labelspath)
        deprecation(txt)

    return labels, affine_labels

def connectomes_to_excel(connectome,ROI_excel,output_path):

    df = pd.read_excel(ROI_excel, sheet_name='Sheet1')
    structure = df['Structure']

    workbook = xlsxwriter.Workbook(output_path)
    worksheet = workbook.add_worksheet()

    num = 1
    for struct in structure:
        worksheet.write(0, num, struct)
        worksheet.write(num, 0, struct)
        num += 1

    row=0
    for col, data in enumerate(connectome):
        worksheet.write_column(row+1, col+1, data)

    workbook.close()

    return

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


def tract_connectome_analysis(dwipath, trkpath, str_identifier, outpath, subject, ROI_excel, bvec_orient, verbose=None):

    trkfilepath = gettrkpath(trkpath, subject, str_identifier, verbose)
    trkprunepath = gettrkpath(trkpath, subject, str_identifier+"_pruned", verbose)
    labelmask, _ = getlabelmask(dwipath, subject, verbose)
    fa_data, _, vox_size, hdr, header = getfa(dwipath, subject, bvec_orient, verbose)

    import numpy as np
    prunesave = True
    pruneforcestart = False

    if (trkfilepath is not None and trkprunepath is None and prunesave) or pruneforcestart:

        trkdata = load_trk(trkfilepath, "same")
        print(trkfilepath)
        affine = trkdata._affine
        trkdata.to_vox()
        trkstreamlines = trkdata.streamlines
        trkprunepath = os.path.dirname(trkfilepath) + '/' + subject + str_identifier + '_pruned.trk'
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
    elif trkprunepath is not None:
        trkprunedata = load_trk(trkprunepath, "same")
        trkprunedata.to_vox()
        pruned_streamlines_SL = trkprunedata.streamlines
        del(trkprunedata)

    affine_streams = np.eye(4)

    M, grouping = utils.connectivity_matrix(pruned_streamlines_SL, affine_streams, labelmask,
                                            return_mapping=True,
                                            mapping_as_streamlines=True)

    print("The nunmber of tracts associated with the label 0 for subject " + subject+" is " + str(np.sum(M[0,:])))
    M = np.delete(M, 0, 0)
    M = np.delete(M, 0, 1)

    picklepath_connect = outpath + subject + str_identifier + '_connectomes.p'
    picklepath_grouping = outpath + subject + str_identifier + '_grouping.p'
    pickle.dump(M, open(picklepath_connect,"wb"))

    if verbose:
        txt= ("The connectomes were saved at "+picklepath_connect)
        send_mail(txt, subject="Pickle save")
        print(txt)

    excel_path = outpath + subject + str_identifier + "_connectomes.xlsx"
    connectomes_to_excel(M, ROI_excel, excel_path)
    if verbose:
        txt= ("The excelfile was saved at "+excel_path)
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




def denoise_pick(data,affine,hdr,outpath,mask,type_denoise='macenko', processes = 1, savedenoise= True, verbose=False, display=None):

    allowed_strings=['mpca','yes','all','gibbs','none', 'macenko']
    string_inclusion(type_denoise, allowed_strings, "type_denoise")

    if type_denoise == 'macenko' or 'mpca' or type_denoise == 'yes' or type_denoise == 'all':
        # data, snr = marcenko_denoise(data, False, verbose=verbose)
        t = time()
        denoised_arr, sigma = mppca(data, patch_radius=2, return_sigma=True, processes=processes, verbose=verbose)
        outpath_mpca = outpath + '_mpca.nii.gz'
        save_nifti(outpath_mpca, denoised_arr, affine, hdr=hdr)
        print("Saved image at " + outpath_mpca)

        mean_sigma = np.mean(sigma[mask])
        b0 = denoised_arr[..., 0]

        mean_signal = np.mean(b0[mask])

        snr = mean_signal / mean_sigma

        if verbose:
            print("Time taken for local MP-PCA ", -t +
                  time())
            print("The SNR of the b0 image appears to be at " + str(snr))
        if display:
            marcenko_denoise_fig(data, denoised_arr, type='macenko')

        data = denoised_arr

    if type_denoise == 'gibbs' or type_denoise =='all':
        t = time()
        data_corrected = gibbs_removal(data, slice_axis=2)
        outpath_gibbs = outpath + '_gibbs.nii.gz'
        save_nifti(outpath_gibbs, denoised_arr, affine, hdr=hdr)
        if verbose:
            print("Time taken for the gibbs removal " - t + time())
        if display:
            denoise_fig(data,data_corrected,type='gibbs')

        data=data_corrected
        
    if type_denoise == 'none':
        print('No denoising was done')

    return data


def dwi_preprocessing(dwipath,outpath,subject, bvec_orient, denoise="none",savefa="yes",processes=1, labelslist=None, strproperty="", verbose = False):

    bvec_orient=[1,2,3]
    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, header, hdr, _ = getdwidata(dwipath, subject, bvec_orient, verbose)
    #fdwi_data, affine, vox_size = load_nifti(fdwi, return_voxsize=True) (mypath, subject, bvec_orient=[1,2,3], verbose=None)

    #labels = glob.glob(mypath + '/' + subject + '*labels*nii.gz') #ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    #mask = glob.glob(mypath + '/' + subject + '*mask*nii.gz')
    #labelmask = labels + mask #Whether the image is a labels or a mask, we output a mask for mpca
        
    # Build Brain Mask
    print(labelslist)
    if isempty(labelslist):
        mask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at " + dwipath + " for subject " + subject)
        mask = np.zeros(np.shape(labelmask), dtype=int)
        for label in labelslist:
            mask = mask + (labelmask == label)
        fdwi_data = fdwi_data * np.repeat(mask[:, :, :, None], np.shape(fdwi_data)[3], axis=3)
    if verbose:
        print('Running the ' + subject + ' file')

    fdwi_data = fdwi_data * np.repeat(mask[:, :, :, None], np.shape(fdwi_data)[3], axis=3)

    outpathdenoise= outpath + subject + '_nii4D_RAS'
    print(denoise)
    fdwi_data = denoise_pick(fdwi_data, affine,hdr, outpathdenoise, mask, denoise, processes=processes, verbose=verbose) #accepts mpca, gibbs, all, none

    #testsnr => not yet fully implemented
    print(savefa)
    if savefa == "yes" or savefa == "y" or savefa == 1 or savefa is True or savefa == "all":
        outpathbmfa = make_tensorfit(fdwi_data,mask,gtab,affine,subject,outpath=outpath,strproperty=strproperty,verbose = verbose)
    else:
        print('FA was not calculated')
        outpathbmfa=None

def create_tracts(dwipath,outpath,subject,step_size,peak_processes,strproperty="",ratio=1,save_fa="yes",
                      labelslist = None, bvec_orient=[1,2,3], doprune=False, overwrite="no", get_params = False,
                  verbose=None):

    print("Do prune is "+str(doprune))
    if verbose:
        print('Running the ' + subject + ' file')

    #outpath_subject = outpathsubject + saved_str + '_stepsize_' + str(step_size) + '.trk'

    fdwi_data, affine, gtab, mask, vox_size, fdwipath, hdr, header = getdwidata(dwipath, subject, bvec_orient)
    if np.mean(fdwi_data) == 0:
        print("The subject " + subject + "could not be found at " + dwipath)
        return

    if save_fa == "yes" or save_fa == "y" or save_fa == 1 or save_fa is True or save_fa == "all" or save_fa == "only":
        make_tensorfit(fdwi_data,mask,gtab,affine,subject,outpath=dwipath,strproperty=strproperty,verbose=verbose)

    else:
        print('FA was not calculated')
        outpathbmfa = None

    print(verbose)
    if verbose:
        txt = ("The QCSA Tractmake is ready to launch for subject " + subject)
        print(txt)
        send_mail(txt,subject="QCSA main function start")
        print("email sent")

    outpathtrk, trkstreamlines, params = QCSA_tractmake(fdwi_data, affine, vox_size, gtab, mask, header, step_size,
                                                        peak_processes, outpath, subject, strproperty, ratio,
                                                        overwrite, get_params, doprune, verbose=verbose)

    if get_params is True and params is None:
        numtracts, minlength, maxlength, meanlength, stdlength = get_trk_params(trkstreamlines, verbose)
        params = [numtracts, minlength, maxlength, meanlength, stdlength]


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

    _, affine, gtab, labelmask, vox_size, fdwipath, _, _ = getdwidata(dwipath, subject)

    if isempty(labelslist):
        if labelmask is None:
            roimask = (fdwi_data[:, :, :, 0] > 0)
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

    
def evaluate_tracts(dwipath,trkpath,subject,stepsize, tractsize, labelslist=None, outpathpickle=None, outpathfig=None, processes=1, allsave=False, display=True, strproperty = "", verbose=None):

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

    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, _, header = getdwidata(dwipath, subject)
    if isempty(labelslist):
        if labelmask is None:
            roimask = (fdwi_data[:, :, :, 0] > 0)
        else:
            roimask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at "+dwipath+ " for subject " + subject)
        roimask = np.zeros(np.shape(labelmask),dtype=int)
        for label in labelslist:
            roimask = roimask + (labelmask == label)

    if roimask.dtype == "bool":
        roimask=roimask.astype(int)

    save_nifti(dwipath+"/"+subject+"_roimask.nii.gz",roimask,affine)
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
            trkstreamlines=prune_streamlines(list(trkstreamlines), fdwi_data[:, :, :, 0], cutoff=cutoff, verbose=verbose)
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
    model_error, mean_error = LiFEvaluation(fdwi_data, trkstreamlines, gtab, subject=subject, header=header,
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
