#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:39:41 2020

@author: Jacques Stout
"""

import numpy as np
import os
from time import time
from dipy.io.image import save_nifti
from dipy.reconst.shm import CsaOdfModel
from dipy.direction import peaks_from_model
from dipy.tracking.local_tracking import LocalTracking
from dipy.direction import peaks
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
from dipy.tracking import utils
from BIAC_tools import send_mail, isempty
from tract_save import save_trk_heavy_duty
from dipy.io.utils import create_tractogram_header
from dipy.io.streamline import load_trk
import tract_save
from tract_handler import get_trk_params, get_tract_params
import glob
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion, ThresholdStoppingCriterion

#from dipy.denoise.localpca import mppca
#import dipy.tracking.life as life
#import dipy.reconst.dti as dti
#from dipy.reconst.dti import fractional_anisotropy

from dipy.io.image import load_nifti
import matplotlib.pyplot as plt

from tract_handler import target, prune_streamlines
import nibabel as nib
from dipy.tracking.streamline import Streamlines

def make_tensorfit(data,mask,gtab,affine,subject,outpath,strproperty, overwrite="yes", forcestart = False, verbose=None):


    from dipy.reconst.dti import TensorModel
    if os.path.isdir(outpath):
        outpathbmfa = outpath + '/bmfa' + subject + '.nii.gz'
    elif os.path.isfile(outpath):
        outpathbmfa = os.path.dirname(outpath) + '/bmfa' + subject + '.nii.gz'
    if os.path.isfile(outpathbmfa) and overwrite == 'no':
        print("FA has already been saved at " + outpathbmfa + ", no overwrite")
        return outpathbmfa
    if verbose:
        print('Calculating the tensor model from bval/bvec values of ', subject)

    outpathbmfa = os.path.join(outpath, subject + '_bmfa.nii.gz')
    if os.path.exists(outpathbmfa) and not forcestart:
        if verbose:
            txt = "FA already computed at "+outpathbmfa
            print(txt)
        return(outpathbmfa)
    else:
        tensor_model = TensorModel(gtab)

        t1 = time()
        if len(mask.shape) == 4:
            mask = mask[...,0]
        tensor_fit = tensor_model.fit(data, mask)

        duration1 = time() - t1
        if verbose:
            print(subject + ' DTI duration %.3f' % (duration1,))

        outpathbmfa = outpath + 'bmfa' + subject + strproperty + '.nii.gz'
        print(outpathbmfa)
        save_nifti(outpathbmfa, tensor_fit.fa, affine)
        if verbose:
            print('Saving subject'+ subject+ ' at ' + outpathbmfa)
            #print('Warning: Failure to save the fa as nifti')
            #outpathbmfa = None

        #fa = tensor_fit.fa
        return outpathbmfa

# strproperty, trkpath, subject, affine, header
#(trkfile, roislist, roisexcel, labelmask, ratio)

def save_roisubset(trkfile, roislist, roisexcel, labelmask):

    trkdata = load_trk(trkfile, 'same')
    trkdata.to_vox()
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes
    trkstreamlines = trkdata.streamlines
    import pandas as pd
    df = pd.read_excel(roisexcel, sheet_name='Sheet1')
    df['Structure'] = df['Structure'].str.lower()

    for rois in roislist:

        labelslist = []
        for roi in rois:
            rslt_df = df.loc[df['Structure'] == roi.lower()]
            if rois[0].lower() == "wholebrain" or rois[0].lower() == "brain":
                labelslist = None
            else:
                labelslist = np.concatenate((labelslist, np.array(rslt_df.index2)))
        print(labelslist)
        if isempty(labelslist) and roi.lower() != "wholebrain" and roi.lower() != "brain":
            txt = "Warning: Unrecognized roi, will take whole brain as ROI. The roi specified was: " + roi
            print(txt)

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

        trkroipath = trkfile.replace(".trk", "_" + rois + ".trk")
        if not os.path.exists(trkroipath):
            affinetemp = np.eye(4)
            trkroistreamlines = target(trkstreamlines, affinetemp, roimask, include=True, strict="longstring")
            trkroistreamlines = Streamlines(trkroistreamlines)
            myheader = create_tractogram_header(trkroipath, *header)
            roi_sl = lambda: (s for s in trkroistreamlines)
            tract_save.save_trk_heavy_duty(trkroipath, streamlines=roi_sl,
                                           affine=header[0], header=myheader)


def QCSA_tractmake(data, affine, vox_size, gtab, mask, header, step_size, peak_processes, outpathtrk, subject='NA',
                   ratio=1, overwrite=False, get_params=False, doprune = False, pathfa = None, verbose=None):
    # Compute odfs in Brain Mask
    t2 = time()

    csa_model = CsaOdfModel(gtab, 6)
    if peak_processes < 2:
        parallel = False
    else:
        parallel = True
    if verbose:
        send_mail("Starting calculation of Constant solid angle model for subject " + subject,subject="CSA model start")

    wholemask = np.where(mask == 0, False, True)

    csa_peaks = peaks_from_model(model=csa_model,
                                 data=data,
                                 sphere=peaks.default_sphere,  # issue with complete sphere
                                 mask=wholemask,
                                 relative_peak_threshold=.5,
                                 min_separation_angle=25,
                                 parallel=parallel,
                                 nbr_processes=peak_processes)

    duration = time() - t2
    if verbose:
        print(subject + ' CSA duration %.3f' % (duration,))

    t3 = time()


    if verbose:
        send_mail('Computing classifier for local tracking for subject ' + subject +
                  ',it has been ' + str(round(duration)) + 'seconds since the start of tractmaker',subject="Seed computation" )

        print('Computing classifier for local tracking for subject ' + subject)

    if pathfa is None:
        classifier = BinaryStoppingCriterion(wholemask)
    else:
        #tensor_model = dti.TensorModel(gtab)
        #tenfit = tensor_model.fit(data, mask=labels > 0)
        #FA = fractional_anisotropy(tenfit.evals)
        FA_nifti = load_nifti(pathfa)
        FA = FA_nifti[0]
        classifier = ThresholdStoppingCriterion(FA, .2)

        fig = plt.figure()
        mask_fa = FA.copy()
        mask_fa[mask_fa < 0.2] = 0
        plt.xticks([])
        plt.yticks([])
        plt.imshow(mask_fa[:, :, data.shape[2] // 2].T, cmap='gray', origin='lower',
                   interpolation='nearest')
        fig.tight_layout()
        fig.savefig('/Users/alex/whiston_figures/whiston_methodtesting/threshold_fa.png')

    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    # seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))

    if verbose:
        print('Computing seeds')
    seeds = utils.seeds_from_mask(wholemask, density=1,
                                  affine=np.eye(4))

    ##streamlines_generator = local_tracking.local_tracker(csa_peaks,classifier,seeds,affine=np.eye(4),step_size=.5)
    if verbose:
        print('Computing the local tracking')
        duration = time() - t2
        send_mail('Start of the local tracking ' + ',it has been ' + str(round(duration)) +
                  'seconds since the start of tractmaker', subject="Seed computation")

    stringstep = str(step_size)
    stringstep = stringstep.replace(".", "_")
    print("stringstep is "+stringstep)



    classifier = BinaryStoppingCriterion(wholemask)

    """
    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=1)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=0.5)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)
    """

    cutoff = 2
    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_pruned.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=1)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_pruned.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=0.5)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_pruned.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)







    FA_nifti = load_nifti(pathfa)
    FA = FA_nifti[0]
    classifier = ThresholdStoppingCriterion(FA, .2)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=1)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=0.5)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)



    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_pruned_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=1)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_pruned_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=0.5)
    streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                              verbose=verbose)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_pruned_fa.trk"
    myheader = create_tractogram_header(outpathtrk, *header)
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)



    """
    if verbose:
        duration = time() - t2
        txt = 'About to save streamlines at ' + outpathtrk + ',it has been ' + str(round(duration)) + \
              'seconds since the start of tractmaker',
        send_mail(txt,subject="Tract saving" )

    cutoff = 2
    if doprune:
        streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                                  verbose=verbose)
        myheader = create_tractogram_header(outpathtrk, *header)
        sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=myheader,
                            shape=mask.shape, vox_size=vox_size)
    else:
        sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
        myheader = create_tractogram_header(outpathtrk, *header)
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=myheader,
                            shape=mask.shape, vox_size=vox_size)
    if verbose:
        duration = time() - t2
        txt = "Tract files were saved at "+outpathtrk + ',it has been ' + str(round(duration)) + \
              'seconds since the start of tractmaker'
        print(txt)
        send_mail(txt,subject="Tract saving" )

    # save everything - will generate a 20+ GBytes of data - hard to manipulate

    # possibly add parameter in csv file or other to decide whether to save large tractogram file
    # outpathfile=outpath+subject+"bmCSA_detr"+stringstep+".trk"
    # myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))
    duration3 = time() - t2
    if verbose:
        print(duration3)
        print(subject + ' Tracking duration %.3f' % (duration3,))
        send_mail("Finished file save at "+outpathtrk+" with tracking duration of " + str(duration3) + "seconds",
                  subject="file save update" )
"""
    if get_params:
        numtracts, minlength, maxlength, meanlength, stdlength = get_trk_params(streamlines_generator, verbose)
        params = [numtracts, minlength, maxlength, meanlength, stdlength]
        if verbose:
            print("For subject " + str(subject) + " the number of tracts is " + str(numtracts) + ", the minimum length is " +
                  str(minlength) + ", the maximum length is " + str(maxlength) + ", the mean length is " +
                  str(meanlength) + ", the std is " + str(stdlength))
    else:
        params = None
    return outpathtrk, streamlines_generator, params

def QCSA_tractmake(data, affine, vox_size, gtab, mask, header, step_size, peak_processes, outpathtrk, subject='NA',
                   ratio=1, overwrite=False, get_params=False, pathfig = None, doprune = False, pathfa = None, verbose=None):
    # Compute odfs in Brain Mask
    t2 = time()

    csa_model = CsaOdfModel(gtab, 6)
    if peak_processes < 2:
        parallel = False
    else:
        parallel = True
    if verbose:
        send_mail("Starting calculation of Constant solid angle model for subject " + subject,subject="CSA model start")

    wholemask = np.where(mask == 0, False, True)

    csa_peaks = peaks_from_model(model=csa_model,
                                 data=data,
                                 sphere=peaks.default_sphere,  # issue with complete sphere
                                 mask=wholemask,
                                 relative_peak_threshold=.5,
                                 min_separation_angle=25,
                                 parallel=parallel,
                                 nbr_processes=peak_processes)

    duration = time() - t2
    if verbose:
        print(subject + ' CSA duration %.3f' % (duration,))

    t3 = time()


    if verbose:
        send_mail('Computing classifier for local tracking for subject ' + subject +
                  ',it has been ' + str(round(duration)) + 'seconds since the start of tractmaker',subject="Seed computation" )

        print('Computing classifier for local tracking for subject ' + subject)

    if pathfa is None:
        classifier = BinaryStoppingCriterion(wholemask)
    else:
        #tensor_model = dti.TensorModel(gtab)
        #tenfit = tensor_model.fit(data, mask=labels > 0)
        #FA = fractional_anisotropy(tenfit.evals)
        FA_nifti = load_nifti(pathfa)
        FA = FA_nifti[0]
        classifier = ThresholdStoppingCriterion(FA, .2)

        if pathfig is not None:
            fig = plt.figure()
            mask_fa = FA.copy()
            mask_fa[mask_fa < 0.2] = 0
            plt.xticks([])
            plt.yticks([])
            plt.imshow(mask_fa[:, :, data.shape[2] // 2].T, cmap='gray', origin='lower',
                       interpolation='nearest')
            fig.tight_layout()
            fig.savefig(pathfig + 'threshold_fa.png')

    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    # seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))

    if verbose:
        print('Computing seeds')
    seeds = utils.seeds_from_mask(wholemask, density=1,
                                  affine=np.eye(4))

    #streamlines_generator = local_tracking.local_tracker(csa_peaks,classifier,seeds,affine=np.eye(4),step_size=step_size)
    if verbose:
        print('Computing the local tracking')
        duration = time() - t2
        send_mail('Start of the local tracking ' + ',it has been ' + str(round(duration)) +
                  'seconds since the start of tractmaker', subject="Seed computation")

    stringstep = str(step_size)
    stringstep = stringstep.replace(".", "_")
    print("stringstep is "+stringstep)

    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)

    if verbose:
        duration = time() - t2
        txt = 'About to save streamlines at ' + outpathtrk + ',it has been ' + str(round(duration)) + \
              'seconds since the start of tractmaker',
        send_mail(txt,subject="Tract saving" )

    cutoff = 2
    if doprune:
        streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                                  verbose=verbose)
        myheader = create_tractogram_header(outpathtrk, *header)
        sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=myheader,
                            shape=mask.shape, vox_size=vox_size)
    else:
        sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
        myheader = create_tractogram_header(outpathtrk, *header)
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=myheader,
                            shape=mask.shape, vox_size=vox_size)
    if verbose:
        duration = time() - t2
        txt = "Tract files were saved at "+outpathtrk + ',it has been ' + str(round(duration)) + \
              'seconds since the start of tractmaker'
        print(txt)
        send_mail(txt,subject="Tract saving" )

    # save everything - will generate a 20+ GBytes of data - hard to manipulate

    # possibly add parameter in csv file or other to decide whether to save large tractogram file
    # outpathfile=outpath+subject+"bmCSA_detr"+stringstep+".trk"
    # myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))
    duration3 = time() - t2
    if verbose:
        print(duration3)
        print(subject + ' Tracking duration %.3f' % (duration3,))
        send_mail("Finished file save at "+outpathtrk+" with tracking duration of " + str(duration3) + "seconds",
                  subject="file save update" )

    if get_params:
        numtracts, minlength, maxlength, meanlength, stdlength = get_trk_params(streamlines_generator, verbose)
        params = [numtracts, minlength, maxlength, meanlength, stdlength]
        if verbose:
            print("For subject " + str(subject) + " the number of tracts is " + str(numtracts) + ", the minimum length is " +
                  str(minlength) + ", the maximum length is " + str(maxlength) + ", the mean length is " +
                  str(meanlength) + ", the std is " + str(stdlength))
    else:
        params = None

    return outpathtrk, streamlines_generator, params


"""
streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=step_size)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=1)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=0.5)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)
"""

"""
cutoff = 2
streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=step_size)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_pruned.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=1)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_pruned.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=0.5)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_pruned.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

FA_nifti = load_nifti(pathfa)
FA = FA_nifti[0]
classifier = ThresholdStoppingCriterion(FA, .2)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=step_size)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=1)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=0.5)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

"""
"""
streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=step_size)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_2_ratio_100_wholebrain_pruned_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=1)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_1_ratio_100_wholebrain_pruned_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

streamlines_generator = LocalTracking(csa_peaks, classifier,
                                      seeds, affine=np.eye(4), step_size=0.5)
streamlines_generator = prune_streamlines(list(streamlines_generator), data[:, :, :, 0], cutoff=cutoff,
                                          verbose=verbose)
sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
outpathtrk = outpathtrk = "/Users/alex/whiston_figures/whiston_methodtesting/H26637_stepsize_0_5_ratio_100_wholebrain_pruned_fa.trk"
myheader = create_tractogram_header(outpathtrk, *header)
save_trk_heavy_duty(outpathtrk, streamlines=sg,
                    affine=affine, header=myheader,
                    shape=mask.shape, vox_size=vox_size)

"""