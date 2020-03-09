#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:39:41 2020

@author: Jacques Stout
"""

import numpy as np

from time import time
from dipy.io.image import save_nifti
from dipy.reconst.shm import CsaOdfModel
from dipy.direction import peaks_from_model
from dipy.tracking.local_tracking import LocalTracking
from dipy.direction import peaks
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
from dipy.tracking import utils
from BIAC_tools import send_mail
from tract_save import save_trk_heavy_duty

#from dipy.denoise.localpca import mppca
#import dipy.tracking.life as life

from tract_handler import target, prune_streamlines

def make_tensorfit(data,mask,gtab,affine,subject,outpath,strproperty,verbose=None):


    from dipy.reconst.dti import TensorModel

    if verbose:
        print('Calculating the tensor model from bval/bvec values of ', subject)
    tensor_model = TensorModel(gtab)

    t1 = time()
    print(np.shape(data))
    print(np.shape(mask))
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

def QCSA_tractmake(data,affine,vox_size,gtab,mask,trkheader,step_size,peak_processes,outpathsubject,ratio,verbose=None,subject = 'NA'):
    # Compute odfs in Brain Mask
    t2 = time()

    csa_model = CsaOdfModel(gtab, 6)
    if peak_processes < 2:
        parallel = False
    else:
        parallel = True
    if verbose:
        send_mail("Starting calculation of Constant solid angle model for subject " + subject,subject="CSA model start")
        #print("Starting calculation of Constant solid angle model for subject " + subject)
        #headers="From: %s\r\nTo: %s\r\nSubject:Reached the point where we start the CSA model for subject %s!\r\n\r\n" % (useremail,useremail,subject)
        #text="""The CSA model is about to be calculated""" % (warning_time,max_run_hours,os.getpid())  
        #message=headers+text 
        #mailServer=smtplib.SMTP(serverURL) 
        #mailServer.sendmail(useremail,useremail,message) 
        #mailServer.quit() 
    wholemask = np.where(mask == 0, False, True)
    csa_peaks = peaks_from_model(model=csa_model,
                                 data=data,
                                 sphere=peaks.default_sphere,  # issue with complete sphere
                                 mask=wholemask,
                                 relative_peak_threshold=.5,
                                 min_separation_angle=25,
                                 parallel=parallel,
                                 nbr_processes=peak_processes)

    duration2 = time() - t2
    if verbose:
        print(duration2) \

    if verbose:
        print(subject + ' CSA duration %.3f' % (duration2,))

    t3 = time()

    from dipy.tracking.stopping_criterion import BinaryStoppingCriterion

    if verbose:
        send_mail('Computing classifier for local tracking for subject ' + subject,subject="Seed computation" )
        print('Computing classifier for local tracking for subject ' + subject)
        #headers="From: %s\r\nTo: %s\r\nSubject:Seed computation\r\n\r\n" % (useremail,useremail)
        #text="""About to start binary stopping criterion, duration of CSA was %.2f""" % (duration2)  
        #message=headers+text 
        #mailServer=smtplib.SMTP(serverURL) 
        #mailServer.sendmail(useremail,useremail,message) 
        #mailServer.quit() 
    classifier = BinaryStoppingCriterion(wholemask)

    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    # seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))
    # why are those not binary?
    if verbose:
        print('Computing seeds')
    seeds = utils.seeds_from_mask(wholemask, density=1,
                                  affine=np.eye(4))

    ##streamlines_generator = local_tracking.local_tracker(csa_peaks,classifier,seeds,affine=np.eye(4),step_size=.5)
    if verbose:
        print('Computing the local tracking')
        send_mail('Start of the local tracking ', subject="Seed computation")
        #headers="From: %s\r\nTo: %s\r\nSubject:Reached the point where we start the local tracking!\r\n\r\n" % (useremail,useremail)
        #text="""Seeds have been computed, about to start the local tracking"""
        #message=headers+text 
        #mailServer=smtplib.SMTP(serverURL) 
        #mailServer.sendmail(useremail,useremail,message) 
        #mailServer.quit() 

    stringstep = str(step_size)
    stringstep = stringstep.replace(".", "_")
    # stringstep=""
    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)

    if verbose:
        #headers="From: %s\r\nTo: %s\r\nSubject:Reached the point where we start saving the file!\r\n\r\n" % (useremail,useremail)
        txt = 'About to save streamlines at ' + outpathsubject
        send_mail(txt,subject="Seed computation" )

    ratios=str(ratio)
    sg = lambda: (s for i, s in enumerate(streamlines_generator) if i % ratio == 0)
    outpathtrk = outpathsubject + "_bmCSA_detr_ratio"+ ratios +"_stepsize_" + stringstep + ".trk"
    save_trk_heavy_duty(outpathtrk, streamlines=sg,
                affine=affine, header=trkheader,
                shape=mask.shape, vox_size=vox_size)
    txt="Tract files were saved at " + outpathtrk

    # save everything - will generate a 20+ GBytes of data - hard to manipulate

    # possibly add parameter in csv file or other to decide whether to save large tractogram file
    # outpathfile=outpath+subject+"bmCSA_detr"+stringstep+".trk"
    # myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))

    doprune=True
    cutoff = 2
    trkpath = outpathsubject
    if doprune:
        trkprunefile = trkpath + '/' + subject + '_ratio_' + ratios + '_stepsize_' + stepsize + '_pruned.trk'
        trkstreamlines = sg
        trkstreamlines=prune_streamlines(list(trkstreamlines), fdwi_data[:, :, :, 0], cutoff=cutoff, verbose=verbose)
        myheader = create_tractogram_header(trkprunefile, *header)
        prune_sl = lambda: (s for s in trkstreamlines)
        tract_save.save_trk_heavy_duty(trkprunefile, streamlines=prune_sl,
                                           affine=affine, header=myheader)
        trkfile = trkprunefile

    labelslist= [59, 1059, 62, 1062]
    labelmask=mask
    roiname = "_hyptsept_"
    strproperty = roiname

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
    
    if not isempty(labelslist):
        trkroipath = trkpath + '/' + subject + "_" + roiname + "_stepsize_" + stepsize + '.trk'
        if not os.path.exists(trkroipath):
            affinetemp=np.eye(4)
            trkstreamlines = target(trkorigstreamlines, affinetemp, roimask, include=True, strict="longstring")
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
            trkroiminipath = trkpath + '/' + subject + '_ratio_' + ratios + roiname + "_stepsize_" + stepsize + '.trk'
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


    labelslist= [121, 1121]
    labelmask=mask
    roiname = "_corpus_callosum_"
    strproperty = roiname

    if isempty(labelslist):
        if labelmask is None:
            roimask = (fdwi_data[:, :, :, 0] > 0)
        else:
            roimask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at " + dwipath+ " for subject " + subject)
        roimask = np.zeros(np.shape(labelmask),dtype=int)
        for label in labelslist:
            roimask = roimask + (labelmask == label)
    
    if not isempty(labelslist):
        trkroipath = trkpath + '/' + subject + "_" + roiname + "_stepsize_" + stepsize + '.trk'
        if not os.path.exists(trkroipath):
            affinetemp=np.eye(4)
            trkstreamlines = target(trkorigstreamlines, affinetemp, roimask, include=True, strict="longstring")
            trkstreamlines = Streamlines(trkstreamlines)
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
            trkroiminipath = trkpath + '/' + subject + '_ratio_' + ratios + roiname + "_stepsize_" + stepsize + '.trk'
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

    labelslist= [19, 1019]
    labelmask=mask
    roiname = "_primary_motor_cortex_"
    strproperty = roiname

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
    
    if not isempty(labelslist):
        trkroipath = trkpath + '/' + subject + "_" + roiname + "_stepsize_" + stepsize + '.trk'
        if not os.path.exists(trkroipath):
            affinetemp=np.eye(4)
            trkstreamlines = target(trkorigstreamlines, affinetemp, roimask, include=True, strict="longstring")
            trkstreamlines = Streamlines(trkstreamlines)
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
            trkroiminipath = trkpath + '/' + subject + '_ratio_' + ratios + roiname + "_stepsize_" + stepsize + '.trk'
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

    labelslist= [120, 1120]
    labelmask=mask
    roiname = "_fimbria2_"
    strproperty = roiname

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
    
    if not isempty(labelslist):
        trkroipath = trkpath + '/' + subject + "_" + roiname + "_stepsize_" + stepsize + '.trk'
        if not os.path.exists(trkroipath):
            affinetemp=np.eye(4)
            trkstreamlines = target(trkorigstreamlines, affinetemp, roimask, include=True, strict="longstring")
            trkstreamlines = Streamlines(trkstreamlines)
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
            trkroiminipath = trkpath + '/' + subject + '_ratio_' + ratios + roiname + "_stepsize_" + stepsize + '.trk'
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


    # the function above will bring all streamlines in memory
    # streamlines = Streamlines(streamlines_generator)

    # save a smaller part by only keeping one in 10 streamlines
    """
    if verbose:
        #headers="From: %s\r\nTo: %s\r\nSubject:Reached the point where we start saving the file!\r\n\r\n" % (useremail,useremail)
        txt = 'About to save streamlines at ' + outpathsubject
        send_mail(txt,subject="Seed computation" )

    if saved_tracts == "small" or saved_tracts == "both":
        sg_small = lambda: (s for i, s in enumerate(streamlines_generator) if i % 10 == 0)
        outpathtrk = outpathsubject + "_bmCSA_detr_small_" + stringstep + ".trk"
        save_trk_heavy_duty(outpathtrk, streamlines=sg_small,
                            affine=affine, header=trkheader,
                            shape=mask.shape, vox_size=vox_size)
        txt="Tract files were saved at " + outpathtrk

    else:
        outpathtrk = None
    if saved_tracts == "large" or saved_tracts == "both" or saved_tracts == "all":
        sg = lambda: (s for s in streamlines_generator)
        outpathtrk = outpathsubject + "bmCSA_detr_all_" + stringstep + ".trk"
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=trkheader,
                            shape=mask.shape, vox_size=vox_size)
        if verbose:
            print("Tract files were saved at "+outpathtrk)
    if saved_tracts == "none" or saved_tracts is None:
        print("Tract files were not saved")

    # save everything - will generate a 20+ GBytes of data - hard to manipulate

    # possibly add parameter in csv file or other to decide whether to save large tractogram file
    # outpathfile=outpath+subject+"bmCSA_detr"+stringstep+".trk"
    # myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))
   """
    duration3 = time() - t2
    if verbose:
        print(duration3)
        print(subject + ' Tracking duration %.3f' % (duration3,))
        send_mail("Finished file save at "+outpathtrk+" with tracking duration of " + str(duration3) + "seconds",
                  subject="file save update" )


    return outpathtrk
