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

#from dipy.denoise.localpca import mppca
#import dipy.tracking.life as life

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

def QCSA_tractmake(data,affine,vox_size,gtab,mask,trkheader,step_size,peak_processes,outpathsubject,saved_tracts="small",verbose=None,subject = 'NA'):
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
    csa_peaks = peaks_from_model(model=csa_model,
                                 data=data,
                                 sphere=peaks.default_sphere,  # issue with complete sphere
                                 mask=mask,
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
    classifier = BinaryStoppingCriterion(mask)

    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    # seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))
    # why are those not binary?
    if verbose:
        print('Computing seeds')
    seeds = utils.seeds_from_mask(mask, density=1,
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

    # the function above will bring all streamlines in memory
    # streamlines = Streamlines(streamlines_generator)

    # save a smaller part by only keeping one in 10 streamlines
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
        print("Tract files were saved at " + outpathtrk)
    else:
        outpathtrk = None
    if saved_tracts == "large" or saved_tracts == "both" or saved_tracts == "all":
        sg = lambda: (s for s in streamlines_generator)
        outpathtrk = outpathsubject + "bmCSA_detr_all" + stringstep + ".trk"
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

    duration3 = time() - t3
    if verbose:
        print(duration3)
        print(subject + ' Tracking duration %.3f' % (duration3,))
        send_mail("Finished file save at "+outpathtrk+" with tracking duration of " + str(duration3) + "seconds",
                  subject="file save update" )


    return outpathtrk