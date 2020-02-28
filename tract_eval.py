#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:17:04 2020

@author: Jacques Stout

"""


import numpy as np
import pickle

from dipy.io.utils import (create_tractogram_header)
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
from dipy.tracking import utils
import matplotlib.pyplot as plt

#import dipy.tracking.life as life
import JSdipy.tracking.life as life
import dipy.core.optimize as opt

from BIAC_tools import send_mail, getsize
from tract_save import save_trk_heavy_duty
from figures_handler import LifEcreate_fig

def bundle_coherence(data,t1_data,hardi_img, gtab, labels_img,affine):

    # Enables/disables interactive visualization
    interactive = False
    # Fix seed
    np.random.seed(1)

    # Read data
    #hardi_img, gtab, labels_img = read_stanford_labels()
    #affine = hardi_img.affine
    #t1 = read_stanford_t1()
    #t1_data = t1.get_data()

    # Select a relevant part of the data (left hemisphere)
    # Coordinates given in x bounds, y bounds, z bounds
    dshape = data.shape[:-1]
    xa, xb, ya, yb, za, zb = [15, 42, 10, 65, 18, 65]
    data_small = data[xa:xb, ya:yb, za:zb]
    selectionmask = np.zeros(dshape, 'bool')
    selectionmask[xa:xb, ya:yb, za:zb] = True

    # Perform CSA
    from dipy.reconst.shm import CsaOdfModel
    from dipy.data import default_sphere
    from dipy.direction import peaks_from_model

    csa_model = CsaOdfModel(gtab, sh_order=6)
    csa_peaks = peaks_from_model(csa_model, data, default_sphere,
                                 relative_peak_threshold=.6,
                                 min_separation_angle=45,
                                 mask=selectionmask)

    # Stopping Criterion
    from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion

    stopping_criterion = ThresholdStoppingCriterion(csa_peaks.gfa, 0.25)

    # Perform CSD on the original data
    from dipy.reconst.csdeconv import auto_response
    from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel

    response, ratio = auto_response(gtab, data, roi_radius=10, fa_thr=0.7)
    csd_model = ConstrainedSphericalDeconvModel(gtab, response)
    csd_fit = csd_model.fit(data_small)
    csd_fit_shm = np.lib.pad(csd_fit.shm_coeff, ((xa, dshape[0]-xb),
                                                 (ya, dshape[1]-yb),
                                                 (za, dshape[2]-zb),
                                                 (0, 0)), 'constant')

    # Probabilistic direction getting for fiber tracking
    from dipy.direction import ProbabilisticDirectionGetter

    prob_dg = ProbabilisticDirectionGetter.from_shcoeff(csd_fit_shm,
                                                        max_angle=30.,
                                                        sphere=default_sphere)


    # Set a seed region region for tractography.

    mask = np.zeros(data.shape[:-1], 'bool')
    rad = 3
    mask[26-rad:26+rad, 29-rad:29+rad, 31-rad:31+rad] = True
    seeds = utils.seeds_from_mask(mask, affine, density=[4, 4, 4])

    # Perform tracking using Local Tracking
    from dipy.tracking.local_tracking import LocalTracking

    streamlines_generator = LocalTracking(prob_dg, stopping_criterion, seeds,
                                          affine, step_size=.5)

    # Compute streamlines.
    from dipy.tracking.streamline import Streamlines
    streamlines = Streamlines(streamlines_generator)

    # Set a mask for the lateral geniculate nucleus (LGN)
    mask_lgn = np.zeros(data.shape[:-1], 'bool')
    rad = 5
    mask_lgn[35-rad:35+rad, 42-rad:42+rad, 28-rad:28+rad] = True

    # Select all the fibers that enter the LGN and discard all others
    filtered_fibers2 = utils.near_roi(streamlines, affine, mask_lgn, tol=1.8)

    sfil = []
    for i in range(len(streamlines)):
        if filtered_fibers2[i]:
            sfil.append(streamlines[i])
    streamlines = Streamlines(sfil)

    # Compute lookup table
    from dipy.denoise.enhancement_kernel import EnhancementKernel

    D33 = 1.0
    D44 = 0.02
    t = 1
    k = EnhancementKernel(D33, D44, t)

    # Apply FBC measures
    from dipy.tracking.fbcmeasures import FBCMeasures

    fbc = FBCMeasures(streamlines, k)

    # Calculate LFBC for original fibers
    fbc_sl_orig, clrs_orig, rfbc_orig = \
      fbc.get_points_rfbc_thresholded(0, emphasis=0.01)

    # Apply a threshold on the RFBC to remove spurious fibers
    fbc_sl_thres, clrs_thres, rfbc_thres = \
      fbc.get_points_rfbc_thresholded(0.125, emphasis=0.01)

    # Visualize the results
    from dipy.viz import window, actor

    # Create renderer
    ren = window.Renderer()

    # Original lines colored by LFBC
    lineactor = actor.line(fbc_sl_orig, clrs_orig, linewidth=0.2)
    ren.add(lineactor)

    # Horizontal (axial) slice of T1 data
    vol_actor1 = actor.slicer(t1_data, affine=affine)
    vol_actor1.display(z=20)
    ren.add(vol_actor1)

    # Vertical (sagittal) slice of T1 data
    vol_actor2 = actor.slicer(t1_data, affine=affine)
    vol_actor2.display(x=35)
    ren.add(vol_actor2)

    # Show original fibers
    ren.set_camera(position=(-264, 285, 155),
                   focal_point=(0, -14, 9),
                   view_up=(0, 0, 1))
    window.record(ren, n_frames=1, out_path='OR_before.png', size=(900, 900))
    if interactive:
        window.show(ren)

    # Show thresholded fibers
    ren.rm(lineactor)
    ren.add(actor.line(fbc_sl_thres, clrs_thres, linewidth=0.2))
    window.record(ren, n_frames=1, out_path='OR_after.png', size=(900, 900))
    if interactive:
        window.show(ren)

def LiFEvaluation(dwidata, trk_streamlines, gtab, subject="lifesubj", header=None, roimask=None, affine=None,
                  display = True, outpathpickle = None, outpathtrk = None, processes = 1,
                  outpathfig=None, verbose = None):

    """     Implementation of Linear Fascicle Evaluation, outputs histograms, evals

    Parameters
    ----------
    dwidata : array
        output trk filename
    trkdata : array
    gtab : array og bval & bvec table
    outpath: string
    folder location for resulting values and figures
    display : boolean, optional
    Condition to display the results (default = False)
    savefig: boolean, optional
    Condition to save the results in outpath (default = True)

    Defined by Pestilli, F., Yeatman, J, Rokem, A. Kay, K. and Wandell B.A. (2014).
    Validation and statistical inference in living connectomes and recreated by Dipy

    :param dwidata: array of diffusion data
    :param trkdata: array of tractography data obtained from dwi
    :param gtab: bval & bvec table
    :param outpath: location to save analysis outputs
    :param display:
    :param savefig:
    :return:
    """

    """""""""
    if not op.exists('lr-superiorfrontal.trk'):
        
    else:
        # We'll need to know where the corpus callosum is from these variables:
        from dipy.data import (read_stanford_labels,
                               fetch_stanford_t1,
                               read_stanford_t1)
        hardi_img, gtab, labels_img = read_stanford_labels()
        labels = labels_img.get_data()
        cc_slice = labels == 2
        fetch_stanford_t1()
        t1 = read_stanford_t1()
        t1_data = t1.get_data()
        data = hardi_img.get_data()
    """""""""
    # Read the candidates from file in voxel space:

    if roimask is None:
        roimask = dwidata > 0
    else:
        dwidataroi = dwidata * np.repeat(roimask[:, :, :, None], np.shape(dwidata)[3], axis=3)

    print("verbose: "+str(verbose) +" outpathpickle: "+ str(outpathpickle))
    fiber_model = life.FiberModel(gtab)
    # inv_affine must be used if the streamlines are in the world space, and thus we must useapply the inverse affine of dwi
    #when comparing the diffusion directions gtab and the voxels of trk
    #inv_affine = np.linalg.inv(hardi_img.affine)

    #fiber_fit will fit the streamlines to the original diffusion data and
    if verbose:
        txt="Begin the evaluation over "+str(np.size(trk_streamlines))+" streamlines"
        print(txt)
        send_mail(txt,subject="LifE start msg ")

    fiber_fit = fiber_model.fit(dwidata, trk_streamlines, affine=np.eye(4), processes=processes, verbose=verbose)
    #fiber_fit_roi = fiber_model.fit(dwidataroi, trk_streamlines, affine=np.eye(4), processes=processes, verbose=verbose)
    optimized_sl = list(np.array(trk_streamlines)[np.where(fiber_fit.beta > 0)[0]])
    plt.ioff()
    if verbose:
        txt="End of the evaluation over "+str(np.size(trk_streamlines))
        print(txt)
        send_mail(txt,subject="LifE status msg ")
    if outpathtrk is not None:
        outpathfile = str(outpathtrk) + subject + strproperty + "_lifeopt_test.trk"
        myheader = create_tractogram_header(outpathfile, *header)
        optimized_sl_gen = lambda: (s for s in optimized_sl)
        save_trk_heavy_duty(outpathfile, streamlines=optimized_sl_gen,
                                affine=affine, header=myheader)
        txt = ("Saved final trk at "+ outpathfile)
        print(txt)
        send_mail(txt,subject="LifE save msg ")
        """
        except TypeError:
            txt=('Could not save new tractogram file, header of original trk file not properly implemented into '
                  'LifEvaluation')
            print(txt)
            send_mail(txt,subject="LifE error msg ")
        """
    """
    if interactive:
        ren = window.Renderer()
        ren.add(actor.streamtube(optimized_sl, cmap.line_colors(optimized_sl)))
        ren.add(ROI_actor)
        #ren.add(vol_actor)
        if interactive:
            window.show(ren)      
        if outpathfig is not None:
            print("reached windowrecord")
            window.record(ren, n_frames=1, out_path=outpathfig +'_life_optimized.png',
                size=(800, 800))
            print("did window record")
    """
    maxsize_var = 20525023825

    sizebeta = getsize(fiber_fit.beta)
    if sizebeta<maxsize_var:
        picklepath = outpathpickle + subject + '_beta.p'
        txt=("fiber_fit.beta saved at "+picklepath)
        pickle.dump(fiber_fit.beta, open(picklepath, "wb"))
        if verbose:
            print(txt)
            send_mail(txt,subject="LifE save msg ")
    else:
        txt=("Object fiber_fit.beta exceeded the imposed the 20GB limit with a size of: "+str(sizebeta/(10^9))+ "GB")
        print(txt)
        send_mail(txt,subject="LifE error msg")

    sizecoords=getsize(fiber_fit.vox_coords)
    if sizecoords<maxsize_var:
        picklepath = outpathpickle + subject + '_voxcoords.p'
        txt=("fiber_fit.voxcoords saved at "+picklepath)
        pickle.dump(fiber_fit.vox_coords, open(picklepath, "wb"))
        if verbose:
            print(txt)
            send_mail(txt,subject="LifE save msg ")
    else:
        txt=("Object fiber_fit.beta exceeded the imposed the 20GB limit with a size of: "+str(sizebeta/(10^9))+ "GB")
        print(txt)
        send_mail(txt,subject="LifE error msg")

    #predict diffusion data based on new model
    model_predict = fiber_fit.predict() #possible to predict based on different gtab or base signal (change gtab, S0)
    model_error = model_predict - fiber_fit.data    #compare original dwi data and the model fit, calculate error
    model_rmse = np.sqrt(np.mean(model_error[:, 10:] ** 2, -1))     #this is good, but must check ways to avoid overfitting
                                                            #how does the model get built? add lasso? JS

    beta_baseline = np.zeros(fiber_fit.beta.shape[0]) #baseline assumption where the streamlines weight is 0
    pred_weighted = np.reshape(opt.spdot(fiber_fit.life_matrix, beta_baseline),
                               (fiber_fit.vox_coords.shape[0],
                                np.sum(~gtab.b0s_mask)))
    mean_pred = np.empty((fiber_fit.vox_coords.shape[0], gtab.bvals.shape[0]))
    S0 = fiber_fit.b0_signal

    mean_pred[..., gtab.b0s_mask] = S0[:, None]
    mean_pred[..., ~gtab.b0s_mask] = \
        (pred_weighted + fiber_fit.mean_signal[:, None]) * S0[:, None]
    mean_error = mean_pred - fiber_fit.data
    mean_rmse = np.sqrt(np.mean(mean_error ** 2, -1))

    size_meanrmse=getsize(mean_rmse)
    if size_meanrmse<maxsize_var:
        picklepath = outpathpickle + subject + '_mean_rmse.p'
        txt=("mean_rmse saved at "+picklepath)
        pickle.dump(mean_rmse, open(picklepath, "wb"))
        if verbose:
            print(txt)
            send_mail(txt,subject="LifE save msg ")
    else:
        txt=("Object mean_rmse exceeded the imposed the 20GB limit with a size of: "+str(size_meanrmse/(10^9)) +" GB")
        print(txt)
        send_mail(txt,subject="LifE error msg")   

    size_modelrmse=getsize(model_rmse)
    if size_modelrmse<maxsize_var:
        picklepath = outpathpickle + subject + '_model_rmse.p'
        txt=("model_rmse saved at "+picklepath)
        pickle.dump(model_rmse, open(picklepath, "wb"))
        if verbose:
            print(txt)
            send_mail(txt,subject="LifE save msg ")
    else:
        txt=("Object model_rmse exceeded the imposed the 20GB limit with a size of: "+str(size_modelrmse/(10^9)) +" GB")
        print(txt)
        send_mail(txt,subject="LifE error msg")   

    if outpathfig is not None:
        LifEcreate_fig(fiber_fit.beta, mean_rmse, model_rmse, fiber_fit.vox_coords, dwidata, subject, t1_data = dwidata[:,:,:,0], outpathfig=outpathfig, interactive=False, verbose=verbose)
    return model_error, mean_error