#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:14:39 2020

@author: Jacques Stout

Part of DTC pipeline
Creates and generally handles results of said pipeline onto figures
"""


import numpy as np
from dipy.viz import window, actor

# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.

import matplotlib.pyplot as plt


from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib
#import dipy.tracking.life as life
from dipy.viz import colormap as cmap
from BIAC_tools import send_mail
from dipy.reconst.shore import ShoreModel
from dipy.io.image import load_nifti, save_nifti
import os
from dipy.segment.clustering import QuickBundles
from dipy.segment.bundles import RecoBundles
from itertools import combinations, groupby
from dipy.data.fetcher import fetch_bundles_2_subjects, read_bundles_2_subjects
from dipy.tracking.streamline import Streamlines


def win_callback(obj, event):
    global size
    if size != obj.GetSize():
        size_old = size
        size = obj.GetSize()
        size_change = [size[0] - size_old[0], 0]
        panel.re_align(size_change)


def denoise_fig(data,denoised_arr,type='macenko',savefigpath='none'):
    """ Sets up figure that shows differences between two data arrays
    with assumption that one is standard and other is after denoising

    Parameters
    data: arr
    Initial data
    denoised_arr: arr
    denoised data array
    type: str
    Type of denoising
    savefigpath: str
    path of saved output of figure, if none no figure is saved
    """
    sli = data.shape[2] // 2
    gra = data.shape[3] - 1
    orig = data[:, :, sli, gra]
    den = denoised_arr[:, :, sli, gra]
    rms_diff = np.sqrt((orig - den) ** 2)

    if show_orig_denoised:
        fig1, ax = plt.subplots(1, 3, figsize=(12, 6),
                                subplot_kw={'xticks': [], 'yticks': []})

        fig1.subplots_adjust(hspace=0.3, wspace=0.05)

        ax.flat[0].imshow(orig.T, cmap='gray', interpolation='none',
                          origin='lower')
        ax.flat[0].set_title('Original')
        ax.flat[1].imshow(den.T, cmap='gray', interpolation='none',
                          origin='lower')
        ax.flat[1].set_title('Denoised Output')
        if type == 'macenko':
            ax.flat[2].imshow(rms_diff.T, cmap='gray', interpolation='none',
                          origin='lower')
            ax.flat[2].set_title('Residuals')
        if type == 'gibbs':
            ax.flat[2].imshow(data_corrected[:, :, 0, 0].T - data_slices[:, :, 0, 0].T,
                              cmap='gray', origin='lower', vmin=-500, vmax=500)
            ax.flat[2].set_title('Gibbs residuals')
        if savefigpath.lower=='none':
            pass
        else:
            fig1.savefig(savefigpath)

def viewstreamlines_anat(streamlines_full, anat_path, affine, ratio = 1, threshold = 10., verbose = False):

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)

    #colors = ['white', 'cadmium_red_deep', 'misty_rose', 'slate_grey_dark', 'ivory_black', 'chartreuse']
    colors = [window.colors.white, window.colors.cadmium_red_deep, window.colors.misty_rose, window.colors.slate_grey_dark, window.colors.ivory_black, window.colors.chartreuse]
    streamline_cut = []
    i = 0
    if ratio != 1:
        for streamline in streamlines_full:
            if i % ratio == 0:
                streamline_cut.append(streamline)
            i += 1
    else:
        streamline_cut = streamlines_full
    qb = QuickBundles(threshold=threshold)
    clusters = qb.cluster(streamline_cut)

    if verbose:
        print("Nb. clusters:", len(clusters))
        print("Cluster sizes:", map(len, clusters))
        print("Small clusters:", clusters < 10)
        print("Streamlines indices of the first cluster:\n", clusters[0].indices)
        print("Centroid of the last clustker:\n", clusters[-1].centroid)

    j=0
    scene = window.Scene()
    scene.add(actor.streamtube(streamline_cut, colors[j]))
    slicer_opacity = 0.6
    j += 1

    if isinstance(anat_path, str) and os.path.exists(anat_path):
        anat_nifti = load_nifti(anat_path)
        try:
            data = anat_nifti.data
        except AttributeError:
            data = anat_nifti[0]
        if affine is None:
            try:
                affine = anat_nifti.affine
            except AttributeError:
                affine = anat_nifti[1]
    else:
        data = anat_path


    shape = np.shape(data)
    if np.size(shape)==4:
        data=data[:,:,:,0]
    image_actor_z = actor.slicer(data, affine)
    image_actor_z.opacity(slicer_opacity)

    image_actor_x = image_actor_z.copy()
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display_extent(x_midpoint,
                                 x_midpoint, 0,
                                 shape[1] - 1,
                                 0,
                                 shape[2] - 1)

    image_actor_y = image_actor_z.copy()
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display_extent(0, shape[0] - 1,
                                 y_midpoint,
                                 y_midpoint,
                                 0,
                                 shape[2] - 1)

    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)
    global size
    size = scene.GetSize()
    show_m = window.ShowManager(scene, size=(1200, 900))
    show_m.initialize()

    interactive = True
    interactive = False
    if interactive:

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()
    else:
        window.record(scene, out_path='bundles_and_3_slices.png', size=(1200, 900),
                      reset_camera=False)


def connective_streamlines_figuremaker(allstreamlines, ROI_streamlines, ROI_names, anat_path, threshold=10., verbose=False):

    #streamlines = Streamlines(res['af.left'])
    #streamlines.extend(res['cst.right'])
    #streamlines.extend(res['cc_1'])
    world_coords = True


    # Cluster sizes: [64, 191, 47, 1]

    # Small clusters: array([False, False, False, True], dtype=bool)

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)

    colors = ['white', 'cadmium_red_deep', 'misty_rose', 'slate_grey_dark', 'ivory_black', 'chartreuse']
    colors = [window.colors.white, window.colors.cadmium_red_deep, window.colors.misty_rose, window.colors.slate_grey_dark, window.colors.ivory_black, window.colors.chartreuse]
    i = 0
    for ROI in ROI_streamlines:
        ROI_streamline = allstreamlines[ROI]
        qb = QuickBundles(threshold=threshold)
        clusters = qb.cluster(ROI_streamline)
        if verbose:
            print("Nb. clusters:", len(clusters))
            print("Cluster sizes:", map(len, clusters))
            print("Small clusters:", clusters < 10)
            print("Streamlines indices of the first cluster:\n", clusters[0].indices)
            print("Centroid of the last cluster:\n", clusters[-1].centroid)
        #if not world_coords:
        #    from dipy.tracking.streamline import transform_streamlines
        #    streamlines = transform_streamlines(ROI_streamline, np.linalg.inv(affine))
        scene = window.Scene()
        #stream_actor = actor.line(ROI_streamline)
        #scene.add(actor.streamtube(ROI_streamline, window.colors.misty_rose))
        scene.add(actor.streamtube(ROI_streamline, colors[i]))

        #if not world_coords:
        #    image_actor_z = actor.slicer(data, affine=np.eye(4))
        #else:
        #    image_actor_z = actor.slicer(data, affine)

        slicer_opacity = 0.6
        i = i + 1

    anat_nifti = load_nifti(anat_path)
    try:
        data = anat_nifti.data
    except AttributeError:
        data = anat_nifti[0]
    try:
        affine = anat_nifti.affine
    except AttributeError:
        affine = anat_nifti[1]
    shape = np.shape(data)
    image_actor_z = actor.slicer(data[:,:,:,0], affine)
    image_actor_z.opacity(slicer_opacity)

    image_actor_x = image_actor_z.copy()
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display_extent(x_midpoint,
                                 x_midpoint, 0,
                                 shape[1] - 1,
                                 0,
                                 shape[2] - 1)

    image_actor_y = image_actor_z.copy()
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display_extent(0, shape[0] - 1,
                                 y_midpoint,
                                 y_midpoint,
                                 0,
                                 shape[2] - 1)

    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)
    global size
    size = scene.GetSize()
    show_m = window.ShowManager(scene, size=(1200, 900))
    show_m.initialize()

    interactive = True
    interactive = False
    if interactive:

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()
    else:
        window.record(scene, out_path='bundles_and_3_slices.png', size=(1200, 900),
                      reset_camera=False)


def shore_scalarmaps(data, gtab, outpath_fig, verbose = None):

    #bvecs[1:] = (bvecs[1:] /
    #             np.sqrt(np.sum(bvecs[1:] * bvecs[1:], axis=1))[:, None])
    #gtab = gradient_table(bvals, bvecs)
    if verbose:
        print('data.shape (%d, %d, %d, %d)' % data.shape)
    asm = ShoreModel(gtab)
    #Let’s just use only one slice only from the data.

    dataslice = data[30:70, 20:80, data.shape[2] // 2]
    #Fit the signal with the model and calculate the SHORE coefficients.

    asmfit = asm.fit(dataslice)
    #Calculate the analytical RTOP on the signal that corresponds to the integral of the signal.
    if verbose:
        print('Calculating... rtop_signal')
    rtop_signal = asmfit.rtop_signal()
    #Now we calculate the analytical RTOP on the propagator, that corresponds to its central value.

    if verbose:
        print('Calculating... rtop_pdf')
    rtop_pdf = asmfit.rtop_pdf()
    #In theory, these two measures must be equal, to show that we calculate the mean square error on this two measures.

    mse = np.sum((rtop_signal - rtop_pdf) ** 2) / rtop_signal.size
    if verbose:
        print("MSE = %f" % mse)
    MSE = 0.000000

    #Let’s calculate the analytical mean square displacement on the propagator.
    if verbose:
        print('Calculating... msd')
    msd = asmfit.msd()
    #Show the maps and save them to a file.

    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_subplot(2, 2, 1, title='rtop_signal')
    ax1.set_axis_off()
    ind = ax1.imshow(rtop_signal.T, interpolation='nearest', origin='lower')
    plt.colorbar(ind)
    ax2 = fig.add_subplot(2, 2, 2, title='rtop_pdf')
    ax2.set_axis_off()
    ind = ax2.imshow(rtop_pdf.T, interpolation='nearest', origin='lower')
    plt.colorbar(ind)
    ax3 = fig.add_subplot(2, 2, 3, title='msd')
    ax3.set_axis_off()
    ind = ax3.imshow(msd.T, interpolation='nearest', origin='lower', vmin=0)
    plt.colorbar(ind)
    print("about to save")
    plt.savefig(outpath_fig)
    print("save done")


def show_bundles(bundles, colors=None, show=True, fname=None,fa = False, str_tube = False, size=(900,900)):
    """ Displays bundles
    
    Parameters
    ---------
    bundles: bundles object
    """
    ren = window.Renderer()
    ren.SetBackground(1., 1, 1)
    if str_tube:
        bundle_actor = actor.streamtube(bundles, colors, linewidth=0.5)
        ren.add(bundle_actor)
    else:
        for (i, bundle) in enumerate(bundles):
            color = colors[i]
    #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05

            lines_actor = actor.line(bundle, color,linewidth=2.5)
            #lines_actor.RotateX(-90)
            #lines_actor.RotateZ(90)
            ren.add(lines_actor)
        
    if fa:
        fa, affine_fa= load_nifti('/Users/alex/code/Wenlin/data/wenlin_results/bmfaN54900.nii.gz')
        fa_actor = actor.slicer(fa, affine_fa)
        ren.add(fa_actor)
    
    if show:
        window.show(ren)
    if fname is not None:
        sleep(1)
        window.record(ren, n_frames=1, out_path=fname, size=size)

def window_show_test(bundles, mask_roi, anat, interactive= True, outpath=None):
    """

    :param bundles:
    :param mask_roi:
    :param anat:
    :param interactive:
    :param outpath:
    :return:
    """

    candidate_streamlines_actor = actor.streamtube(bundles,
                                                   cmap.line_colors(candidate_sl))
    ROI_actor = actor.contour_from_roi(mask_roi, color=(1., 1., 0.),
                                          opacity=0.5)

    ren = window.Renderer()
    if anat:
        vol_actor = actor.slicer(anat)
        vol_actor.display(x=40)
        vol_actor2 = vol_actor.copy()
        vol_actor2.display(z=35)

    # Add display objects to canvas
    ren.add(candidate_streamlines_actor)
    ren.add(ROI_actor)
    ren.add(vol_actor)
    ren.add(vol_actor2)
    if outpath is not None:
        window.record(ren, n_frames=1,
                  out_path=outpath,
                  size=(800, 800))
    if interactive:
        window.show(ren)

def LifEcreate_fig(fiber_fit_beta,mean_rmse,model_rmse, vox_coords, dwidata, subject, t1_data=None, outpathfig=None, interactive=False, strproperty="_", verbose=False):

    #fiber_fit_beta_path = glob.glob(pickles_folder + '/*beta.p')[0]
    #mean_rmse_path = glob.glob(pickles_folder + '/*mean_rmse.p')[0]
    #model_rmse_path = glob.glob(pickles_folder + '/*model_rmse.p')[0]
    #fiber_fit_beta = pickle.load(open(fiber_fit_beta_path, "rb"))
    #mean_rmse = pickle.load(open(mean_rmse_path, "rb"))
    #model_rmse = pickle.load(open(model_rmse_path, "rb"))

    fig, ax = plt.subplots(1)
    ax.hist(fiber_fit_beta, bins=100, histtype='step')
    ax.set_xlabel('Fiber weights')
    ax.set_ylabel('# fibers')
    #ROI_actor = actor.contour_from_roi(roimask, color=(1., 1., 0.),
    #                                      opacity=0.5)
    #sizebeta=getsize(fiber_fit_beta)
    if interactive:
        plt.show()
    if outpathfig is not None:
        histofig_path = (outpathfig + subject + strproperty + "_beta_histogram.png")
        fig.savefig(histofig_path)
        if verbose:
            txt="file saved at "+histofig_path
            print(txt)
            send_mail(txt,subject="LifE save msg ")

    """
    vol_actor = actor.slicer(t1_data)
    
    vol_actor.display(x=40)
    vol_actor2 = vol_actor.copy()
    vol_actor2.display(z=35)        
    
    """

    fig, ax = plt.subplots(1)
    ax.hist(mean_rmse - model_rmse, bins=100, histtype='step')
    ax.text(0.2, 0.9, 'Median RMSE, mean model: %.2f' % np.median(mean_rmse),
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax.transAxes)
    ax.text(0.2, 0.8, 'Median RMSE, LiFE: %.2f' % np.median(model_rmse),
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax.transAxes)
    ax.set_xlabel('RMS Error')
    ax.set_ylabel('# voxels')
    if interactive:
        plt.show()
    if outpathfig is not None:
        errorhistofig_path=(outpathfig + subject + strproperty + "_error_histograms.png")
        fig.savefig(errorhistofig_path)
        if verbose:
            txt="file saved at "+errorhistofig_path
            print(txt)
            send_mail(txt,subject="LifE save msg ")

    runspatialerrors=True
    try:
        dwidata.shape[:3]
    except AttributeError:
        runspatialerrors=False

    if runspatialerrors:
        vol_model = np.ones(dwidata.shape[:3]) * np.nan
        vol_model[vox_coords[:, 0],
              vox_coords[:, 1],
              vox_coords[:, 2]] = model_rmse
        vol_mean = np.ones(dwidata.shape[:3]) * np.nan
        vol_mean[vox_coords[:, 0],
            vox_coords[:, 1],
            vox_coords[:, 2]] = mean_rmse
        vol_improve = np.ones(dwidata.shape[:3]) * np.nan
        vol_improve[vox_coords[:, 0],
            vox_coords[:, 1],
            vox_coords[:, 2]] = mean_rmse - model_rmse
        sl_idx = 49
        fig = plt.figure()
        fig.subplots_adjust(left=0.05, right=0.95)
        ax = AxesGrid(fig, 111,
            nrows_ncols=(1, 3),
            label_mode="1",
            share_all=True,
            cbar_location="top",
            cbar_mode="each",
            cbar_size="10%",
            cbar_pad="5%")

        ax[0].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
        im = ax[0].matshow(np.rot90(vol_model[sl_idx, :, :]), cmap=matplotlib.cm.hot)
        ax.cbar_axes[0].colorbar(im)
        ax[1].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
        im = ax[1].matshow(np.rot90(vol_mean[sl_idx, :, :]), cmap=matplotlib.cm.hot)
        ax.cbar_axes[1].colorbar(im)
        ax[2].matshow(np.rot90(t1_data[sl_idx, :, :]), cmap=matplotlib.cm.bone)
        im = ax[2].matshow(np.rot90(vol_improve[sl_idx, :, :]),
            cmap=matplotlib.cm.RdBu)
        ax.cbar_axes[2].colorbar(im)
        for lax in ax:
            lax.set_xticks([])
            lax.set_yticks([])
        if outpathfig is not None:
            histofig_path=(outpathfig+ subject+ strproperty + "_spatial_errors.png")
            fig.savefig(histofig_path)
        if verbose:
            txt="spatial errors figure saved at " + histofig_path
            print(txt)
            send_mail(txt,subject="LifE save msg ")        
