#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:14:39 2020

@author: Jacques Stout
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
        histofig_path = (outpathfig + subject + strproperty + "beta_histogram.png")
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
        errorhistofig_path=(outpathfig + subject + strproperty + "error_histograms.png")
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
            histofig_path=(outpathfig+ subject+ strproperty + "spatial_errors.png")
            fig.savefig(histofig_path)
        if verbose:
            txt="spatial errors figure saved at " + histofig_path
            print(txt)
            send_mail(txt,subject="LifE save msg ")        
