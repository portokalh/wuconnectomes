#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:17:04 2020

@author: Jacques Stout
Functions for evaluating the results of other processes
Is location for Life (inhouse parallel processing version), bundle coherence, others
A bit left behind, has some WIP functions. May need fixing in places.
"""


import numpy as np
import pickle

from dipy.io.utils import create_tractogram_header
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
from dipy.tracking import utils
import matplotlib.pyplot as plt
from numpy import ravel_multi_index

#import dipy.tracking.life as life
#import JSdipy.tracking.life as life
import dipy.core.optimize as opt

from BIAC_tools import send_mail, getsize
from tract_save import save_trk_heavy_duty
from figures_handler import LifEcreate_fig
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
from collections import defaultdict, OrderedDict

from dipy.denoise.enhancement_kernel import EnhancementKernel
from dipy.tracking.fbcmeasures import FBCMeasures
from dipy.viz import window, actor
from dipy.segment.clustering import QuickBundles
from dipy.segment.bundles import RecoBundles
from itertools import combinations, groupby


def ndbincount(x, weights=None, shape=None):
    """Like bincount, but for nd-indices.

    Parameters
    ----------
    x : array_like (N, M)
        M indices to a an Nd-array
    weights : array_like (M,), optional
        Weights associated with indices
    shape : optional
        the shape of the output
    """
    x = np.asarray(x)
    if shape is None:
        shape = x.max(1) + 1

    x = ravel_multi_index(x, shape)
    out = np.bincount(x, weights, minlength=np.prod(shape))
    out.shape = shape

    return out

def connectivity_selection_getsl(streamlines, affine, label_volume,
                            symmetric=True, return_mapping=True,
                        mapping_as_streamlines=False):
    """Counts the streamlines that start and end at each label pair.

    Parameters
    ----------
    streamlines : sequence
        A sequence of streamlines.
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline coordinates.
        The voxel_to_rasmm matrix, typically from a NIFTI file.
    label_volume : ndarray
        An image volume with an integer data type, where the intensities in the
        volume map to anatomical structures.
    labels : tuple (2,1)
        The labels that are to be isolated
    inclusive: bool
        Whether to analyze the entire streamline, as opposed to just the
        endpoints. Allowing this will increase calculation time and mapping
        size, especially if mapping_as_streamlines is True. False by default.
    symmetric : bool, True by default
        Symmetric means we don't distinguish between start and end points. If
        symmetric is True, ``matrix[i, j] == matrix[j, i]``.
    return_mapping : bool, False by default
        If True, a mapping is returned which maps matrix indices to
        streamlines.
    mapping_as_streamlines : bool, False by default
        If True voxel indices map to lists of streamline objects. Otherwise
        voxel indices map to lists of integers.

    Returns
    -------
    matrix : ndarray
        The number of connection between each pair of regions in
        `label_volume`.
    mapping : defaultdict(list)
        ``mapping[i, j]`` returns all the streamlines that connect region `i`
        to region `j`. If `symmetric` is True mapping will only have one key
        for each start end pair such that if ``i < j`` mapping will have key
        ``(i, j)`` but not key ``(j, i)``.

    """
    # Error checking on label_volume
    kind = label_volume.dtype.kind
    labels_positive = ((kind == 'u') or
                       ((kind == 'i') and (label_volume.min() >= 0)))
    valid_label_volume = (labels_positive and label_volume.ndim == 3)
    if not valid_label_volume:
        raise ValueError("label_volume must be a 3d integer array with"
                         "non-negative label values")

    # If streamlines is an iterator
    if return_mapping and mapping_as_streamlines:
        streamlines = list(streamlines)

    label_dict = {}
    #singlecase = np.size(np.shape(label_vals)) == 1

    matrix_sl = np.empty((3,), dtype=object)
    for i, v in enumerate(matrix_sl):
        matrix_sl[i] = [v, i]
    for v in matrix_sl:
        v.append(34)

    if inclusive:
        # Create ndarray to store streamline connections
        edges = np.ndarray(shape=(3, 0), dtype=int)
        lin_T, offset = _mapping_to_voxel(affine)
        for sl, _ in enumerate(streamlines):
            # Convert streamline to voxel coordinates
            entire = _to_voxel_coordinates(streamlines[sl], lin_T, offset)
            i, j, k = entire.T
            if symmetric:
                # Create list of all labels streamline passes through
                entirelabels = list(OrderedDict.fromkeys(label_volume[i, j, k]))
                # Append all connection combinations with streamline number
                for comb in combinations(entirelabels, 2):
                    if singlecase:
                        if (comb == label_vals).all():
                            label_dict[tuple(label_vals)].append(sl)
                    else:
                        for label in label_vals:
                            if (comb == label).all():
                                label_dict[tuple(label)].append(sl)
                    edges = np.append(edges, [[comb[0]], [comb[1]], [sl]],
                                      axis=1)
            else:
                # Create list of all labels streamline passes through, keeping
                # order and whether a label was entered multiple times
                entirelabels = list(groupby(label_volume[i, j, k]))
                # Append connection combinations along with streamline number,
                # removing duplicates and connections from a label to itself
                combs = set(combinations([z[0] for z in entirelabels], 2))
                for comb in combs:
                    if comb[0] == comb[1]:
                        pass
                    else:
                        edges = np.append(edges, [[comb[0]], [comb[1]], [sl]],
                                          axis=1)
        if symmetric:
            edges[0:2].sort(0)
        mx = label_volume.max() + 1
        matrix = ndbincount(edges[0:2], shape=(mx, mx))

        if symmetric:
            matrix = np.maximum(matrix, matrix.T)
        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b, c) in enumerate(edges.T):
                mapping[a, b].append(c)
            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            return matrix, mapping

        return matrix
    else:
        # take the first and last point of each streamline
        endpoints = [sl[0::len(sl)-1] for sl in streamlines]

        # Map the streamlines coordinates to voxel coordinates
        lin_T, offset = _mapping_to_voxel(affine)
        endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)

        # get labels for label_volume
        i, j, k = endpoints.T
        endlabels = label_volume[i, j, k]
        if symmetric:
            endlabels.sort(0)
        mx = label_volume.max() + 1
        matrix = ndbincount(endlabels, shape=(mx, mx))
        if symmetric:
            matrix = np.maximum(matrix, matrix.T)

        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b) in enumerate(endlabels.T):
                mapping[a, b].append(i)

            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            # Return the mapping matrix and the mapping
            return matrix, mapping

        return matrix




def connectivity_selection(streamlines, affine, label_volume, label_vals,
                            symmetric=True, return_mapping=True,
                        mapping_as_streamlines=False):
    """Counts the streamlines that start and end at each label pair.

    Parameters
    ----------
    streamlines : sequence
        A sequence of streamlines.
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline coordinates.
        The voxel_to_rasmm matrix, typically from a NIFTI file.
    label_volume : ndarray
        An image volume with an integer data type, where the intensities in the
        volume map to anatomical structures.
    labels : tuple (2,1)
        The labels that are to be isolated
    inclusive: bool
        Whether to analyze the entire streamline, as opposed to just the
        endpoints. Allowing this will increase calculation time and mapping
        size, especially if mapping_as_streamlines is True. False by default.
    symmetric : bool, True by default
        Symmetric means we don't distinguish between start and end points. If
        symmetric is True, ``matrix[i, j] == matrix[j, i]``.
    return_mapping : bool, False by default
        If True, a mapping is returned which maps matrix indices to
        streamlines.
    mapping_as_streamlines : bool, False by default
        If True voxel indices map to lists of streamline objects. Otherwise
        voxel indices map to lists of integers.

    Returns
    -------
    matrix : ndarray
        The number of connection between each pair of regions in
        `label_volume`.
    mapping : defaultdict(list)
        ``mapping[i, j]`` returns all the streamlines that connect region `i`
        to region `j`. If `symmetric` is True mapping will only have one key
        for each start end pair such that if ``i < j`` mapping will have key
        ``(i, j)`` but not key ``(j, i)``.

    """
    # Error checking on label_volume
    kind = label_volume.dtype.kind
    labels_positive = ((kind == 'u') or
                       ((kind == 'i') and (label_volume.min() >= 0)))
    valid_label_volume = (labels_positive and label_volume.ndim == 3)
    if not valid_label_volume:
        raise ValueError("label_volume must be a 3d integer array with"
                         "non-negative label values")

    # If streamlines is an iterator
    if return_mapping and mapping_as_streamlines:
        streamlines = list(streamlines)

    label_dict = {}
    singlecase = np.size(np.shape(label_vals)) == 1

    if singlecase:
        label_dict[tuple(label_vals)] = []
    else:
        for label in label_vals:
            label_dict[label] = []

    if inclusive:
        # Create ndarray to store streamline connections
        edges = np.ndarray(shape=(3, 0), dtype=int)
        lin_T, offset = _mapping_to_voxel(affine)
        for sl, _ in enumerate(streamlines):
            # Convert streamline to voxel coordinates
            entire = _to_voxel_coordinates(streamlines[sl], lin_T, offset)
            i, j, k = entire.T
            if symmetric:
                # Create list of all labels streamline passes through
                entirelabels = list(OrderedDict.fromkeys(label_volume[i, j, k]))
                # Append all connection combinations with streamline number
                for comb in combinations(entirelabels, 2):
                    if singlecase:
                        if (comb == label_vals).all():
                            label_dict[tuple(label_vals)].append(sl)
                    else:
                        for label in label_vals:
                            if (comb == label).all():
                                label_dict[tuple(label)].append(sl)
                    edges = np.append(edges, [[comb[0]], [comb[1]], [sl]],
                                      axis=1)
            else:
                # Create list of all labels streamline passes through, keeping
                # order and whether a label was entered multiple times
                entirelabels = list(groupby(label_volume[i, j, k]))
                # Append connection combinations along with streamline number,
                # removing duplicates and connections from a label to itself
                combs = set(combinations([z[0] for z in entirelabels], 2))
                for comb in combs:
                    if comb[0] == comb[1]:
                        pass
                    else:
                        edges = np.append(edges, [[comb[0]], [comb[1]], [sl]],
                                          axis=1)
        if symmetric:
            edges[0:2].sort(0)
        mx = label_volume.max() + 1
        matrix = ndbincount(edges[0:2], shape=(mx, mx))

        if symmetric:
            matrix = np.maximum(matrix, matrix.T)
        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b, c) in enumerate(edges.T):
                mapping[a, b].append(c)
            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            return matrix, mapping

        return matrix
    else:
        # take the first and last point of each streamline
        endpoints = [sl[0::len(sl)-1] for sl in streamlines]

        # Map the streamlines coordinates to voxel coordinates
        lin_T, offset = _mapping_to_voxel(affine)
        endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)

        # get labels for label_volume
        i, j, k = endpoints.T
        endlabels = label_volume[i, j, k]
        if symmetric:
            endlabels.sort(0)
        mx = label_volume.max() + 1
        matrix = ndbincount(endlabels, shape=(mx, mx))
        if symmetric:
            matrix = np.maximum(matrix, matrix.T)

        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b) in enumerate(endlabels.T):
                mapping[a, b].append(i)

            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            # Return the mapping matrix and the mapping
            return matrix, mapping

        return matrix


def bundle_coherence(streamlines, affine, k, t1_data=None,interactive=False):

    # Compute lookup table

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
    if t1_data is not None:
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

    return k

def LiFEvaluation(dwidata, trk_streamlines, gtab, subject="lifesubj", header=None, roimask=None, affine=None,
                  display = True, outpathpickle = None, outpathtrk = None, processes = 1,
                  outpathfig=None, strproperty="", verbose = None):

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
        picklepath = outpathpickle + subject + strproperty + '_beta.p'
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
        picklepath = outpathpickle + subject + strproperty + '_voxcoords.p'
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
        picklepath = outpathpickle + subject + strproperty + '_mean_rmse.p'
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
        picklepath = outpathpickle + subject + strproperty + '_model_rmse.p'
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
        try:
            import matplotlib.pyplot as myplot
            fig, ax = plt.subplots(1)
            ax.hist(fiber_fit.beta, bins=100, histtype='step')
            LifEcreate_fig(fiber_fit.beta, mean_rmse, model_rmse, fiber_fit.vox_coords, dwidata, subject, t1_data = dwidata[:,:,:,0], outpathfig=outpathfig, interactive=False, strproperty=strproperty, verbose=verbose)
        except:
            print("Coult not launch life create fig, possibly qsub location (this is a template warning, to be improved upon")
    return model_error, mean_error



def launch_quickbundles(streamlines, outpath, ROIname="all", threshold = 10., labelmask = None, affine = np.eye(4), interactive = False):

    #qb = QuickBundles(threshold=10.)
    qb = QuickBundles(threshold=threshold)
    clusters = qb.cluster(streamlines)

    print("Nb. clusters:", len(clusters))
    print("Cluster sizes:", map(len, clusters))
    print("Small clusters:", clusters < 10)
    print("Streamlines indices of the first cluster:\n", clusters[0].indices)
    print("Centroid of the last cluster:\n", clusters[-1].centroid)

    # Cluster sizes: [64, 191, 47, 1]

    # Small clusters: array([False, False, False, True], dtype=bool)

    scene = window.Scene()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(streamlines, window.colors.misty_rose))
    if labelmask is not None:
        shape = labelmask.shape
        image_actor_z = actor.slicer(labelmask, affine)
        slicer_opacity = 0.6
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

    window.record(scene, out_path=outpath + ROIname + '_initial.png', size=(600, 600))
    if interactive:
        window.show(scene)

    colormap = actor.create_colormap(np.arange(len(clusters)))

    scene.clear()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(streamlines, window.colors.white, opacity=0.05))
    scene.add(actor.streamtube(clusters.centroids, colormap, linewidth=0.4))
    if labelmask is not None:
        image_actor_z = actor.slicer(labelmask, affine)
    window.record(scene, out_path=outpath + ROIname + '_centroids.png', size=(600, 600))
    if interactive:
        window.show(scene)

    colormap_full = np.ones((len(streamlines), 3))
    for cluster, color in zip(clusters, colormap):
        colormap_full[cluster.indices] = color

    scene.clear()
    scene.SetBackground(1, 1, 1)
    scene.add(actor.streamtube(streamlines, colormap_full))
    window.record(scene, out_path=outpath + ROIname + '_clusters.png', size=(600, 600))
    if interactive:
        window.show(scene)
