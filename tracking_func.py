#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:38:04 2019

@author: Jacques Stout
File to add any repeated python functions in tractography processing
(Might be eventually deleted and combined with previous SAMBA function files)
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import nibabel as nib
import numpy as np
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from dipy.io.streamline import load_trk

from time import time
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.local_tracking import LocalTracking
from dipy.direction import peaks
import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes
from nibabel.streamlines import detect_format
from dipy.io.utils import (create_tractogram_header,
                           get_reference_info)
from dipy.viz import window, actor

import multiprocessing
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
import multiprocessing.pool

from random import randint


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


def save_trk_heavy_duty(fname, streamlines, affine, vox_size=None, shape=None, header=None):
    """ Saves tractogram files (*.trk)

    Parameters
    ----------
    fname : str
        output trk filename
    streamlines : list of 2D arrays, generator or ArraySequence
        Each 2D array represents a sequence of 3D points (points, 3).
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline points.
    vox_size : array_like (3,), optional
        The sizes of the voxels in the reference image (default: None)
    shape : array, shape (dim,), optional
        The shape of the reference image (default: None)
    header : dict, optional
        Metadata associated to the tractogram file(*.trk). (default: None)
    """
    if vox_size is not None and shape is not None:
        if not isinstance(header, dict):
            header = {}
        header[Field.VOXEL_TO_RASMM] = affine.copy()
        header[Field.VOXEL_SIZES] = vox_size
        header[Field.DIMENSIONS] = shape
        header[Field.VOXEL_ORDER] = "".join(aff2axcodes(affine))

    tractogram = nib.streamlines.LazyTractogram(streamlines)
    tractogram.affine_to_rasmm = affine
    trk_file = nib.streamlines.TrkFile(tractogram, header=header)
    nib.streamlines.save(trk_file, fname)


def unload_trk(tractogram_path,reference='same'):
    """ Similar functionality as the older version of load_trk, as it directly
    extracts the streams and header instead of returning a Tractogram object
    
    Parameters
    ----------
    tractogram_path: the file path of the tractogram data ( path/tract.trk )
    """
    
    if reference.lower()=="same":
        print("Reference taken directly from file")
    tract_obj=load_trk(tractogram_path,reference)
    streams_control=tract_obj.streamlines
    hdr_control=tract_obj.space_attribute
    return streams_control,hdr_control,tract_obj
    

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


def create_tracts(mypath,outpath,subject,step_size,peak_processes,saved_tracts="small",save_fa="yes",verbose=None):

    fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'

    ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    fbvals = mypath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'

    fbvecs = mypath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'

    labels, affine_labels = load_nifti(ffalabels)

    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    allowed_strings=["small","large","all","both","none"]
    try:
        saved_tracts=saved_tracts.lower()
    except AttributeError:
        pass

    if not any(saved_tracts == x for x in allowed_strings):
        raise ValueError("Unrecognized string, please check your input for 'saved tracts' ")
    if saved_tracts == "None" or saved_tracts is None:
        raise Warning("Saved_tracts stated as None value, no tracts will be saved for this run")
    if verbose:
        print('Running the ' + subject + ' file')
    # Correct flipping issue
    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm

    sphere = get_sphere('repulsion724')

    from dipy.reconst.dti import TensorModel

    if verbose:
        print('Calculating the tensor model from bval/bvec values of ', subject)
    tensor_model = TensorModel(gtab)

    t1 = time()
    tensor_fit = tensor_model.fit(data, mask)
    #    save_nifti('bmfa.nii.gz', tensor_fit.fa, affine)
    #   wenlin make this change-adress name to each animal
    try:
        save_fa=save_fa.lower()
    except AttributeError:
        pass
    if save_fa == "yes" or save_fa == "y" or save_fa == 1 or save_fa is True or save_fa == "all":
        outpathbmfa = outpath + 'bmfa' + subject + '.nii.gz'

    if verbose:
        print('Saving subject'+ subject+ ' at ' + outpathbmfa)
    save_nifti(outpathbmfa, tensor_fit.fa, affine)
    fa = tensor_fit.fa
    duration1 = time() - t1
    # wenlin make this change-adress name to each animal
    #    print('DTI duration %.3f' % (duration1,))
    if verbose:
        print(subject + ' DTI duration %.3f' % (duration1,))

    # Compute odfs in Brain Mask
    t2 = time()

    csa_model = CsaOdfModel(gtab, 6)
    if peak_processes < 2:
        parallel=False
    else:
        parallel=True
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
        print('Computing classifier for local tracking')
    classifier = BinaryStoppingCriterion(bm)
    from dipy.tracking import utils

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

    stringstep = str(step_size)
    stringstep = "_" + stringstep.replace(".", "_")
    # stringstep=""
    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)

    # the function above will bring all streamlines in memory
    # streamlines = Streamlines(streamlines_generator)

    # save a smaller part by only keeping one in 10 streamlines

    if saved_tracts == "small" or saved_tracts == "both":
        sg_small = lambda: (s for i, s in enumerate(streamlines_generator) if i % 10 == 0)
        outpathtrk = outpath + subject + "_bmCSA_detr_small" + stringstep + ".trk"
        myheader = create_tractogram_header(outpathtrk, *get_reference_info(fdwi))
        save_trk_heavy_duty(outpathtrk, streamlines=sg_small,
                affine=affine, header=myheader,
                shape=mask.shape, vox_size=vox_size)
    if saved_tracts == "large" or saved_tracts == "both" or saved_tracts == "all":
        sg = lambda: (s for s in streamlines_generator)
        outpathtrk=outpath+subject+"bmCSA_detr"+stringstep+".trk"
        myheader=create_tractogram_header(outpathtrk,*get_reference_info(fdwi))
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                affine=affine, header=myheader,
                shape=mask.shape, vox_size=vox_size)
    if saved_tracts == "none" or saved_tracts is None:
        print("Tract files were not saved")



    # save everything - will generate a 20+ GBytes of data - hard to manipulate

    # possibly add parameter in csv file or other to decide whether to save large tractogram file
    # outpathfile=outpath+subject+"bmCSA_detr"+stringstep+".trk"
    # myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))

    # save_trk_heavy_duty(outpathfile, streamlines=sg_small,
    #                    affine=affine, header=myheader,
    #                    shape=mask.shape, vox_size=vox_size)

    duration3 = time() - t3
    if verbose:
        print(duration3)
    # wenlin make this change-adress name to each animal
    #    print('Tracking duration %.3f' % (duration3,))
    if verbose:
        print(subject + ' Tracking duration %.3f' % (duration3,))

    return(subject,outpathtrk,outpathbmfa)
    
def almicedf_fix(df,verbose=None):
    #masterFile='/Users/alex/AlexBadea_MyPapers/FIGURES/mwm/mwm_master_organized.csv'

    df=df.replace({'runno': {"N54716/N54915" : "N54915", "N54714/N54916": "N54916", 
                                "N54715/N54917" : "N54917", "N54718//N54721" : "N54718", 
                                "N54760/N54754/N54755" : "N54760", "N54757/N54759" : "N54759",
                                "N54805/N54800" : "N54805", "N54891/N54900LRspecific" : "N54891"}})
    df=df.replace({'Day': {"ProbeTrial1": "ProbeTrial", "ProbeTrial2": "ProbeTrial"}})
    df=df.dropna(subset=['runno'])
    
    alldays=df.Day.unique()
    df=df.dropna()
    
    if verbose:
        for day in alldays:
            df_day=df.loc[df['Day'] == day]
            sf_day=df_day.groupby('runno')['Acq Number'].nunique()
            print('For the Day: '+ day + ' those animals had less than the standard number of tests:')
            print(sf_day.where(sf_day<np.max(sf_day)).dropna())
    
    
    return(df)
    #na_columns=df.columns[df.isna().any()]
    #df_withna=df[df.isnull().any(axis=1)][:].head()    
    #df = mice_day2.groupby('runno')['Acq Number'].nunique()
        
