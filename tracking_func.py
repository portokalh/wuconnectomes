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
from os.path import splitext
import pickle

import os, re, struct
import glob
from dipy.tracking.utils import unique_rows

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

from dipy.segment.mask import segment_from_cfa
from dipy.segment.mask import bounding_box

import multiprocessing
# We must import this explicitly, it is not imported by the top-level
# multiprocessing module.
import multiprocessing.pool

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
import dipy.tracking.life as life
from dipy.viz import window, actor, colormap as cmap
import dipy.core.optimize as opt


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
    #This whole class allows us to run sub processes within subprocesses (subprocesses cannot exceed max count though)
    Process = NoDaemonProcess
    
def string_inclusion(string_option,allowed_strings,option_name):

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
    if string == 'any':
        return '' #if the string is any, that means that it is unspecified in filepath, therefore null
    else:
        try:
            string = str(string)
            string = string.replace(".", "_")
            return(string)
        except AttributeError:
            raise AttributeError("strfile error: not a usable number or string ")

def fix_bvals_bvecs(fbvals, fbvecs, b0_threshold=50, atol=1e-2):
    """
    Read b-values and b-vectors from disk

    Parameters
    ----------
    fbvals : str
       Full path to file with b-values. None to not read bvals.
    fbvecs : str
       Full path of file with b-vectors. None to not read bvecs.

    Returns
    -------
    bvals : array, (N,) or None
    bvecs : array, (N, 3) or None

    Notes
    -----
    Files can be either '.bvals'/'.bvecs' or '.txt' or '.npy' (containing
    arrays stored with the appropriate values).
    """

    # Loop over the provided inputs, reading each one in turn and adding them
    # to this list:
    vals = []

    from nibabel.tmpdirs import InTemporaryDirectory

    for this_fname in [fbvals, fbvecs]:
        # If the input was None or empty string, we don't read anything and
        # move on:
        if this_fname is None or not this_fname:
            vals.append(None)
        else:
            if isinstance(this_fname, str):
                base, ext = splitext(this_fname)
                if ext in ['.bvals', '.bval', '.bvecs', '.bvec', '.txt', '.eddy_rotated_bvecs', '']:
                    with open(this_fname, 'r') as f:
                        content = f.read()
                    # We replace coma and tab delimiter by space
                    with InTemporaryDirectory():
                        tmp_fname = "tmp_bvals_bvecs.txt"
                        with open(tmp_fname, 'w') as f:
                            f.write(re.sub(r'(\t|,)', ' ', content))
                        vals.append(np.squeeze(np.loadtxt(tmp_fname)))
                elif ext == '.npy':
                    vals.append(np.squeeze(np.load(this_fname)))
                else:
                    e_s = "File type %s is not recognized" % ext
                    raise ValueError(e_s)
            else:
                raise ValueError('String with full path to file is required')

    # Once out of the loop, unpack them:
    bvals, bvecs = vals[0], vals[1]

    # If bvecs is None, you can just return now w/o making more checks:
    if bvecs is None:
        return bvals, bvecs

    if bvecs.ndim != 2:
        raise IOError('bvec file should be saved as a two dimensional array')
    if bvecs.shape[1] > bvecs.shape[0]:
        bvecs = bvecs.T

    if bvecs.shape[1] == 4:
        if np.max(bvecs[:,0]) > b0_threshold:
            if bvals is None:
                bvals = bvec[0,:]
            bvecs = np.delete(bvecs,0,1)

    # If bvals is None, you don't need to check that they have the same shape:
    if bvals is None:
        return bvals, bvecs

    if len(bvals.shape) > 1:
        raise IOError('bval file should have one row')

    if max(bvals.shape) != max(bvecs.shape):
            raise IOError('b-values and b-vectors shapes do not correspond')

    from dipy.core.geometry import vector_norm

    bvecs_close_to_1 = abs(vector_norm(bvecs) - 1) <= atol

    if bvecs.shape[1] != 3:
        raise ValueError("bvecs should be (N, 3)")
    dwi_mask = bvals > b0_threshold
    if not np.all(bvecs_close_to_1[dwi_mask]):
        correctvals = [i for i,val in enumerate(bvecs_close_to_1) if val and dwi_mask[i]]
        incorrectvals = [i for i,val in enumerate(bvecs_close_to_1) if not val and dwi_mask[i]]
        baseline_bval = bvals[correctvals[0]]
        for i in incorrectvals:
            if dwi_mask[i]:
                bvecs[i,:] = bvecs[i,:] / np.sqrt(bvals[i]/baseline_bval)
        bvecs_close_to_1 = abs(vector_norm(bvecs) - 1) <= atol
        if not np.all(bvecs_close_to_1[dwi_mask]):
            incorrectvals = [i for i, val in enumerate(bvecs_close_to_1) if not val and dwi_mask[i]]
            raise ValueError("The vectors in bvecs should be unit (The tolerance "
                             "can be modified as an input parameter)")

    base, ext = splitext(fbvals)
    fbvals = base + '_fix' + ext
    np.savetxt(fbvals, bvals)
#    with open(fbvals, 'w') as f:
#        f.write(str(bvals))
#    os.remove(fbvec)
    base, ext = splitext(fbvecs)
    fbvecs = base + '_fix' + ext
    np.savetxt(fbvecs, bvecs)
#    with open(fbvecs, 'w') as f:
#        f.write(str(bvec))

    return fbvals, fbvecs

def getdwidata(mypath, subject):

    #fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    fdwi = mypath + '/nii4D_' + subject + '.nii'
    fdwi_data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    #ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'
    ffalabels = mypath + '/mask.nii.gz'
    labels, affine_labels = load_nifti(ffalabels)

    #fbvals = mypath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'
    #fbvecs = mypath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'
    fbvals = mypath + '/' + subject + '_bvals.txt'
    fbvecs = mypath + '/' + subject + '_bvec.txt'
    fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm
    
    return(fdwi_data,affine,gtab,labels)

def extractbvec_fromheader(source_file,basepath=None,save=None,verbose=True):    

    bvals = dsl = dpe = dro = None
    if save is not None:
        if basepath==None:
            basepath=os.path.split(source_file)[0]+"/"
        if os.path.isdir(basepath):
            basepath=basepath+"/"
        else:
            if os.path.isdir(os.path.split(basepath)[0]):
                basepath=basepath+"_"
            else:
                basepath=os.path.join(sourcepath,basepath)
                basepath=basepath+"_"
    with open(source_file, 'rb') as source:
        if verbose: print('INFO    : Extracting acquisition parameters')
        header_size=source.read(4)
        header_size=struct.unpack('I', header_size)
        if verbose: print('INFO    : Header size = ',int(header_size[0]))
        i=0
        stopsign=200
        for line in source: 
        
            # pattern1 = '<ParamLong."BaseResolution">  {*'
            # rx1 = re.compile(pattern1, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            
            #pattern1 = 'z_Agilent_bvalue_m00='
            pattern1 = 'z_Agilent_bvalue'
            rx1 = re.compile(pattern1, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            pattern2 = 'z_Agilent_dsl'
            rx2 = re.compile(pattern2, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            pattern3 = 'z_Agilent_dpe'
            rx3 = re.compile(pattern3, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            pattern4 = 'z_Agilent_dro'
            rx4 = re.compile(pattern4, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            i+=1
            if i==stopsign:
                print("hi")
            for a in rx1.findall(str(line)):
                bvals_all=str(line).split(',')[1]
                bvals=bvals_all.split('\\')[0]
            for a in rx2.findall(str(line)):
                dsl_all=str(line).split(',')[1]
                dsl=dsl_all.split('\\')[0]
                #dsl=([float(s) for s in dsl_all.split() if s.isnumeric()])
            for a in rx3.findall(str(line)):
                dpe_all=str(line).split(',')[1]
                dpe=dpe_all.split('\\')[0]
            for a in rx4.findall(str(line)):
                dro_all=str(line).split(',')[1]
                dro=dro_all.split('\\')[0]
    if save=="all":
        bval_file=basepath+"bvals.txt"
        File_object = open(bval_file,"w")
        File_object.write(bvals)
        File_object.close()

        dsl_file=basepath+"dsl.txt"
        File_object = open(dsl_file,"w")
        File_object.write(dsl)
        File_object.close()

        dpe_file=basepath+"dpe.txt"
        File_object = open(dpe_file,"w")
        File_object.write(bvals)
        File_object.close()
        
        dro_file=basepath+"dro.txt"
        File_object = open(dro_file,"w")
        File_object.write(bvals)
        File_object.close()      
                         
    return bvals,dsl,dpe,dro 

def prune_streamlines(streamline, mask, cutoff=2, harshcut=None, verbose=None):
    """
    gets rid of extraneous and spurious streamlines by getting rid of voxels outside of mask
    streamline: the list of streamlines
    mask: the mask usually coressponding to a brain extracted diffusion (b0) image
    cutoff: the minimum number of voxels necessary for any one streamline to be included. Must be at least 2 for
    any evaluation of the streamlines via Life
    harshcut: if set at True, harshcut will delete any voxels of a streamline outside of a mask as well as any following
    voxel, if set to None or False, harshcut will only delete voxels outside of mask and treat subsequent voxels as normal

    returns the pruned streamline
    """
    delstream=[]
    if verbose:
        print("Starting the pruning")
    duration = time()
    if mask is not None:
        for idx,s in enumerate(streamline): #iterate through all streams
            s_vox=unique_rows(np.round(s).astype(np.intp))
            cutlist=[]
            for vox in range(np.shape(s_vox)[0]):
                if not mask[tuple(s_vox[vox])]:
                    cutlist.append(vox)         #if local mask is 0, add voxel to list of voxels to cut
            if harshcut:                        #if harshcut, destroy all voxels folloxing the out of mask voxel in streamline
                startcut = np.min(cutlist)
                cutlist = range(startcut, len(np.shape(s_vox)[0]))
                s = np.delete(s, cutlist)
            else:
                s = np.delete(s, cutlist, axis=0)   #else, cut all voxels of streamlines not in mask
            streamline[idx] = s                       #replace old streamline with new spurious streamline
            streamshape = np.shape(np.asarray(s))
            if streamshape[0] < cutoff: #minimum number of voxels required in streamline, default and minimum at 2
                delstream.append(idx)       #if number of voxels in streamline too low, cut it
                if verbose == "oververbose":
                    print("Skipped stream " + str(idx) + " out of " + str(len(streamline)) + " streamlines for being too small")
    else:
        for idx, s in enumerate(streamline):
            streamshape = np.shape(np.asarray(s))
            if streamshape[0] < cutoff:     #minimum number of voxels required in streamline, default and minimum at 2
                delstream.append(idx)
                if verbose:
                    print("Skipped stream " + str(idx) + " out of " + str(len(streamline)) + " streamlines for being too small")
    if verbose:
        print("Obtaining fiber signal process done in " + str(time() - duration) + "s")
        if len(delstream) > 0:
            print("Skipped " + str(len(delstream)) + " out of " + str(
                len(streamline)) + " due to size constraints (tiny streamlines)")
    for idx in reversed(delstream):
        streamline.pop(idx)
    return streamline

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



def unload_trk(tractogram_path,reference='same'):
    """ Similar functionality as the older version of load_trk, as it directly
    extracts the streams and header instead of returning a Tractogram object
    
    Parameters
    ----------
    tractogram_path: the file path of the tractogram data ( path/tract.trk )
    reference: the file used for the header information. if 'same', use the hdr from tractogram file
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


def testsnr():

    corpus_mask = np.where(labels == 121, 1, 0) + np.where(labels == 1121, 1, 0)
    #meed to change threshold, possibly ROI, that better corresponds to mice (pick area with high FA)
    #probably should investigate area with
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
    from dipy.tracking import utils

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

def LiFEvaluation(dwidata, trk_streamlines, gtab, header=None, roimask=None, affine=None, display = True,
                  outpathfig = None, outpathtrk = None, processes = 1, interactive = True, verbose = None):

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
        roimask = dwidata>0

    fiber_model = life.FiberModel(gtab)
    # inv_affine must be used if the streamlines are in the world space, and thus we must useapply the inverse affine of dwi
    #when comparing the diffusion directions gtab and the voxels of trk
    #inv_affine = np.linalg.inv(hardi_img.affine)

    #fiber_fit will fit the streamlines to the original diffusion data and
    fiber_fit = fiber_model.fit(dwidata, trk_streamlines, affine=np.eye(4), processes=processes, verbose=verbose)
    optimized_sl = list(np.array(trk_streamlines)[np.where(fiber_fit.beta > 0)[0]])
    plt.ioff()

    fig, ax = plt.subplots(1)
    ax.hist(fiber_fit.beta, bins=100, histtype='step')
    ax.set_xlabel('Fiber weights')
    ax.set_ylabel('# fibers')
    ROI_actor = actor.contour_from_roi(roimask, color=(1., 1., 0.),
                                          opacity=0.5)

    if outpathfig is not None:
        fig.savefig(outpathfig+'_beta_histogram.png')
    """
    vol_actor = actor.slicer(t1_data)
    
    vol_actor.display(x=40)
    vol_actor2 = vol_actor.copy()
    vol_actor2.display(z=35)        
    
    """


    if display:
        plt.show()

    ren = window.Renderer()
    ren.add(actor.streamtube(optimized_sl, cmap.line_colors(optimized_sl)))
    ren.add(ROI_actor)
    #ren.add(vol_actor)

    if outpathfig is not None:
        window.record(ren, n_frames=1, out_path=outpathfig +'_life_optimized.png',
                      size=(800, 800))
    if outpathtrk is not None:
        outpathfile = str(outpathtrk) + "_lifeopt.trk"
        try:
            myheader = create_tractogram_header(outpathfile, *header)
            optimized_sl_gen = lambda: (s for s in optimized_sl)
            save_trk_heavy_duty(outpathfile, streamlines=optimized_sl_gen,
                                affine=affine, header=myheader)
        except TypeError:
            print('Could not save new tractogram file, header of original trk file not properly implemented into '
                  'LifEvaluation')
    if interactive:
        window.show(ren)

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
    if outpathfig is not None:
        fig.savefig(outpathfig + '_error_histograms.png')

    if display:
        plt.show()

    vol_model = np.ones(dwidata.shape[:3]) * np.nan
    vol_model[fiber_fit.vox_coords[:, 0],
              fiber_fit.vox_coords[:, 1],
              fiber_fit.vox_coords[:, 2]] = model_rmse
    vol_mean = np.ones(dwidata.shape[:3]) * np.nan
    vol_mean[fiber_fit.vox_coords[:, 0],
             fiber_fit.vox_coords[:, 1],
             fiber_fit.vox_coords[:, 2]] = mean_rmse
    vol_improve = np.ones(dwidata.shape[:3]) * np.nan
    vol_improve[fiber_fit.vox_coords[:, 0],
                fiber_fit.vox_coords[:, 1],
                fiber_fit.vox_coords[:, 2]] = mean_rmse - model_rmse
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

    if display:
        plt.show()
        """
        
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
        if savefig:
            fig.savefig( outpath + "/spatial_errors.png")        
        
        """
    return model_error, mean_error


def denoise_pick(data,affine,hdr,outpath,mask,type_denoise='macenko', processes = 1, savedenoise= True, verbose=False, display=None):

    allowed_strings=['mpca','yes','all','gibbs','none', 'macenko']
    string_inclusion(type_denoise, allowed_strings, "type_denoise")

    if type_denoise == 'macenko' or 'mpca' or type_denoise == 'yes' or type_denoise == 'all':
        #data, snr = marcenko_denoise(data, False, verbose=verbose)
        t = time()
        denoised_arr, sigma = mppca(data, patch_radius=2, return_sigma=True, processes=processes, verbose=verbose)
        outpath_mpca = outpath + '_mpca_paralleltest_2.nii.gz'
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
            print("Time taken for the gibbs removal " -t + time())
        if display:
            denoise_fig(data,data_corrected,type='gibbs')

        data=data_corrected
        
    if type_denoise == 'none':
        print('No denoising was done')

    return data

def make_tensorfit(data,mask,gtab,outpath,verbose=None):


    from dipy.reconst.dti import TensorModel

    if verbose:
        print('Calculating the tensor model from bval/bvec values of ', subject)
    tensor_model = TensorModel(gtab)

    t1 = time()
    tensor_fit = tensor_model.fit(data, mask)

    duration1 = time() - t1
    if verbose:
        print(subject + ' DTI duration %.3f' % (duration1,))

    outpathbmfa = outpath + 'bmfa' + subject + '.nii.gz'
    try:
        save_nifti(outpathbmfa, tensor_fit.fa, affine)
        if verbose:
            print('Saving subject'+ subject+ ' at ' + outpathbmfa)
    except:
        print('Warning: Failutre to save the fa as nifti')
        outpathbmfa = None

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
        print("Starting calculation of Constant solid angle model for subject " + subject)

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
        print('Computing classifier for local tracking for subject ' + subject)
    classifier = BinaryStoppingCriterion(mask)
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
        outpathtrk = outpathsubject + "_bmCSA_detr_small_" + stringstep + "_v2.trk"
        save_trk_heavy_duty(outpathtrk, streamlines=sg_small,
                            affine=affine, header=trkheader,
                            shape=mask.shape, vox_size=vox_size)
    else:
        outpathtrk = None
    if saved_tracts == "large" or saved_tracts == "both" or saved_tracts == "all":
        sg = lambda: (s for s in streamlines_generator)
        outpathtrk = outpathsubject + "bmCSA_detr_all" + stringstep + ".trk"
        save_trk_heavy_duty(outpathtrk, streamlines=sg,
                            affine=affine, header=trkheader,
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

    return outpathtrk

def dwi_preprocessing(mypath,outpath,subject,denoise="none",savedenoise=True,savefa="yes",processes=1, verbose = False):

    if os.path.isdir(mypath):
        fdwi = glob.glob(mypath + '/' + subject + '*.nii*')[0]
    elif os.path.isfile(mypath):
        fdwi = mypath
    else:
        raise(AttributeError)

    fdwi_img = nib.load(fdwi)
    fdwi_data = fdwi_img.get_data()
    affine = fdwi_img.affine
    hdr = fdwi_img.header
    #fdwi_data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    labels = glob.glob(mypath + '/' + subject + '*labels*nii.gz') #ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    mask = glob.glob(mypath + '/' + subject + '*mask*nii.gz')
    labelmask = labels + mask #Whether the image is a labels or a mask, we output a mask for mpca
    try:
        labelmask, affine_labmask = load_nifti(labelmask[0])
        # Build Brain Mask
        bm = np.where(labelmask == 0, False, True)
        mask = bm
        if verbose:
            print("Mask detected, proceeding as normal")
    except IndexError:
        print("Warning: no mask detected (not recommended)")
        mask = np.ones(fdwi_data.shape)

    if os.path.isdir(mypath):
        fbvals = glob.glob(mypath + '/' + subject + '*bval*')[0]
        fbvecs = glob.glob(mypath + '/' + subject + '*bvec*')[0]
    else:
        print("Error, havent implemented multiple bvals bvec options")

    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    if verbose:
        print('Running the ' + subject + ' file')

    #preprocessing section (still have to test denoising functions)
    #data = denoise_pick(data, mask, 'macenko', display=None) #data = denoise_pick(data, mask, 'gibbs', display=None)
    #fdwi_data = denoise_pick(fdwi_data, mask, 'all', display=None)
    outpathdenoise= outpath + subject + '_nii4D_RAS'
    fdwi_data = denoise_pick(fdwi_data, affine,hdr, outpathdenoise, mask, denoise, savedenoise=savedenoise, processes=processes, verbose=verbose) #accepts mpca, gibbs, all, none

    #testsnr => not yet fully implemented
    if savefa == "yes" or savefa == "y" or savefa == 1 or savefa is True or savefa == "all":
        outpathbmfa = make_tensorfit(fdwi_data,mask,gtab,verbose = verbose)
    else:
        print('FA was not calculated')
        outpathbmfa=None

def dwi_create_tracts(mypath,outpath,subject,step_size,peak_processes,saved_tracts="small",save_fa="yes",
                      denoise="none",verbose=None):

    if verbose:
        print('Running the ' + subject + ' file')

    """
    fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    fdwi = mypath + '/nii4D_' + subject + '.nii'
    fdwi_data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'
    ffalabels = mypath + '/mask.nii.gz'
    labels, affine_labels = load_nifti(ffalabels)

    fbvals = mypath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'
    fbvecs = mypath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'
    fbvals = mypath + '/' + subject + '_bvals.txt'
    fbvecs = mypath + '/' + subject + '_bvec.txt'
    fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm
    """

    fdwi, affine, gtab, labels = getdwidata(dwipath, subject)
    #preprocessing section (still have to test denoising functions)
    #data = denoise_pick(data, mask, 'macenko', display=None) #data = denoise_pick(data, mask, 'gibbs', display=None)
    #fdwi_data = denoise_pick(fdwi_data, mask, 'all', display=None)
    #fdwi_data = denoise_pick(fdwi_data, mask, denoise, verbose) #accepts mpca, gibbs, all, none
    #testsnr => not yet fully implemented
    if save_fa == "yes" or save_fa == "y" or save_fa == 1 or save_fa is True or save_fa == "all":
        outpathbmfa = make_tensorfit(fdwi_data,mask,gtab,verbose=verbose)
    else:
        print('FA was not calculated')
        outpathbmfa=None

    allowed_strings=["small","large","all","both","none"]
    string_inclusion(saved_tracts, allowed_strings, "saved_tracts")

    outpathsubject = outpath + subject

    trkheader = create_tractogram_header("place.trk", *get_reference_info(fdwi))

    #if multishell_split: #possible idea to creatr tracts from either one bval or another before doing it on all
    outpathtrk = QCSA_tractmake(fdwi_data,affine,vox_size,gtab,mask,trkheader,step_size,peak_processes,outpathsubject,saved_tracts=saved_tracts,verbose=verbose,subject=subject)
    
    return subject, outpathbmfa, outpathtrk

def evaluate_tracts(dwipath,trkpath,subject,stepsize, tractsize, outpathfig=None, processes=1, doprune=True,
                    display=True, verbose=None):

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
    dwidata, affine, gtab, labels = getdwidata(dwipath, subject)
    roimask = (labels == 163) + (labels == 1163)  # + (labels == 120) + (labels == 1120)
    roimask = labels > 1
    outpathfig = outpathfig+'/'+subject
    print("Beginning Tract Evaluation of " + subject)
    stepsize = strfile(stepsize)
    trkpaths = glob.glob(trkpath+'/'+subject+'*'+tractsize+'*'+stepsize+'*.trk')
    trkfile = trkpaths[0]
    outpathtrk = trkpath+'/'+subject
    if len(trkpaths) > 1:
        print("Warning, multiple files detected, only applying pipeline on " + trkfile)
    trkdata = load_trk(trkfile, 'same')
    trkdata.to_vox()
    header = trkdata.space_attribute
    trkstreamlines = trkdata.streamlines

    duration=time()

    """
    ministream = []
    for idx, stream in enumerate(trkstreamlines):
        if (idx % 1000) == 0:
            ministream.append(stream)
    trkstreamlines = ministream
    """

    if doprune:
        cutoff = 2
        if np.min(dwidata[:, :, :, 0]) == 0:
            mask = (dwidata[:, :, :, 0] > 0)
        else:
            mask = None
        trkstreamlines=prune_streamlines(list(trkstreamlines), mask, cutoff=cutoff, verbose=verbose)
        #trkstreamlines=prune_streamlines(trkstreamlines, mask, cutoff=cutoff, verbose=verbose)


    #if display:
    #    window_show_test(trkstreamlines, mask_roi, anat, interactive= True, outpath=None)


    #trk_streamlines = [s[0] for s in nib.trackvis.read(trkfile, points_space='voxel')[0]]
    #len_sl=len(trk_streamlines)
    #% memit fiber_fit = fiber_model.fit(data, trk_streamlines[2 ** 12], affine=np.eye(4))
    interactive = False
    model_error, mean_error = LiFEvaluation(dwidata, trkstreamlines, gtab, header=header, roimask=roimask,
                                                     affine=affine,display=display, outpathfig=outpathfig,
                                                     outpathtrk=outpathtrk, processes=processes,
                                                     interactive=interactive, verbose=verbose)
    picklepath = '/Users/alex/jacques/test_pickle_subj'+subject+'.p'
    results=[outpathtrk,model_error,mean_error]
    if subject == "N54859":
        pickle.dump(results, open(picklepath, "wb"))

    print("Finished life evaluation of subject " + (subject)+ ", whole process took " + str(time()-duration) + " s")
    #picklepath = trkpath+subject+'lifevals.p'
    #pickle.dump(tracteval_results, open(picklepath,"wb"))
    return [outpathfig, model_error, mean_error]
"""
def create_tracts(mypath,outpath,subject,step_size,peak_processes=1,saved_tracts="small",save_fa="yes",denoise="mpca",verbose=None):

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

    try:
        denoise=denoise.lower()
    except AttributeError:
        pass




    if denoise == 'mpca' or denoise == 'yes' or denoise == 'all':
        #data, snr = marcenko_denoise(data, False, verbose=verbose)
        t = time()
        denoised_arr, sigma = mppca(data, patch_radius=2, return_sigma=True)

        mean_sigma = np.mean(sigma[mask])
        b0 = denoised_arr[..., 0]

        mean_signal = np.mean(b0[mask])

        snr = mean_signal / mean_sigma

        if verbose:
            print("Time taken for local MP-PCA ", -t +
                  time())
            print("The SNR of the b0 image appears to be at " + str(snr))
        if display:
            marcenko_denoise_fig(data, denoised_arr, 'None')

        data=denoised_arr

    if denoise == 'gibbs' or denoise =='all':
        data_corrected = gibbs_removal(data_slices, slice_axis=2)

        data=data_corrected

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm

    sphere = get_sphere('repulsion724')

    from dipy.reconst.dti import TensorModel

    if verbose:
        print('Calculating the tensor model from bval/bvec values of ', subject)
    tensor_model = TensorModel(gtab)

    t1 = time()
    #tensor_fit = tensor_model.fit(data, mask)
    import pickle
    picklepath = '/Users/alex/jacques/tensor4589.p'
    #pickle.dump(tensor_fit, open(picklepath, "wb"))
    tensor_fit = pickle.load(open(picklepath, "rb"))
    testsnr=False
    if testsnr:
        corpus_mask = np.where(labels == 121, 1, 0) + np.where(labels == 1121, 1, 0)
        #meed to change threshold, possibly ROI, that better corresponds to mice (pick area with high FA)
        #probably should investigate area with
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

    try:
        save_fa=save_fa.lower()
    except AttributeError:
        pass
    if save_fa == "yes" or save_fa == "y" or save_fa == 1 or save_fa is True or save_fa == "all":
        outpathbmfa = outpath + 'bmfa' + subject + '.nii.gz'
        save_nifti(outpathbmfa, tensor_fit.fa, affine)
        if verbose:
            print('Saving subject'+ subject+ ' at ' + outpathbmfa)
    else:
        outpathbmfa = None
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
        outpathtrk = outpath + subject + "_bmCSA_detr_small_" + stringstep + "_v3.trk"
        myheader = create_tractogram_header(outpathtrk, *get_reference_info(fdwi))
        save_trk_heavy_duty(outpathtrk, streamlines=sg_small,
                affine=affine, header=myheader,
                shape=mask.shape, vox_size=vox_size)
    else:
        outpathtrk = None
    if saved_tracts == "large" or saved_tracts == "both" or saved_tracts == "all":
        sg = lambda: (s for s in streamlines_generator)
        outpathtrk = outpath+subject+"bmCSA_detr_all"+stringstep+"_v1.trk"
        myheader = create_tractogram_header(outpathtrk,*get_reference_info(fdwi))
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

    return subject, outpathtrk, outpathbmfa
"""