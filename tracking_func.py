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
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
import pickle

from types import ModuleType, FunctionType
from gc import get_referents

import smtplib 

import os, re, sys, io, struct, socket, datetime
from email.mime.text import MIMEText
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
#import dipy.tracking.life as life
import JSdipy.tracking.life as life
from dipy.viz import window, actor, colormap as cmap
import dipy.core.optimize as opt
from functools import wraps

from tractfile_handler import target, prune_streamlines
from bvec_handler import fix_bvals_bvecs#, extractbvec_fromheader
from figures_handler import denoise_fig, show_bundles, window_show_test, LifEcreate_fig
from tract_eval import bundle_coherence, LiFEvaluation
from dif_to_trk import make_tensorfit, QCSA_tractmake
from daemon_tools import MyPool

mylogin = "jas297"

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
    #Converts strings into more usable 'file strings (mostly takes out . and turns it into _
    if string == 'any':
        return '' #if the string is any, that means that it is unspecified in filepath, therefore null
    else:
        try:
            string = str(string)
            string = string.replace(".", "_")
            return(string)
        except AttributeError:
            raise AttributeError("strfile error: not a usable number or string ")

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    #sum size of object & all members within that object.
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


def send_mail(msg_txt,subject="Cluster message"):
    #Send mail with specified txt (and subject) to default address specified in global variable
    msg_content = io.StringIO()
    msg_content.write("Datetime : %s \n\n" % datetime.datetime.now())
    msg_content.write("JobID : %d \n\n" % os.getpid() )
    msg_content.write("Message : %s \n\n" % msg_txt)
    msg = MIMEText(msg_content.getvalue())
    msg_content.close()
    to_addr = "%s@duke.edu" % mylogin
    from_addr = "%s@duke.edu" % mylogin
    msg['Subject'] = subject
    msg['from'] = from_addr
    msg['to'] = to_addr
    s = smtplib.SMTP('smtpgw.duhs.duke.edu')
    s.sendmail(from_addr, [to_addr], msg.as_string())
    s.quit()

def getdwidata(mypath, subject, verbose=None):

    """



    """

    #fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    #fdwipath = mypath + '/nii4D_' + subject + '.nii'
    if os.path.exists(mypath + '/Reg_' + subject + '_nii4D.nii.gz'):
        fdwipath = mypath + '/Reg_' + subject + '_nii4D.nii.gz'
    elif os.path.exists(mypath + '/nii4D_' + subject + '.nii'):
        fdwipath = mypath + '/nii4D_' + subject + '.nii'
    elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'
    #fdwi_data, affine, vox_size = load_nifti(fdwipath, return_voxsize=True)

    if verbose:
        txt="Extracting information from the dwi file located at "+ fdwipath
        print(txt)
        send_mail(txt,subject="Begin data extraction")

    img = nib.load(fdwipath)
    fdwi_data = img.get_data()
    vox_size = img.header.get_zooms()[:3]
    affine = img.affine
    hdr = img.header
    del(img)

    #ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    if os.path.exists(mypath+'/Reg_'+subject+'_nii4D_brain_mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath+'/Reg_'+subject+'_nii4D_brain_mask.nii.gz')
    elif os.path.exists(mypath+'/'+subject+'_chass_symmetric2_labels_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath+'/'+subject+'_chass_symmetric2_labels_RAS.nii.gz')
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii.gz')
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject +'_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject +'_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz')
    else:
        print('mask not found, taking all non null values in nii file instead (not recommended for complex operations)')
        labels = np.ones(fdwi_data.shape[0:3])
        affine_labels=affine
        

    #fbvals = mypath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'
    #fbvecs = mypath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'
    fbvals = glob.glob(mypath + '*/*' + subject + '*_bval*.txt')[0]
    fbvecs = glob.glob(mypath + '*/*' + subject + '*_bvec*.txt')[0]
    #fbvecs = mypath + '/' + subject + '_bvec.txt'
    fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    #bvecs = np.c_[bvecs[:, 0], -bvecs[:, 1], bvecs[:, 2]]  # FOR RAS according to Alex
    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]] #FOR RAS

    #bvecs = np.c_[-bvecs[:, 0], bvecs[:, 1], bvecs[:, 2]] #FOR ARI (original form of C57)

    #bvecs = np.c_[bvecs[:, 1], bvecs[:, 0], -bvecs[:, 2]]
    #bvecs = np.c_[-bvecs[:, 1], bvecs[:, 0], bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm
    
    return(fdwi_data,affine,gtab,labels,vox_size, fdwipath, hdr)

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

def denoise_pick(data,affine,hdr,outpath,mask,type_denoise='macenko', processes = 1, savedenoise= True, verbose=False, display=None):

    allowed_strings=['mpca','yes','all','gibbs','none', 'macenko']
    string_inclusion(type_denoise, allowed_strings, "type_denoise")

    if type_denoise == 'macenko' or 'mpca' or type_denoise == 'yes' or type_denoise == 'all':
        #data, snr = marcenko_denoise(data, False, verbose=verbose)
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
            print("Time taken for the gibbs removal " -t + time())
        if display:
            denoise_fig(data,data_corrected,type='gibbs')

        data=data_corrected
        
    if type_denoise == 'none':
        print('No denoising was done')

    return data


def dwi_preprocessing(dwipath,outpath,subject,denoise="none",savefa="yes",processes=1, labelslist=None, strproperty="", verbose = False):

    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, hdr = getdwidata(dwipath, subject, verbose)
    #fdwi_data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    #labels = glob.glob(mypath + '/' + subject + '*labels*nii.gz') #ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    #mask = glob.glob(mypath + '/' + subject + '*mask*nii.gz')
    #labelmask = labels + mask #Whether the image is a labels or a mask, we output a mask for mpca
        
    # Build Brain Mask
    
    if labelslist is None:
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

    mask= np.zeros(np.shape(mask),dtype=bool)
    mask[90:110,90:110,90:110]
    fdwi_data = fdwi_data * np.repeat(mask[:, :, :, None], np.shape(fdwi_data)[3], axis=3)

    outpathdenoise= outpath + subject + '_nii4D_RAS'
    #fdwi_data = denoise_pick(fdwi_data, affine,hdr, outpathdenoise, mask, denoise, processes=processes, verbose=verbose) #accepts mpca, gibbs, all, none

    #testsnr => not yet fully implemented
    print(savefa)
    if savefa == "yes" or savefa == "y" or savefa == 1 or savefa is True or savefa == "all":
        outpathbmfa = make_tensorfit(fdwi_data,mask,gtab,affine,subject,outpath=outpath,strproperty=strproperty,verbose = verbose)
    else:
        print('FA was not calculated')
        outpathbmfa=None

def create_tracts(dwipath,outpath,subject,step_size,peak_processes,strproperty="",saved_tracts="small",save_fa="yes",
                      labelslist = None, verbose=None):

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

    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, _ = getdwidata(dwipath, subject)

    if labelslist is None:
        mask = np.where(labelmask == 0, False, True)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at " + dwipath + " for subject " + subject)
        mask = np.zeros(np.shape(labelmask), dtype=int)
        for label in labelslist:
            mask = mask + (labelmask == label)
        fdwi_data = fdwi_data * np.repeat(mask[:, :, :, None], np.shape(fdwi_data)[3], axis=3)

    #preprocessing section (still have to test denoising functions)
    #data = denoise_pick(data, mask, 'macenko', display=None) #data = denoise_pick(data, mask, 'gibbs', display=None)
    #fdwi_data = denoise_pick(fdwi_data, mask, 'all', display=None)
    #fdwi_data = denoise_pick(fdwi_data, mask, denoise, verbose) #accepts mpca, gibbs, all, none
    #testsnr => not yet fully implemented


    if save_fa == "yes" or save_fa == "y" or save_fa == 1 or save_fa is True or save_fa == "all":
        outpathbmfa = make_tensorfit(fdwi_data,mask,gtab,outpath=outpath,strproperty=strproperty,verbose=verbose)
    else:
        print('FA was not calculated')
        outpathbmfa=None

    allowed_strings=["small","large","all","both","none"]
    string_inclusion(saved_tracts, allowed_strings, "saved_tracts")
    outpathsubject = outpath + subject + strproperty

    trkheader = create_tractogram_header("place.trk", *get_reference_info(fdwipath))
    #if multishell_split: #possible idea to creatr tracts from either one bval or another before doing it on all
    print(verbose)
    if verbose:
        print("The QCSA Tractmake is ready to launch for subject " + subject)
        send_mail("The QCSA Tractmake is ready to launch for subject " + subject,subject="QCSA main function start")
        print("email sent")
        #headers="From: %s\r\nTo: %s\r\nSubject:QCSA start\r\n\r\n" % (useremail,useremail)
        #text="""The QCSA Tractmake is ready to launch for subject %s""" % (subject)  
        #message=headers+text 
        #mailServer=smtplib.SMTP(serverURL) 
        #mailServer.sendmail(useremail,useremail,message) 
        #mailServer.quit() 
    outpathtrk = QCSA_tractmake(fdwi_data,affine,vox_size,gtab,mask,trkheader,step_size,peak_processes,outpathsubject,saved_tracts=saved_tracts,verbose=verbose,subject=subject)
    
    return subject, outpathbmfa, outpathtrk

def evaluate_tracts(dwipath,trkpath,subject,stepsize, tractsize, labelslist=None, outpathpickle=None, outpathfig=None, processes=1, doprune=True,
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
    fdwi_data, affine, gtab, labelmask, vox_size, fdwipath, _ = getdwidata(dwipath, subject)
    
    if labelslist is None:
        roimask = (fdwi_data[:, :, :, 0] > 0)
    else:
        if labelmask is None:
            raise ("File not found error: labels requested but labels file could not be found at "+dwipath+ " for subject " + subject)
        roimask = np.zeros(np.shape(labelmask),dtype=int)
        for label in labelslist:
            roimask = roimask + (labelmask == label)
    #(labels == 163) + (labels == 1163) + (labels == 120) + (labels == 1120)
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
    trkpaths = glob.glob(trkpath+'/'+subject+'*'+tractsize+'*'+stepsize+'*.trk')
    trkfile = trkpaths[0]
    outpathtrk = trkpath+'/'
    if len(trkpaths) > 1:
        print("Warning, multiple files detected, only applying pipeline on " + trkfile)
    trkdata = load_trk(trkfile, 'same')
    trkdata.to_vox()
    if hasattr(trkdata, 'space_attribute'):
        header = trkdata.space_attribute
    elif hasattr(trkdata, 'space_attributes'):
        header = trkdata.space_attributes
    trkstreamlines = trkdata.streamlines

    ratio = 1

    if verbose:
        txt="Beginning the evaluation of subject "+ subject + "\n We will evaluate for one in "+str(ratio)+" streamlines"
        print(txt)
        send_mail(txt,subject="LifE start msg ")

    if ratio != 1:
        ministream = []
        for idx, stream in enumerate(trkstreamlines):
            if (idx % ratio) == 0:
                ministream.append(stream)
        trkstreamlines = ministream

    doprune=False
    cutoff = 2
    if doprune:
        trkstreamlines=prune_streamlines(list(trkstreamlines), roimask, cutoff=cutoff, verbose=verbose)
        #fdwi_data = fdwi_data * np.repeat(roimask[:, :, :, None], np.shape(fdwi_data)[3], axis=3)
        #trkstreamlines=prune_streamlines(trkstreamlines, mask, cutoff=cutoff, verbose=verbose)
    else:
        affine=np.eye(4)
        trktempstreamlines = target(trkstreamlines, affine, roimask, include=True, strict="longstring")
        trkstreamlines = Streamlines(trktempstreamlines)
        trkstreamlines = prune_streamlines(list(trkstreamlines), None, cutoff=cutoff, verbose=verbose)
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
            #outpathfile = str(outpathtrk) + "_fimbria.trk"
            #myheader = create_tractogram_header(outpathfile, *header)
            #optimized_sl_gen = lambda: (s for s in trkstreamlines)
            #save_trk_heavy_duty(outpathfile, streamlines=optimized_sl_gen,
            #                    affine=affine, header=myheader)
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
    model_error, mean_error = LiFEvaluation(fdwi_data, trkstreamlines, gtab, subject=subject, header=header, roimask=roimask,
                                            affine=affine,display=display, outpathpickle=outpathpickle,
                                            outpathtrk=outpathtrk, processes=processes, outpathfig=outpathfig,
                                            verbose=verbose)
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
