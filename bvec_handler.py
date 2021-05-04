#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:12:00 2020

@author: alex
"""


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
from nibabel.tmpdirs import InTemporaryDirectory
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

from figures_handler import denoise_fig, show_bundles, window_show_test, LifEcreate_fig
from tract_eval import bundle_coherence, LiFEvaluation
from dif_to_trk import make_tensorfit, QCSA_tractmake

def fix_bvals_bvecs(fbvals, fbvecs, b0_threshold=50, atol=1e-2, outpath=None, identifier = "_fix", format="classic"):
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

    if outpath is None:
        base, ext = splitext(fbvals)
    else:
        base=str(outpath)+os.path.basename(fbvals).replace(".txt","")

    fbvals = base + identifier + ext
    if format == "classic":
        np.savetxt(fbvals, bvals)
    if format=="dsi_format":
        with open(fbvals, 'w') as File_object:
            for bval in bvals:
                if bval>10:
                    bval = int(round(bval))
                else:
                    bval=0
                File_object.write(str(bval) + "\t")
    base, ext = splitext(fbvecs)
    fbvecs = base + identifier + ext
    if format=="classic":
        if not os.path.isfile(fbvecs):
            np.savetxt(fbvecs, bvecs)
    #    with open(fbvecs, 'w') as f:
    #        f.write(str(bvec))
    if format=="dsi_format":
        with open(fbvecs, 'w') as File_object:
            for i in [0,1,2]:
                for j in np.arange(np.shape(bvecs)[0]):
                    if bvecs[j,i]==0:
                        bvec=0
                    else:
                        bvec=round(bvecs[j,i],3)
                    File_object.write(str(bvec)+"\t")
                File_object.write("\n")
            File_object.close()

    return fbvals, fbvecs

#def extract_bxh_bval_bvec(filepath):

def extractbvec_fromheader(source_file,fileoutpath=None,save=None,verbose=True):

    bvals = dsl = dpe = dro = None
    if save is not None:
        """
        if basepath==None:
            basepath=os.path.split(source_file)[0]+"/"
        if os.path.isdir(basepath):
            basepath=basepath+"/"
        else:
            if os.path.isdir(os.path.split(basepath)[0]):
                basepath=basepath
            else:
                basepath=os.path.join(sourcepath,basepath)
                basepath=basepath
        """
    filename = os.path.split(source_file)[1]
    if fileoutpath is None:
        basepath = os.path.split(source_file)[0]
        fileoutpath=basepath+filename.split('.')[0]+"_"
    with open(source_file, 'rb') as source:
        if verbose: print('INFO    : Extracting acquisition parameters')
        header_size=source.read(4)
        header_size=struct.unpack('I', header_size)
        if verbose: print('INFO    : Header size = ',int(header_size[0]))
        i=0
        stopsign=200
        bvec_start = False
        bval_start = False
        bvecs = []
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
            pattern5 = 'diffusiondirection'
            rx5 = re.compile(pattern5, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            pattern6 = 'bvalues'
            rx6 = re.compile(pattern6, re.IGNORECASE|re.MULTILINE|re.DOTALL)
            pattern7 = '<value>'
            i+=1
            if i==stopsign:
                print("hi")

            if bvec_start:
                rx7 = re.compile(pattern7, re.IGNORECASE | re.MULTILINE | re.DOTALL)
                bvec_start = False
                for a in rx7.findall(str(line)):
                    bvec_start = True
                    bvec_all = str(line).split('<value>')[1]
                    bvec = bvec_all.split('</value>')[0]
                    bvec = bvec.split(' ')
                    bvecs.append(bvec)

            for a in rx1.findall(str(line)):
                bvals_all = str(line).split(',')[1]
                bvals = bvals_all.split('\\')[0]
            for a in rx2.findall(str(line)):
                dsl_all = str(line).split(',')[1]
                dsl = dsl_all.split('\\')[0]
                #dsl=([float(s) for s in dsl_all.split() if s.isnumeric()])
            for a in rx3.findall(str(line)):
                dpe_all = str(line).split(',')[1]
                dpe = dpe_all.split('\\')[0]
            for a in rx4.findall(str(line)):
                dro_all = str(line).split(',')[1]
                dro = dro_all.split('\\')[0]
            for a in rx5.findall(str(line)):
                bvec_start = True
            for a in rx6.findall(str(line)):
                if not 'bvals_all' in locals():
                    bvals_all = str(line).split('"bvalues">')[1]
                    bvals_all = str(bvals_all).split('</datapoints>')[0]
                    bvals = bvals_all.split(' ')
                    bvals = "\n".join(bvals)

    if save == "all":
        bval_file=fileoutpath+"_bvals.txt"
        print(bval_file)
        File_object = open(bval_file,"w")
        File_object.write(bvals)
        File_object.close()

        if 'bvecs' in locals():
            bvecs_file=fileoutpath+"_bvec.txt"
            File_object = open(bvecs_file,"w")
            for bvec in bvecs:
                File_object.write(str(bvec[0]) + " " + str(bvec[1]) + " " + str(bvec[2]) + "\n")
            File_object.close()
            dsl = []
            dpe = []
            dro = []
            for bvec in bvecs:
                dsl.append(bvec[0])
                dpe.append(bvec[1])
                dro.append(bvec[2])
            dsl = "\n".join(dsl)
            dpe = "\n".join(dpe)
            dro = "\n".join(dro)

        elif dsl and dpe and dro:
            dsl_file=fileoutpath+"_dsl.txt"
            File_object = open(dsl_file,"w")
            File_object.write(dsl)
            File_object.close()

            dpe_file=fileoutpath+"_dpe.txt"
            File_object = open(dpe_file,"w")
            File_object.write(dpe)
            File_object.close()

            dro_file=fileoutpath+"_dro.txt"
            File_object = open(dro_file,"w")
            File_object.write(dro)
            File_object.close()
                         
            bvecs_file=fileoutpath+"_bvec.txt"
            File_object = open(bvecs_file,"w")
            dpe=dpe.split(" ")
            dsl=dsl.split(" ")
            dro=dro.split(" ")
            print(dpe)
            print(np.shape(dsl))
            for i in range(np.size(dsl)):
                File_object.write(str(dro[i]) + " " + str(dpe[i]) + " " + str(dsl[i]) + "\n")
            File_object.close()
        
        return bval_file, bvecs_file, bvals, dsl, dpe, dro

import shutil

def find_bval_bvecs(subjectpath):

    fbtable = glob.glob(os.path.join(subjectpath,"b_table.txt"))
    if np.size(fbtable)>0:
        bvals=[]
        bvecs=[]
        fbtable=fbtable[0]
        with open(fbtable, 'rb') as source:
            for line in source:
                bvals_all = str(line).split("'")[1]
                bvals_all = bvals_all.split("\\n")[0]
                bvals_all = bvals_all.split("\\t")
                bvals.append(float(bvals_all[0]))
                bvecs.append([float(bvals_all[1]), float(bvals_all[2]), float(bvals_all[3])])
        return np.array(bvals), np.array(bvecs)

    fbvals = glob.glob(os.path.join(subjectpath,"input_bvals.txt"))[0]
    fbvecs = glob.glob(os.path.join(subjectpath, "*input_gradient_matrix*"))[0]
    vals=[]
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
    bvals, bvecs = vals[0], vals[1]
    return bvals,bvecs

def writebfiles(subjectpath, subject, outpath = None, writeformat = "line", fbvals=None, fbvecs=None, overwrite=False):
    if outpath is None:
        outpath = subjectpath
    if fbvals is None or fbvecs is None:
        """
        bval_input = glob.glob(os.path.join(subjectpath, "input_bvals.txt"))
        if np.size(bval_input) > 0:
            with open(bval_input[0], 'rb') as source:
                pattern6 = 'b-value'
                rx1 = re.compile(pattern6, re.IGNORECASE | re.MULTILINE | re.DOTALL)
                for line in source:
                    for a in rx1.findall(str(line)):
                        bval_all = str(line).split('=')[1]
                        bval = bval_all.split('\\')[0]
                        bval = float(bval)
                        bval = int(round(bval))
                        pass
        """
        bvals, bvecs = find_bval_bvecs(subjectpath)

    #bvecs = glob.glob(os.path.join(subjectpath,"*input_gradient_matrix*"))

    bvec_file = os.path.join(outpath, subject+"_bvecs.txt")
    if os.path.exists(bvec_file) and overwrite:
        os.remove(bvec_file)
    if writeformat=="dsi_format":
        with open(bvec_file, 'w') as File_object:
            for i in [0,1,2]:
                for j in np.arange(np.shape(bvecs)[0]):
                    if bvecs[j,i]==0:
                        bvec=0
                    else:
                        bvec=round(bvecs[j,i],3)
                    File_object.write(str(bvec)+"\t")
                File_object.write("\n")
            File_object.close()
    else:
        shutil.copyfile(bvecs[0], bvec_file)
    """
        bval_input = glob.glob(os.path.join(subjectpath,"input_bvals.txt"))
        if np.size(bval_input)>0:
            with open(bval_input[0], 'rb') as source:
                pattern6 = 'b-value'
                rx1 = re.compile(pattern6, re.IGNORECASE | re.MULTILINE | re.DOTALL)
                for line in source:
                    for a in rx1.findall(str(line)):
                        bval_all = str(line).split('=')[1]
                        bval = bval_all.split('\\')[0]
                        bval = float(bval)
                        bval = int(round(bval))
                        pass
    bvals = []
    with open(bvec_file) as source:
        for line in source:
            if "0, 0, 0" in line:
                bvals.append(0)
            else:
                bvals.append(bval)
    """
    bval_file = os.path.join(outpath, subject+"_bvals.txt")
    if os.path.exists(bval_file) and overwrite:
        os.remove(bval_file)
    File_object = open(bval_file, "w")
    for bval in bvals:
        if writeformat == "line":
            File_object.write(str(bval) + "\n")
        if writeformat == "tab":
            File_object.write(str(bval) + "\t")
        if writeformat == "dsi_format":
            bval = int(round(bval))
            File_object.write(str(bval) + "\t")
    File_object.close()
    return bval_file, bvec_file

def extractbvals(dwipath, subject, outpath=None, writeformat="tab", overwrite=False):

    if os.path.isdir(dwipath):
        #subjectpath = os.path.join(dwipath, subject)
        subjectpath = glob.glob(os.path.join(os.path.join(dwipath, "*"+subject+"*")))
        subjectpath = subjectpath[0]
        if outpath is None:
            outpath = subjectpath
        fbvals = np.size(glob.glob(subjectpath + '*_bval*fix*'))
        fbvecs = np.size(glob.glob(subjectpath + '*_bvec*fix*'))
        if (fbvals == 0 and fbvecs == 0) or overwrite:
            #fbvals = (glob.glob(subjectpath + '*_bval*'))
            #fbvecs = (glob.glob(subjectpath + '*_bvec*'))
            fbvals=(glob.glob(os.path.join(subjectpath, '*' + subject + '*_bvals.txt')))
            fbvecs=(glob.glob(os.path.join(subjectpath, '*' + subject + '*_bvecs.txt')))
            if (np.size(fbvals) > 0 and np.size(fbvecs) > 0) and not overwrite:
                fbvals = fbvals[0]
                fbvecs = fbvecs[0]
                fix_bvals_bvecs(fbvals, fbvecs)
            else:
                bval_file, bvec_file = writebfiles(subjectpath, subject, outpath, writeformat=writeformat)
                fbvals, fbvecs = fix_bvals_bvecs(bval_file, bvec_file)
        return fbvals, fbvecs

    elif os.path.isfile(dwipath):
        dwifolder = os.path.dirname(os.path.abspath(dwipath))
        subjectpath = os.path.join(dwifolder, "*" + subject)
        if outpath is None:
            outpath = subjectpath
        fbvals = np.size(glob.glob(subjectpath + '*_bval*fix*'))
        fbvecs = np.size(glob.glob(subjectpath + '*_bvec*fix*'))
        if fbvals == 0 and fbvecs == 0:
            fbvals = np.size(glob.glob(subjectpath + '*_bval*'))
            fbvecs = np.size(glob.glob(subjectpath + '*_bvec*'))
            if (fbvals) == 0 or (fbvecs) == 0:
                bxhpath = dwipath.replace(".nii.gz", ".bxh")
                bxhpath = bxhpath.replace(".nii", ".bxh")
                fbvals, fbvecs, _, _, _, _ = extractbvec_fromheader(bxhpath,
                                                                    fileoutpath=os.path.join(dwifolder, subject),
                                                                    save="all")
                fix_bvals_bvecs(fbvals, fbvecs)

def rewrite_subject_bvalues(dwipath, subject, outpath=None, writeformat="tab", overwrite=False):

    if os.path.isdir(dwipath):
        #subjectpath = os.path.join(dwipath, subject)
        subjectpath = glob.glob(os.path.join(os.path.join(dwipath, "*"+subject+"*")))
        subjectpath = subjectpath[0]
        if outpath is None:
            outpath = subjectpath
        fbvals = np.size(glob.glob(subjectpath + '*_bval*fix*'))
        fbvecs = np.size(glob.glob(subjectpath + '*_bvec*fix*'))
        if (fbvals == 0 and fbvecs == 0) or overwrite:
            #fbvals = (glob.glob(subjectpath + '*_bval*'))
            #fbvecs = (glob.glob(subjectpath + '*_bvec*'))
            fbvals=(glob.glob(os.path.join(subjectpath, '*' + subject + '*_bvals.txt')))
            fbvecs=(glob.glob(os.path.join(subjectpath, '*' + subject + '*_bvecs.txt')))
            if (np.size(fbvals) > 0 and np.size(fbvecs) > 0) and not overwrite:
                fbvals = fbvals[0]
                fbvecs = fbvecs[0]
                #fix_bvals_bvecs(fbvals, fbvecs)
            else:
                bval_file, bvec_file = writebfiles(subjectpath, subject, outpath, writeformat=writeformat, overwrite=overwrite)
                #fbvals, fbvecs = fix_bvals_bvecs(bval_file, bvec_file)
        return fbvals, fbvecs

    elif os.path.isfile(dwipath):
        dwifolder = os.path.dirname(os.path.abspath(dwipath))
        subjectpath = os.path.join(dwifolder, "*" + subject)
        if outpath is None:
            outpath = subjectpath
        fbvals = np.size(glob.glob(subjectpath + '*_bval*fix*'))
        fbvecs = np.size(glob.glob(subjectpath + '*_bvec*fix*'))
        if fbvals == 0 and fbvecs == 0:
            fbvals = np.size(glob.glob(subjectpath + '*_bval*'))
            fbvecs = np.size(glob.glob(subjectpath + '*_bvec*'))
            if (fbvals) == 0 or (fbvecs) == 0:
                bxhpath = dwipath.replace(".nii.gz", ".bxh")
                bxhpath = bxhpath.replace(".nii", ".bxh")
                fbvals, fbvecs, _, _, _, _ = extractbvec_fromheader(bxhpath,
                                                                    fileoutpath=os.path.join(dwifolder, subject),
                                                                    save="all")
                fix_bvals_bvecs(fbvals, fbvecs)



def checkbxh(source_file, verbose = False):
    bxhtype = None
    with open(source_file, 'rb') as source:
        for line in source:

            pattern1 = "psdname"
            rx1 = re.compile(pattern1, re.IGNORECASE | re.MULTILINE | re.DOTALL)

            for a in rx1.findall(str(line)):
                sequencetype = str(line).split('<psdname>')[1]
                sequencetype = str(sequencetype).split('</psdname>')[0]
                break
    if sequencetype == "epi2muse":
        bxhtype = "dwi"
    elif sequencetype == "BRAVO":
        bxhtype = "T1"

    return bxhtype

def same_axis(first,second):
    if first==second:
        return 1
    if first=="A" and second=="P" or first=="P" and second=="A" or first==" ":
        print('hi')

"""
def bvec_reorient(bvecs,orig_orient="ARI",new_orient="RAS"):

    for i in range(3):
        #if new_orient[i]
        echo('jo')

    new_bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]]
"""