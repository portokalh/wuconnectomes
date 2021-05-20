

from time import time
import numpy as np
import os

import multiprocessing as mp
import pickle

from tract_manager import getdwidata
from bvec_handler import extractbvec_fromheader
from BIAC_tools import send_mail
from Daemonprocess import MyPool

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from BIAC_tools import isempty
from dipy.io.streamline import load_trk
import tract_save
from dipy.io.image import load_nifti, save_nifti
from dipy.io.utils import create_tractogram_header

from dif_to_trk import save_roisubset


l = ['N57433', 'N57434', 'N57435', 'N57436', 'N57437', 'N57440']

BIGGUS_DISKUS = "/Volumes/Badea/Lab/mouse"

dwipath = BIGGUS_DISKUS + "/C57_JS/DWI_RAS_40subj/"
outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS_40subj/"
outtrkpath = BIGGUS_DISKUS + "/C57_JS/TRK_RAS/"

subject = l[0]

trkfile = "/Volumes/Badea/Lab/mouse/C57_JS/TRK_RAS/N57433_smaller_stepsize_2.trk"

trkdata = load_trk(trkfile, 'same')
trkdata.to_vox()
if hasattr(trkdata, 'space_attribute'):
    header = trkdata.space_attribute
elif hasattr(trkdata, 'space_attributes'):
    header = trkdata.space_attributes
trkstreamlines = trkdata.streamlines


rois = ["wholebrain"]
roislist = [['hypothalamus', 'septum'], ['fimbria'], ['corpus_callosum'], ['primary_motor_cortex']]

roisexcel = BIGGUS_DISKUS + "/atlases/CHASSSYMM3AtlasLegends.xlsx"

labelfile = "/Volumes/Badea/Lab/mouse/C57_JS/DWI_RAS/N57433_chass_symmetric2_labels_RAS.nii.gz"
labelmask, affine_labels = load_nifti(labelfile)

stringstep="2"
stepsize = 2

#ratios = 100
trkpath = outtrkpath


"""
ratio=10
ministream = []
trkminipath = "/Volumes/Badea/Lab/mouse/C57_JS/TRK_RAS/N57433_smaller_stepsize_2.trk"
for idx, stream in enumerate(trkstreamlines):
    if (idx % ratio) == 0:
        ministream.append(stream)
myheader = create_tractogram_header(trkminipath, *header)
ratioed_sl_gen = lambda: (s for s in ministream)
tract_save.save_trk_heavy_duty(trkminipath, streamlines=ratioed_sl_gen,
                                   affine=affine, header=myheader)
"""
strproperty=""

affine=affine_labels
save_roisubset(trkfile, roislist, roisexcel, labelmask)
