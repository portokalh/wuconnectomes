#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eleftherios and Serge

Wenlin make some changes to track on the whole brain
Wenlin add for loop to run all the animals 2018-20-25
"""

from time import time
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.streamline import Streamlines
from dipy.tracking.local import LocalTracking
from dipy.reconst import peaks
from dipy.io.streamline import save_trk
import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes

from tracking_func import save_trk_heavy_duty
from dipy.io.utils import create_tractogram_header, get_reference_info

l=['N54717','N54718','N54719','N54720','N54722','N54759','N54760','N54761','N54762','N54763','N54764','N54765','N54766','N54770','N54771','N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806','N54807','N54818','N54824','N54825','N54826','N54837','N54838','N54843','N54844','N54856','N54857','N54858','N54859','N54860','N54861','N54873','N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892','N54893','N54897','N54898','N54899','N54900','N54915','N54916','N54917']
l=['CHASS']

# please set the parameter here

mypath = '/Users/alex/brain_data/E3E4/wenlin/'  # wenlin make this change
#mypath = ''

outpath = mypath + 'results/'  # wenlin make this change


#---------------------------------------------------------
tall = time()
for j in range(1):
    print(j+1)
    runno = l[j]
#    fdwi = mypath + 'N54900_nii4D_RAS.nii.gz'
#
#    ffalabels = mypath + 'fa_labels_warp_N54900_RAS.nii.gz'
#
#    fbvals = mypath + 'N54900_RAS_ecc_bvals.txt'
#
#    fbvecs = mypath + 'N54900_RAS_ecc_bvecs.txt'

#   wenlin make this change
    fdwi = '/Users/alex/brain_data/atlases/chass_symmetric3/bedpost_ESR120/data.nii.gz'
    
    ffalabels = '/Users/alex/brain_data/atlases/chass_symmetric3/bedpost_ESR120/nodif_brain_mask.nii'
    
    fbvecs = '/Users/alex/brain_data/atlases/chass_symmetric3/bedpost_ESR120/bvecs'
    
    fbvals = '/Users/alex/brain_data/atlases/chass_symmetric3/bedpost_ESR120/bvals'

    labels, affine_labels = load_nifti(ffalabels)

    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    # Correct flipping issue
    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[: ,2]]

    gtab = gradient_table(bvals, bvecs)

    data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)

    # Build Brain Mask
    bm = np.where(labels == 0, False, True)
    mask = bm

    sphere = get_sphere('repulsion724')


    from dipy.reconst.dti import TensorModel

    tensor_model = TensorModel(gtab)

    t1 = time()
    tensor_fit = tensor_model.fit(data, mask)
#    save_nifti('chassbmfa.nii.gz', tensor_fit.fa, affine)
#   wenlin make this change-adress name to each animal
    save_nifti(outpath+'bmfa'+runno+'.nii.gz', tensor_fit.fa, affine)
    fa = tensor_fit.fa
    duration1 = time() - t1
    #wenlin make this change-adress name to each animal
#    print('DTI duration %.3f' % (duration1,))
    print(runno+' DTI duration %.3f' % (duration1,))

    # Compute odfs in Brain Mask
    t2 = time()

    csa_model = CsaOdfModel(gtab, 6)

    csa_peaks = peaks_from_model(model=csa_model,
                                 data=data,
                                 sphere=peaks.default_sphere,  # issue with complete sphere
                                 mask=mask,
                                 relative_peak_threshold=.5,
                                 min_separation_angle=25,
                                 parallel=True,
                                 nbr_processes=4)

    duration2 = time() - t2
    print(duration2)\
    #wenlin make this change-adress name to each animal
#    print('CSA duration %.3f' % (duration2,))
    print(runno+' CSA duration %.3f' % (duration2,))


    t3 = time()

    # from dipy.tracking.local import ThresholdTissueClassifier
    from dipy.tracking.local import BinaryTissueClassifier # Wenlin Make this change

    # classifier = ThresholdTissueClassifier(csa_peaks.gfa, .25)
    classifier = BinaryTissueClassifier(bm) # Wenlin Make this change

    from dipy.tracking import utils

    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    #seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))

    seeds = utils.seeds_from_mask(mask, density=1,
                                  affine=np.eye(4))  # Wenlin make this change

    step_size=0.5
    stringstep=str(step_size); stringstep="_"+stringstep.replace(".","_")
    #stringstep=""
    streamlines_generator = LocalTracking(csa_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)

    # the function above will bring all streamlines in memory
    # streamlines = Streamlines(streamlines_generator)


    # save a smaller part by only keeping one in 10 streamlines
    sg_small = lambda: (s for i, s in enumerate(streamlines_generator) if i % 10 == 0)

    #wenlin make this change-adress name to each animal
    outpathfile=outpath+runno+"_bmCSA_detr_small"+stringstep+".trk"
    myheader=create_tractogram_header(outpathfile,*get_reference_info(fdwi))
    
    save_trk_heavy_duty(outpathfile, streamlines=sg_small,
                        affine=affine, header=myheader,
                        shape=mask.shape, vox_size=vox_size)


    # save everything - will generate a 20+ GBytes of data - hard to manipulate
    sg = lambda: (s for s in streamlines_generator)
    #wenlin make this change-adress name to each animal
#    save_trk_heavy_duty(outpath+"bmCSA_detr.trk", streamlines=sg, affine=affine,
#                        shape=mask.shape,
#                        vox_size=vox_size)
#    save_trk_heavy_duty(outpath+runno+"_bmCSA_detr.trk", streamlines=sg, affine=affine,
#                        shape=mask.shape,
#                        vox_size=vox_size)

    duration3 = time() - t3
    print(duration3)
    #wenlin make this change-adress name to each animal
#    print('Tracking duration %.3f' % (duration3,))
    print(runno+' Tracking duration %.3f' % (duration3,))

duration_all = time() -  tall
print('All animals tracking finished, running time is {}'.format(duration_all))
