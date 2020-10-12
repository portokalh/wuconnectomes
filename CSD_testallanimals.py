#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eleftherios and Serge

Wenlin make some changes to track on the whole brain
Wenlin add for loop to run all the animals 2018-20-25
"""
#do save response in order to reload and reuse at a later time
# for i in `pgrep -f python*`; do kill $i; done

from time import time
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.tracking.streamline import Streamlines
from dipy.reconst import peaks
from dipy.io.streamline import save_trk
import nibabel as nib
from skimage.transform import rescale
from scipy.ndimage import zoom
from dipy.reconst.csdeconv import (ConstrainedSphericalDeconvModel,
                                   auto_response)
from dipy.tracking import utils
from dipy.segment.mask import applymask
from tracking_func import save_trk_heavy_duty
from dipy.io.utils import create_tractogram_header, get_reference_info


import multiprocessing
from dipy.tracking.local import LocalTracking, BinaryTissueClassifier


radius=10 #20

l=['N54900', 'N54717','N54718','N54719','N54720','N54722','N54759','N54760','N54761','N54762','N54763','N54764','N54765','N54766','N54770','N54771','N54772','N54798','N54801','N54802','N54803','N54804','N54805','N54806','N54807','N54818','N54824','N54825','N54826','N54837','N54838','N54843','N54844','N54856','N54857','N54858','N54859','N54860','N54861','N54873','N54874','N54875','N54876','N54877','N54879','N54880','N54891','N54892','N54893','N54897','N54898','N54899','N54900','N54915','N54916','N54917']




# please set the parameter here

#mypath = '/Users/alex/brain_data/E3E4/wenlin/'  # wenlin make this change
mypath = '/Users/alex/dipy_workshop/DIPYWORKSHOP_ABDATA/'

outpath = mypath + 'results/'  # wenlin make this change


#---------------------------------------------------------
tall = time()
#for j in range(55):

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
    fdwi = mypath+runno+'_nii4D_RAS.nii.gz'
    
    ffalabels = mypath+'fa_labels_warp_'+runno+'_RAS.nii.gz'
    
    fbvals = mypath+runno+'_RAS_ecc_bvals.txt'
    
    fbvecs = mypath+runno+'_RAS_ecc_bvecs.txt'

    labels, affine_labels = load_nifti(ffalabels)

    print('affine for labels', affine_labels)

    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    # Correct flipping issue
    bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[: ,2]]

    gtab = gradient_table(bvals, bvecs)

    data, affine, vox_size = load_nifti(fdwi, return_voxsize=True)
    scale= [0.5, 0.5, 0.5, 1.]
    data= zoom(data, zoom= scale, order= 1, mode= 'constant')
    labels= zoom(labels, zoom= scale[:3], order= 1, mode= 'constant')
    labels= labels.astype(int)


    # Build Brain Mask
    #bm = np.where(labels == 0, False, True)
    bm = (labels!=0)*1
    mask=bm
    #sphere = get_sphere('repulsion724')
    sphere = get_sphere('repulsion200')


    from dipy.reconst.dti import TensorModel
    tensor_model = TensorModel(gtab)

    t1 = time()
    tensor_fit = tensor_model.fit(data, bm)
#    save_nifti('bmfa.nii.gz', tensor_fit.fa, affine)
#   wenlin make this change-adress name to each animal

    affine= affine @ np.diag(scale)

    save_nifti(outpath+'bmfa'+runno+'.nii.gz', tensor_fit.fa, affine)
    fa = tensor_fit.fa
    duration1 = time() - t1
    #wenlin make this change-adress name to each animal
#    print('DTI duration %.3f' % (duration1,))
    print(runno+' DTI duration %.3f' % (duration1,))
    #response : 3.96154132e-04, 9.23377324e-05, 9.23377324e-05

#replace CSA with CSD
# Build Brain Mask
    t2=time()
    myroi=120
    # bm = np.where(labels == myroi, False, True)
    roimask = (labels == myroi) * 1
    # Compute odfs in Brain Mask
    t2 = time()
    #add fa_thresh=0.5 for wm
    response, ratio, nvl = auto_response(gtab, applymask(data,roimask), fa_thresh=0.5, roi_radius=radius, return_number_of_voxels=True)
    print('We use the roi_radius={},\nand the response is {},\nthe ratio is {},\nusing {} of voxels'
    .format(radius, response, ratio, nvl))
    #np.savetxt(outpath+runno+'.txt',response)
    csd_model = ConstrainedSphericalDeconvModel(gtab, response, sh_order=4)

    #csd peack is 4D with very voxel having a set of sh coefficients
    '''
    csd_peaks = peaks_from_model(csd_model, data, sphere=peaks.default_sphere,
                                 relative_peak_threshold=0.5, min_separation_angle=60,mask=bm, 
                                 return_odf=True,return_sh=True, parallel=True, sh_order=4, sh_basis_type='tournier07', 
                                 nbr_processes=multiprocessing.cpu_count())
                                 '''
    #npeaks?
    #bm4=np.repeat(bm[:, :, :, np.newaxis], shape.data[4], axis=4)
    #data=np.matmult(data,rbm4)
    csd_peaks = peaks_from_model(csd_model, applymask(data,bm), sphere=peaks.default_sphere,
                                 relative_peak_threshold=0.5, min_separation_angle=25, mask=bm,
                                 return_odf=False,return_sh=True, parallel=True, sh_order=4, sh_basis_type='tournier07', 
                                 nbr_processes=3)
                                 #nbr_processes=multiprocessing.cpu_count())
    
    #nib.save(nib.Nifti1Image(csd_peaks.shm_coeff.astype(np.float32), affine), outpath+runno+'fod.nii.gz')
    nib.save(nib.Nifti1Image(csd_peaks.shm_coeff.astype(np.float32), np.eye(4)), outpath+runno+'fod.nii.gz')
    
    duration2 = time() - t2
    print(runno+' CSD duration %.3f' % (duration2,))
    #exit()

    
    #mask and classifier
    seeds = utils.seeds_from_mask(bm, density=1,affine=np.eye(4))  
    classifier = BinaryTissueClassifier(bm)
    #local tracking
    streamlines_generator = LocalTracking(csd_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=.5)

    #end replace CSA with CSD

    """
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
        #end CSA
    """
    t3 = time()

    # from dipy.tracking.local import ThresholdTissueClassifier

    # classifier = ThresholdTissueClassifier(csa_peaks.gfa, .25)
    classifier = BinaryTissueClassifier(bm) # Wenlin Make this change


    # generates about 2 seeds per voxel
    # seeds = utils.random_seeds_from_mask(fa > .2, seeds_count=2,
    #                                      affine=np.eye(4))

    # generates about 2 million streamlines
    #seeds = utils.seeds_from_mask(fa > .2, density=1,
    #                              affine=np.eye(4))

    seeds = utils.seeds_from_mask(bm, density=1,
                                  affine=np.eye(4))  # Wenlin make this change

    step_size=0.5
    stringstep=str(step_size); stringstep="_"+stringstep.replace(".","_")
    #stringstep=""
    streamlines_generator = LocalTracking(csd_peaks, classifier,
                                          seeds, affine=np.eye(4), step_size=step_size)

    # the function above will bring all streamlines in memory
    # streamlines = Streamlines(streamlines_generator)


    # save a smaller part by only keeping one in 10 streamlines
    sg_small = lambda: (s for i, s in enumerate(streamlines_generator) if i % 10 == 0)

    #wenlin make this change-adress name to each animal
#    save_trk_heavy_duty(outpath+"bmCSA_detr_small.trk", streamlines=sg_small,
#                        affine=affine,
#                        shape=mask.shape, vox_size=vox_size)
    outpathfile=outpath+runno+"_bmCSA_detr_small_header"+stringstep+".trk"
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
