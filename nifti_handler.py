
import nibabel as nib
import numpy as np
import os, glob
from dipy.io.gradients import read_bvals_bvecs
from bvec_handler import fix_bvals_bvecs#, extractbvec_fromheader
import pathlib
from BIAC_tools import send_mail
from dipy.core.gradients import gradient_table


def getfa(mypath, subject, bvec_orient=[1, 2, 3], verbose=None):

    # fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    fapath = mypath + '/' + subject + '_fa_RAS.nii.gz'
    if os.path.exists(fapath):
        fapath = mypath + '/' + subject + '_fa_RAS.nii.gz'
    # fdwi_data, affine, vox_size = load_nifti(fdwipath, return_voxsize=True)

    if os.path.exists(mypath + '/' + subject + '_fa_RAS.nii.gz'):
        fapath = (mypath + '/' + subject + '_fa_RAS.nii.gz')
    elif os.path.exists(mypath+'/'+'bmfa' + subject+'_wholebrain_.nii.gz'):
        fapath = (mypath+'/'+'bmfa' + subject+'_wholebrain_.nii.gz')
    else:
        print("Could not find at either "+ (mypath + '/' + subject + '_fa_RAS.nii.gz') + " or " + (mypath+'/'+'bmfa' + subject+'_wholebrain.nii.gz'))

    if verbose:
        txt = "Extracting information from the fa file located at " + fapath
        print(txt)
        send_mail(txt, subject="Begin data extraction")

    if 'fapath' not in locals():
        txt = "The fa of subject " + subject + " was not detected at " + fapath + ", exit"
        print(txt)
        send_mail(txt, subject="Error")
        return (0, 0, 0, 0, 0, 0, 0, 0)

    img = nib.load(fapath)
    fa_data = img.get_data()
    vox_size = img.header.get_zooms()[:3]
    affine = img.affine
    hdr = img.header
    header = get_reference_info(fapath)
    del (img)

    try:
        fbvals = glob.glob(mypath + '/' + subject + '*_bvals_fix.txt')[0]
        fbvecs = glob.glob(mypath + '/' + subject + '*_bvec_fix.txt')[0]
    except IndexError:
        fbvals = glob.glob(mypath + '/' + subject + '*_bvals.txt')[0]
        fbvecs = glob.glob(mypath + '/' + subject + '*_bvec.txt')[0]
        fbvals, fbvecs = fix_bvals_bvecs(fbvals, fbvecs)
    print(fbvecs)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    # bvecs = np.c_[bvecs[:, 0], -bvecs[:, 1], bvecs[:, 2]]  # FOR RAS according to Alex
    # bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]] #FOR RAS

    # bvecs = np.c_[bvecs[:, -], bvecs[:, 0], -bvecs[:, 2]] #estimated for RAS based on headfile info
    bvec_sign = bvec_orient / np.abs(bvec_orient)
    bvecs = np.c_[bvec_sign[0] * bvecs[:, np.abs(bvec_orient[0]) - 1], bvec_sign[1] * bvecs[:, np.abs(bvec_orient[1]) - 1],
        bvec_sign[2] * bvecs[:, np.abs(bvec_orient[2]) - 1]]

    # bvecs = np.c_[bvecs[:, 1], bvecs[:, 0], -bvecs[:, 2]]
    # bvecs = np.c_[-bvecs[:, 1], bvecs[:, 0], bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    # bm = np.where(labels == 0, False, True)
    # mask = bm

    return fa_data, affine, gtab, vox_size, hdr, header



def get_reference_info(reference, affine = np.eye(4).astype(np.float32)):
    """ Will compare the spatial attribute of 2 references

    Parameters
    ----------
    reference : Nifti or Trk filename, Nifti1Image or TrkFile, Nifti1Header or
        trk.header (dict)
        Reference that provides the spatial attribute.

    Returns
    -------
    output : tuple
        - affine ndarray (4,4), np.float32, tranformation of VOX to RASMM
        - dimensions list (3), int, volume shape for each axis
        - voxel_sizes  list (3), float, size of voxel for each axis
        - voxel_order, string, Typically 'RAS' or 'LPS'
    """

    is_nifti = False
    is_trk = False
    is_sft = False
    if isinstance(reference, str):
        try:
            header = nib.load(reference).header
            is_nifti = True
        except nib.filebasedimages.ImageFileError:
            pass
        try:
            header = nib.streamlines.load(reference, lazy_load=True).header
            _, extension = os.path.splitext(reference)
            if extension == '.trk':
                is_trk = True
        except ValueError:
            pass
    elif isinstance(reference, nib.nifti1.Nifti1Image):
        header = reference.header
        is_nifti = True
    elif isinstance(reference, nib.streamlines.trk.TrkFile):
        header = reference.header
        is_trk = True
    elif isinstance(reference, nib.nifti1.Nifti1Header):
        header = reference
        is_nifti = True
    elif isinstance(reference, dict) and 'magic_number' in reference:
        header = reference
        is_trk = True
    elif isinstance(reference, dipy.io.stateful_tractogram.StatefulTractogram):
        is_sft = True

    if is_nifti:
        if np.sum(header['srow_x']) != 0:
            affine[0, 0:4] = header['srow_x']
            affine[1, 0:4] = header['srow_y']
            affine[2, 0:4] = header['srow_z']
        dimensions = header['dim'][1:4]
        voxel_sizes = header['pixdim'][1:4]
        voxel_order = ''.join(nib.aff2axcodes(affine))
    elif is_trk:
        affine = header['voxel_to_rasmm']
        dimensions = header['dimensions']
        voxel_sizes = header['voxel_sizes']
        voxel_order = header['voxel_order']
    elif is_sft:
        affine, dimensions, voxel_sizes, voxel_order = reference.space_attribute
    else:
        raise TypeError('Input reference is not one of the supported format')

    if isinstance(voxel_order, np.bytes_):
        voxel_order = voxel_order.decode('utf-8')

    return affine, dimensions, voxel_sizes, voxel_order


def getdwidata(mypath, subject, bvec_orient=[1,2,3], verbose=None):

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
    elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '/'):
        fdwipath = glob.glob(mypath + '/' + subject + '/' + subject + '*nii4D*.nii*')[0]

    #fdwi_data, affine, vox_size = load_nifti(fdwipath, return_voxsize=True)

    if verbose:
        txt = "Extracting information from the dwi file located at " + fdwipath
        print(txt)
        send_mail(txt,subject="Begin data extraction")

    if 'fdwipath' not in locals():
        txt = "The subject " + subject + " was not detected, exit"
        print(txt)
        send_mail(txt,subject="Error")
        return(0,0,0,0,0,0,0,0)

    mypath = str(pathlib.Path(fdwipath).parent.absolute())

    try:
        fbvals = glob.glob(mypath + '/' + subject + '*_bvals_fix.txt')[0]
        fbvecs = glob.glob(mypath + '/' + subject + '*_bvec_fix.txt')[0]
    except IndexError:
        fbvals = glob.glob(mypath + '/' + subject + '*_bvals.txt')[0]
        fbvecs = glob.glob(mypath + '/' + subject + '*_bvec*.txt')[0]
        fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
    print(fbvecs)
    bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    #bvecs = np.c_[bvecs[:, 0], -bvecs[:, 1], bvecs[:, 2]]  # FOR RAS according to Alex
    #bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]] #FOR RAS

    #bvecs = np.c_[bvecs[:, -], bvecs[:, 0], -bvecs[:, 2]] #estimated for RAS based on headfile info
    bvec_sign = bvec_orient/np.abs(bvec_orient)
    bvecs = np.c_[bvec_sign[0]*bvecs[:, np.abs(bvec_orient[0])-1], bvec_sign[1]*bvecs[:, np.abs(bvec_orient[1])-1],
                  bvec_sign[2]*bvecs[:, np.abs(bvec_orient[2])-1]]


    img = nib.load(fdwipath)
    fdwi_data = img.get_data()
    vox_size = img.header.get_zooms()[:3]
    affine = img.affine
    hdr = img.header
    del(img)

    if os.path.exists(mypath+'/Reg_'+subject+'_nii4D_brain_mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath+'/Reg_'+subject+'_nii4D_brain_mask.nii.gz')
    elif os.path.exists(mypath+'/'+subject+'_chass_symmetric3_labels_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath+'/'+subject+'_chass_symmetric3_labels_RAS.nii.gz')
    elif os.path.exists(mypath+'/'+subject+'_chass_symmetric3_labels_RAS_combined.nii.gz'):
        labels, affine_labels = load_nifti(mypath+'/'+subject+'_chass_symmetric3_labels_RAS_combined.nii.gz')
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject +'_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject +'_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii.gz')
    elif os.path.exists(mypath + '/mask.nii'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii')
    else:
        print('mask not found, taking all non null values in nii file instead (not recommended for complex operations)')
        labels = np.ones(fdwi_data.shape[0:3])
        affine_labels = affine

    #bvecs = np.c_[bvecs[:, 1], bvecs[:, 0], -bvecs[:, 2]]
    #bvecs = np.c_[-bvecs[:, 1], bvecs[:, 0], bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    #bm = np.where(labels == 0, False, True)
    #mask = bm

    header = get_reference_info(fdwipath)
    return fdwi_data, affine, gtab, labels, vox_size, fdwipath, hdr, header
