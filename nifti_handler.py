
import nibabel as nib
import numpy as np
import os, glob
from dipy.io.gradients import read_bvals_bvecs
from bvec_handler import fix_bvals_bvecs, checkbxh, reorient_bvecs
import pathlib
from BIAC_tools import send_mail
from dipy.core.gradients import gradient_table
from diff_preprocessing import make_tensorfit
from dipy.io.image import load_nifti
import shutil
from convert_atlas_mask import convert_labelmask, atlas_converter
import errno
from computer_nav import load_nifti_remote, glob_remote, checkfile_exists_remote, read_bvals_bvecs_remote

def getfa(mypath, subject, bvec_orient, verbose=None):

    # fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    fapath = mypath + '/' + subject + '_fa_RAS.nii.gz'
    if os.path.exists(fapath):
        fapath = mypath + '/' + subject + '_fa_RAS.nii.gz'
    # fdwi_data, affine, vox_size = load_nifti(fdwipath, return_voxsize=True)

    if os.path.exists(mypath + '/' + subject + '_fa_RAS.nii.gz'):
        fapath = (mypath + '/' + subject + '_fa_RAS.nii.gz')
    elif os.path.exists(mypath+'/'+'bmfa' + subject+'_wholebrain_.nii.gz'):
        fapath = (mypath+'/'+'bmfa' + subject+'_wholebrain_.nii.gz')
    elif os.path.exists(mypath + '/' + subject + '/' + 'bmfa' + subject + '.nii.gz'):
        fapath = (mypath + '/' + subject + '/' + 'bmfa' + subject + '.nii.gz')
    else:
        print("Could not find the fa file anywhere")
        print("Will attempt to create new fa file")
        fdiff_data, affine, gtab, mask, vox_size, fdiffpath, hdr, header = getdiffdata(mypath, subject, bvec_orient)
        fapath = make_tensorfit(fdiff_data, mask, gtab, affine, subject, outpath=os.path.dirname(fdiffpath), strproperty="", verbose=verbose)
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

    """
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
    """
    return fa_data, affine, vox_size, hdr, header

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

def getrefpath(mypath, subject, reference = 'fa', verbose=None, sftp=None):

    if sftp is None:
        if os.path.exists(os.path.join(mypath,subject+"_subjspace_"+reference+".nii.gz")):
            refpath = (os.path.join(mypath,subject+"_subjspace_"+reference+".nii.gz"))
        elif os.path.exists(os.path.join(mypath,subject+"_"+reference+"_RAS.nii.gz")):
            refpath = (os.path.join(mypath,subject+"_coreg_RAS.nii.gz"))

        if 'refpath' not in locals():
            txt = "The subject " + subject + " was not detected, exit"
            print(txt)
            send_mail(txt, subject="Error")
            return None
    else:
        refpaths = glob_remote(os.path.join(mypath,subject+"_subjspace_"+reference+"*.nii.gz"),sftp)
        if np.size(refpaths)==1:
            refpath = refpaths[0]
        else:
            for refpath_t in refpaths:
                if 'RAS' not in refpath_t:
                    refpath = refpath_t
    return refpath



def getdiffpath_old(mypath, subject, denoise="", verbose=None):

    if denoise is None:
        denoise=""

    subjfolder = glob.glob(os.path.join(mypath, "*" + subject + "*/"))
    if np.size(subjfolder)==1:
        subjfolder = subjfolder[0]
    else:
        subjfolder = None
    print('hi')
    print(os.path.join(mypath,subject+"_subjspace_coreg.nii.gz"))
    if os.path.isfile(mypath) and os.path.exists(mypath):
        fdiffpath = mypath
    #elif os.path.exists(os.path.join(mypath,subject+"_"+denoise+"_diff.nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_"+denoise+"_diff.nii.gz"))
    #elif os.path.exists(os.path.join(mypath,subject+"_"+denoise+".nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_"+denoise+".nii.gz"))
    #elif os.path.exists(os.path.join(mypath,subject+"_rawnii.nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_rawnii.nii.gz"))
    elif os.path.exists(os.path.join(mypath,subject+"_subjspace_coreg.nii.gz")):
        fdiffpath = (os.path.join(mypath,subject+"_subjspace_coreg.nii.gz"))
    #elif os.path.exists(os.path.join(mypath,subject+"_coreg.nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_coreg.nii.gz"))
    #elif os.path.exists(os.path.join(mypath,subject+"_coreg.nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_coreg.nii.gz"))
    elif os.path.exists(os.path.join(mypath,subject+"_coreg_RAS.nii.gz")):
        fdiffpath = (os.path.join(mypath,subject+"_coreg_RAS.nii.gz"))
    #elif np.size(glob.glob(os.path.join(mypath,subject+"*_dwi.nii.gz"))) == 1:
    #    fdiffpath = glob.glob(os.path.join(mypath,subject+"*_dwi.nii.gz"))[0]
    #elif os.path.exists(mypath + '/Reg_' + subject + '_nii4D.nii.gz'):
    #    fdiffpath = mypath + '/Reg_' + subject + '_nii4D.nii.gz'
    #elif os.path.exists(mypath + '/nii4D_' + subject + '.nii'):
    #    fdiffpath = mypath + '/nii4D_' + subject + '.nii'
    #elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
    #    fdiffpath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    #elif os.path.exists(mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'):
    #    fdiffpath = mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'
    #elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
    #    fdiffpath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    #elif os.path.exists(mypath + '/' + subject + '/') and np.size(glob.glob(os.path.join(subjfolder, subject + '*nii4D*.nii*'))) > 0:
    #    fdiffpath = glob.glob(os.path.join(subjfolder, subject + '*nii4D*.nii*'))[0]
    #elif os.path.exists(os.path.join(mypath,subject+"_dwi.nii.gz")):
    #    fdiffpath = (os.path.join(mypath,subject+"_dwi.nii.gz"))
    elif os.path.exists(mypath) and subjfolder is not None and np.size(glob.glob(os.path.join(subjfolder, "*.bxh"))) > 0:
        subjbxh = glob.glob(os.path.join(subjfolder, "*.bxh"))
        for bxhfile in subjbxh:
            bxhtype = checkbxh(bxhfile, False)
            if bxhtype == "diff":
                fdiffpath = bxhfile.replace(".bxh", ".nii.gz")
                break

    if 'fdiffpath' not in locals():
        txt = "The subject " + subject + " was not detected, exit"
        print(txt)
        send_mail(txt, subject="Error")
        return None

    return(fdiffpath)

def getdiffpath(mypath, subject, denoise="", verbose=None, sftp=None):
    if denoise is None:
        denoise=""

    listoptions = [os.path.join(mypath,subject+"_subjspace_coreg.nii.gz"), os.path.join(mypath,subject+"_subjspace_coreg_RAS.nii.gz"), os.path.join(mypath,subject+"_coreg_RAS.nii.gz"), os.path.join(mypath,subject+"_coreg_diff.nii.gz")]

    if sftp is None:
        for list_option in listoptions:
            if '*' in list_option:
                option = glob.glob(list_option)
                if np.size(option) > 0:
                    fdiffpath = option[0]
                    break
            else:
                if os.path.exists(list_option):
                    fdiffpath = list_option
    else:
        for list_option in listoptions:
            option = glob_remote(list_option, sftp)
            if np.size(option) > 0:
                if np.size(option)>1:
                    raise Warning("too many diffusion fitting parameters!!")
                fdiffpath = option[0]
                break

    if 'fdiffpath' not in locals():
        txt = "The subject " + subject + " was not detected, exit"
        print(txt)
        send_mail(txt, subject="Error")
        return None

    return fdiffpath

def extract_nii_info(path, verbose=None, sftp=None):
    if verbose:
        txt = "Extracting information from the nifti file located at " + path
        print(txt)
        send_mail(txt, subject="Begin data extraction")
    if sftp is None:
        img = nib.load(path)
        data = img.get_data()
        vox_size = img.header.get_zooms()[:3]
        affine = img.affine
        header = img.header
        ref_info = get_reference_info(path)
    else:
        data, affine, vox_size, header, ref_info = load_nifti_remote(path, sftp)
    return data, affine, vox_size, header, ref_info

def getrefdata(mypath, subject, reference, verbose=None, sftp=None):

    ref_fpath = getrefpath(mypath, subject, reference, verbose=verbose, sftp=sftp)
    if ref_fpath is None:
        return None, None, None, None, None, None
    ref_data, affine, vox_size, header, ref_info = extract_nii_info(ref_fpath, verbose,sftp=sftp)

    return ref_data, affine, vox_size, ref_fpath, header, ref_info


def getdiffdata(mypath, subject, denoise="", verbose=None,sftp=None):

    diff_fpath = getdiffpath(mypath, subject, denoise=denoise, verbose=verbose,sftp=sftp)
    diff_data, affine, vox_size, header, ref_info = extract_nii_info(diff_fpath, verbose,sftp=sftp)

    return diff_data, affine, vox_size, diff_fpath, header, ref_info

def get_bvals_bvecs(mypath, subject,sftp=None):
    if sftp is None:
        try:
            fbvals = glob.glob(mypath + '/' + subject + '*_bvals_fix.txt')[0]
            fbvecs = glob.glob(mypath + '/' + subject + '*_bvec_fix.txt')[0]
        except IndexError:
            print(mypath + '/' + subject + '*_bvals.txt')
            fbvals = glob.glob(mypath + '/' + subject + '*_bvals.txt')[0]
            fbvecs = glob.glob(mypath + '/' + subject + '*_bvec*.txt')[0]
            fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
        print(fbvecs)
        bvals, bvecs = read_bvals_bvecs(fbvals, fbvecs)

    else:
        try:
            fbvals = glob_remote(mypath + '/' + subject + '*_bvals_fix.txt', sftp)[0]
            fbvecs = glob_remote(mypath + '/' + subject + '*_bvecs_fix.txt', sftp)[0]
        except IndexError:
            print(mypath + '/' + subject + '*_bvals.txt')
            fbvals = glob_remote(mypath + '/' + subject + '*_bvals.txt', sftp)[0]
            fbvecs = glob_remote(mypath + '/' + subject + '*_bvec*.txt', sftp)[0]
            fbvals, fbvecs = fix_bvals_bvecs(fbvals, fbvecs,sftp=sftp)
        print(fbvecs)
        bvals, bvecs = read_bvals_bvecs_remote(fbvals, fbvecs,sftp=sftp)
    return bvals, bvecs


def getgtab(mypath, subject, bvec_orient=[1,2,3],sftp=None):

    bvals, bvecs = get_bvals_bvecs(mypath, subject,sftp=sftp)
    bvecs = reorient_bvecs(bvecs, bvec_orient)
    #bvec_sign = bvec_orient/np.abs(bvec_orient)
    #bvecs = np.c_[bvec_sign[0]*bvecs[:, np.abs(bvec_orient[0])-1], bvec_sign[1]*bvecs[:, np.abs(bvec_orient[1])-1],
    #              bvec_sign[2]*bvecs[:, np.abs(bvec_orient[2])-1]]

    gtab = gradient_table(bvals, bvecs)

    return gtab

def getb0s(mypath, subject,sftp):
    bvals, _ = get_bvals_bvecs(mypath, subject,sftp)
    b0s = []
    i=0
    for bval in bvals:
        if bval < 10:
            b0s.append(i)
        i += 1
    return(b0s)

def getdiffdata_all(mypath, subject, bvec_orient=[1,2,3], denoise="", verbose=None,sftp=None):

    fdiff_data, affine, vox_size, fdiffpath, header, ref_info = getdiffdata(mypath, subject, denoise=denoise, verbose=verbose,sftp=sftp)
    mypath = str(pathlib.Path(fdiffpath).parent.absolute())

    if bvec_orient is None:
        img = nib.load(fdiffpath)
        fdiff_data = img.get_data()
        vox_size = img.header.get_zooms()[:3]
        affine = img.affine
        hdr = img.header
        header = get_reference_info(fdiffpath)
        gtab = None
        return fdiff_data, affine, gtab, vox_size, fdiffpath, hdr, header

    gtab = getgtab(mypath, subject, bvec_orient,sftp=sftp)

    #bvecs = np.c_[bvecs[:, 0], -bvecs[:, 1], bvecs[:, 2]]  # FOR RAS according to Alex
    #bvecs = np.c_[bvecs[:, 0], bvecs[:, 1], -bvecs[:, 2]] #FOR RAS

    #bvecs = np.c_[bvecs[:, -], bvecs[:, 0], -bvecs[:, 2]] #estimated for RAS based on headfile info

    return fdiff_data, affine, gtab, vox_size, fdiffpath, header, ref_info

def getlabelmask(mypath, subject, verbose=None, sftp=None):

    list_options = [mypath + '/' + subject + '/' + subject + '*labels.nii.gz',
                    mypath + '/*' + subject + '*labels.nii.gz', mypath + '/' + subject + '_labels_RAS.nii.gz', (mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'),
                    (mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'), (mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'),
                    (mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'), (mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'), (mypath + '/mask.nii.gz'),
                    (mypath + '/mask.nii')]

    if sftp is None:
        for list_option in list_options:
            if '*' in list_option:
                labelsoption = glob.glob(list_option)
                if np.size(labelsoption) > 0:
                    labelspath = labelsoption[0]
                    break
            else:
                if os.path.exists(list_option):
                    labelspath=list_option
    else:
        for list_option in list_options:
            labelsoption = glob_remote(list_option, sftp)
            if np.size(labelsoption) > 0:
                labelspath = labelsoption[0]
                break

    """       
    labelsoption = glob.glob(mypath + '/' + subject + '/' + subject + '*labels.nii.gz')
    print(mypath + '/' + subject + '/' + subject + '*labels.nii.gz')
    if np.size(labelsoption)>0:
        labelspath = labelsoption[0]
    labelsoption = glob.glob(mypath + '/*' + subject + '*labels.nii.gz')
    print((mypath + '/' + subject + '_labels_RAS.nii.gz'))
    if np.size(labelsoption)>0:
        labelspath = labelsoption[0]
    elif os.path.exists(mypath + '/' + subject + '_labels_RAS.nii.gz'):
        labelspath = mypath + '/' + subject + '_labels_RAS.nii.gz'
    elif os.path.exists(mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'):
        labelspath = mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'):
        labelspath = mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'):
        labelspath = mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labelspath = mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labelspath = mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labelspath = mypath + '/mask.nii.gz'
    elif os.path.exists(mypath + '/mask.nii'):
        labelspath = mypath + '/mask.nii'
    """

    if 'labelspath' in locals():
        if sftp is None:
            img = nib.load(labelspath)
            labels = np.asanyarray(img.dataobj)
            affine_labels = img.header.get_zooms()[:3]
        else:
            labels, affine_labels, _, _,_ = load_nifti_remote(labelspath, sftp)

        #labels, affine_labels = load_nifti(labelspath)
        if verbose:
            print("Label mask taken from " + labelspath)
    else:
        txt=f"Mask for subject {subject} not found"
        raise Exception(txt)

    return labels, affine_labels, labelspath

def getlabeltypemask(mypath, subject, ROI_legends, labeltype = '', verbose=False, sftp=None):

    labelmask, labelaffine, labelpath = getlabelmask(mypath, subject, verbose=verbose, sftp=sftp)

    converter_lr, converter_comb, index_to_struct_lr, index_to_struct_comb = atlas_converter(ROI_legends)
    if labeltype == 'combined':
        labeloutpath = labelpath.replace('.nii.gz', '_comb.nii.gz')
        if not os.path.isfile(labeloutpath):
            labelmask = convert_labelmask(labelmask, converter_comb, atlas_outpath=labeloutpath,
                                          affine_labels=labelaffine)
        else:
            labelmask, labelaffine = load_nifti(labeloutpath)
        index_to_struct = index_to_struct_comb
    elif labeltype == 'lrordered':
        labeloutpath = labelpath.replace('.nii.gz', '_lr_ordered.nii.gz')
        if not checkfile_exists_remote(labeloutpath,sftp):
            labelmask = convert_labelmask(labelmask, converter_lr, atlas_outpath=labeloutpath,
                                          affine_labels=labelaffine)
        else:
            labelmask, labelaffine, _, _, _ = load_nifti_remote(labeloutpath, sftp)
        index_to_struct = index_to_struct_lr
    else:
        labeloutpath = labelpath
        index_to_struct = None

    if verbose:
        print(f'Final label taken from {labeloutpath}')
    return labelmask, labelaffine, labeloutpath, index_to_struct


def getmask_old(mypath, subject, masktype = "subjspace", verbose=None):
    if os.path.isfile(mypath):
        if mypath.contains(masktype+'binary_mask.nii.gz'):
            mask, affine_mask = load_nifti(mypath)
            if verbose:
                print("Mask taken from " + mypath)
            return(mask, affine_mask)
        else:
            mypath = str(pathlib.Path(mypath).parent.absolute())
    subjectdir = glob.glob(os.path.join(mypath, "*" + subject + "*"))
    if np.size(subjectdir) == 1:
        mypath = subjectdir[0]
    maskpath = glob.glob(os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz'))
    if np.size(maskpath)>0:
        maskpath = maskpath[0]

    if np.size(maskpath) == 1:
        mask, affine_mask = load_nifti(maskpath)
        if verbose:
            print("Mask taken from " + maskpath)
        return mask, affine_mask
    elif np.size(maskpath) == 0:
        print(f"mask not found {os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz')}")
        raise Exception(f"here is what is going on {os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz')}")
        return None, None
    elif np.size(maskpath)>1:
        raise Warning("too many masks fitting parameters!!")


def getmask(mypath, subject, masktype = "subjspace", verbose=None, sftp=None):
    list_options = [os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz'), os.path.join(mypath,subject + '*_mask*.nii.gz')]

    if sftp is None:
        for list_option in list_options:
            if '*' in list_option:
                maskoption = glob.glob(list_option)
                if np.size(maskoption) > 0:
                    maskpath = maskoption[0]
                    break
            else:
                if os.path.exists(list_option):
                    maskpath = list_option
    else:
        for list_option in list_options:
            maskoption = glob_remote(list_option, sftp)
            if np.size(maskoption) > 0 and '.nii' in maskoption[0]:
                if np.size(maskoption)>1:
                    txt = f'Too many masks fitting the parameters for subject{subject}'
                    raise Exception(txt)
                maskpath = maskoption[0]
                break

    if 'maskpath' in locals():
        if sftp is None:
            mask, affine_mask = load_nifti(maskpath)
        else:
            mask, affine_mask,_,_,_ = load_nifti_remote(maskpath, sftp)
        if verbose:
            print("Mask taken from " + maskpath)
        return mask, affine_mask
    else:
        print(f"mask not found {os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz')}")
        return None, None
        #raise Exception(
        #    f"here is what is going on {os.path.join(mypath, subject + '*' + masktype + '*_mask*.nii.gz')}")



def get_diff_ref(label_folder, subject, ref,sftp=None):
    diff_path = os.path.join(label_folder,f'{subject}_{ref}_to_MDT.nii.gz')
    if checkfile_exists_remote(diff_path, sftp):
        return diff_path
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), diff_path)

def move_bvals(mypath, subject, diffpathnew):


    subjfolder = glob.glob(os.path.join(mypath, "*" + subject + "*/"))
    if np.size(subjfolder) == 1 and os.path.isdir(subjfolder[0]):
        subjfolder = subjfolder[0]
        if np.size(glob.glob(os.path.join(subjfolder,"*nii*"))) > 0:
            mypath = subjfolder
    elif np.size(subjfolder) > 1:
        raise Warning
    elif np.size(glob.glob(os.path.join(mypath,subject+"*rawnii*"))) > 0:
        fdiffpath = (glob.glob(os.path.join(mypath,subject+"*rawnii*")))[0]
    elif np.size(glob.glob(os.path.join(mypath,subject+"*dwi*nii*"))) > 0:
        fdiffpath = (glob.glob(os.path.join(mypath,subject+"*dwi*nii*")))[0]
    elif os.path.exists(os.path.join(mypath,subject+"_dwi.nii.gz")):
        fdiffpath = (os.path.join(mypath,subject+"_dwi.nii.gz"))
    elif os.path.exists(mypath + '/Reg_' + subject + '_nii4D.nii.gz'):
        fdiffpath = mypath + '/Reg_' + subject + '_nii4D.nii.gz'
    elif os.path.exists(mypath + '/nii4D_' + subject + '.nii'):
        fdiffpath = mypath + '/nii4D_' + subject + '.nii'
    elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
        fdiffpath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'):
        fdiffpath = mypath + '/4Dnii/'+subject+'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/'+subject+'_nii4D_RAS.nii.gz'):
        fdiffpath = mypath + '/'+subject+'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/' + subject + '/'):
        fdiffpath = glob.glob(mypath + '/' + subject + '/' + subject + '*nii4D*.nii*')[0]

    if os.path.isfile(mypath):
        mypath = str(pathlib.Path(fdiffpath).parent.absolute())
    elif os.path.isdir(mypath):
        mypath = mypath

    fbvals_new = os.path.join(diffpathnew, subject + "_bvals_fix.txt")
    fbvec_new = os.path.join(diffpathnew, subject + "_bvec_fix.txt")

    if not os.path.exists(fbvals_new) and not os.path.exists(fbvec_new):
        try:
            fbvals = glob.glob(os.path.join(mypath, subject + '*_bvals_fix.txt'))[0]
            fbvecs = glob.glob(os.path.join(mypath, subject + '*_bvec_fix.txt'))[0]
        except IndexError:
            fbvals = glob.glob(mypath + '/' + subject + '*_bvals.txt')[0]
            fbvecs = glob.glob(mypath + '/' + subject + '*_bvec*.txt')[0]
            fbvals, fbvecs = fix_bvals_bvecs(fbvals,fbvecs)
        shutil.copyfile(fbvals, fbvals_new)
        shutil.copyfile(fbvecs, fbvec_new)

    return fbvals_new, fbvec_new
