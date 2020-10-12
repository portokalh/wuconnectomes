import glob
import numpy as np


def get_tract_params(mypath, subject, str_identifier, verbose):

    trkpath = gettrkpath(mypath, subject, str_identifier, verbose)
    trkdata = load_trk(trkpath, "same")
    verbose = True
    if verbose:
        print("loaded ")
    # trkdata.to_vox()
    header = trkdata.space_attribute
    affine = trkdata._affine
    lengths = length(trkdata.streamlines)
    del trkdata
    # lengths = list(length(trkstreamlines))
    lengths = list(lengths)
    numtracts = np.size(lengths)
    minlength = np.min(lengths)
    maxlength = np.max(lengths)
    meanlength = np.mean(lengths)
    stdlength = np.std(lengths)
    if verbose:
        print \
            ("For subject " + subject + " the number of tracts is " + numbtracts + ", the minimum length is " + minlength + ", the maximum length is " + maxlength + ", the mean length is " + meanlength + ", the std is " + stdlength)
    return subject, numtracts, minlength, maxlength, meanlength, stdlength, header, affine


def get_anat(mypath, subject, verbose=None):
    filepath =(mypath + '/' + subject + '*_nii4D*.nii*')
    anatpath = glob.glob(filepath)[0]
    return(anatpath)


def get_niftis(mypath, subject, verbose = None):
    filepaths =(mypath + '/' + subject + '*.nii*')
    filepaths = glob.glob(filepaths)
    return filepaths


def get_trks(mypath, subject, verbose = None):
    filepaths = (mypath + '/' + subject + '*.trk*')
    filepaths = glob.glob(filepaths)
    return filepaths


def gettrkpath(trkpath, subject, str_identifier, verbose=False):
    #filepath=(trkpath + '/' + subject + '_wholebrain_' + tractsize + strproperty + 'stepsize_' + str(stepsize) + '.trk')
    filepath=(trkpath + '/' + subject + str_identifier + '.trk')
    trkpaths = glob.glob(filepath)
    if trkpaths:
        trkfile = trkpaths[0]
        if verbose:
            print("Subject " + subject + " was found at " + trkfile)
    else:
        print("Could not find "+filepath)
        return
    return trkfile


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


def getlabelmask(mypath, subject, bvec_orient=[1, 2, 3], verbose=None):

    # ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    if os.path.exists(mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/Reg_' + subject + '_nii4D_brain_mask.nii.gz')
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/' + subject + '_chass_symmetric3_labels_RAS.nii.gz')
    elif os.path.exists(mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/' + subject + '_chass_symmetric3_labels_RAS_combined.nii.gz')
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii.gz')
    elif os.path.exists(mypath + '/mask.nii'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii')
    else:
        print('mask not found, taking all non null values in nii file instead (not recommended for complex operations)')
        labels = np.ones(fdwi_data.shape[0:3])
        affine_labels = affine

    # Build Brain Mask
    # bm = np.where(labels == 0, False, True)
    # mask = bm

    return labels, affine_labels

def getdwidata(mypath, subject, bvec_orient=[1 ,2 ,3], verbose=None):

    # fdwi = mypath + '4Dnii/' + subject + '_nii4D_RAS.nii.gz'
    # fdwipath = mypath + '/nii4D_' + subject + '.nii'
    if os.path.exists(mypath + '/Reg_' + subject + '_nii4D.nii.gz'):
        fdwipath = mypath + '/Reg_' + subject + '_nii4D.nii.gz'
    elif os.path.exists(mypath + '/nii4D_' + subject + '.nii'):
        fdwipath = mypath + '/nii4D_' + subject + '.nii'
    elif os.path.exists(mypath + '/ ' +subject +'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/ ' +subject +'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/4Dnii/ ' +subject +'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/4Dnii/ ' +subject +'_nii4D_RAS.nii.gz'
    elif os.path.exists(mypath + '/ ' +subject +'_nii4D_RAS.nii.gz'):
        fdwipath = mypath + '/ ' +subject +'_nii4D_RAS.nii.gz'
    # fdwi_data, affine, vox_size = load_nifti(fdwipath, return_voxsize=True)

    if verbose:
        txt = "Extracting information from the dwi file located at " + fdwipath
        print(txt)
        send_mail(txt ,subject="Begin data extraction")

    if 'fdwipath' not in locals():
        txt = "The subject " + subject + " was not detected, exit"
        print(txt)
        send_mail(txt ,subject="Error")
        return(0 ,0 ,0 ,0 ,0 ,0 ,0 ,0)

    img = nib.load(fdwipath)
    fdwi_data = img.get_data()
    vox_size = img.header.get_zooms()[:3]
    affine = img.affine
    hdr = img.header
    header = get_reference_info(fdwipath)
    del(img)

    # ffalabels = mypath + 'labels/' + 'fa_labels_warp_' + subject + '_RAS.nii.gz'

    if os.path.exists(mypath +'/Reg_ ' +subject +'_nii4D_brain_mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath +'/Reg_ ' +subject +'_nii4D_brain_mask.nii.gz')
    elif os.path.exists(mypath +'/ ' +subject +'_chass_symmetric3_labels_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath +'/ ' +subject +'_chass_symmetric3_labels_RAS.nii.gz')
    elif os.path.exists(mypath +'/ ' +subject +'_chass_symmetric3_labels_RAS_combined.nii.gz'):
        labels, affine_labels = load_nifti(mypath +'/ ' +subject +'_chass_symmetric3_labels_RAS_combined.nii.gz')
    elif os.path.exists(mypath + '/fa_labels_warp_' + subject +'_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/labels/fa_labels_warp_' + subject + '_RAS.nii.gz')
    elif os.path.exists(mypath + '/mask.nii.gz'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii.gz')
    elif os.path.exists(mypath + '/mask.nii'):
        labels, affine_labels = load_nifti(mypath + '/mask.nii')
    else:
        print('mask not found, taking all non null values in nii file instead (not recommended for complex operations)')
        labels = np.ones(fdwi_data.shape[0:3])
        affine_labels = affine

    # fbvals = mypath + '4Dnii/' + subject + '_RAS_ecc_bvals.txt'
    # fbvecs = mypath + '4Dnii/' + subject + '_RAS_ecc_bvecs.txt'
    # fbvals = glob.glob(mypath + '*/*' + subject + '*_bval*.txt')[0]
    # fbvecs = glob.glob(mypath + '*/*' + subject + '*_bvec*.txt')[0]
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
    bvecs = np.c_[
        bvec_sign[0] * bvecs[:, np.abs(bvec_orient[0]) - 1], bvec_sign[1] * bvecs[:, np.abs(bvec_orient[1]) - 1],
        bvec_sign[2] * bvecs[:, np.abs(bvec_orient[2]) - 1]]

    # bvecs = np.c_[bvecs[:, 1], bvecs[:, 0], -bvecs[:, 2]]
    # bvecs = np.c_[-bvecs[:, 1], bvecs[:, 0], bvecs[:, 2]]

    gtab = gradient_table(bvals, bvecs)

    # Build Brain Mask
    # bm = np.where(labels == 0, False, True)
    # mask = bm

    return fdwi_data, affine, gtab, labels, vox_size, fdwipath, hdr, header
