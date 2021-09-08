
"""
Created by Jacques Stout
Part of the DTC pipeline
Just a variant of the dipy denoiser methods, but small modifications give it a multiprocessing option
"""

import numpy as np
import multiprocessing as mp
from time import time

try:
    from scipy.linalg.lapack import dgesvd as svd
    svd_args = [1, 0]
    # If you have an older version of scipy, we fall back
    # on the standard scipy SVD API:
except ImportError:
    from scipy.linalg import svd
    svd_args = [False]
from scipy.linalg import eigh


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def _pca_classifier(L, nvoxels):
    """ Classifies which PCA eigenvalues are related to noise and estimates the
    noise variance

    Parameters
    ----------
    L : array (n,)
        Array containing the PCA eigenvalues in ascending order.
    nvoxels : int
        Number of voxels used to compute L

    Returns
    -------
    var : float
        Estimation of the noise variance
    ncomps : int
        Number of eigenvalues related to noise

    Notes
    -----
    This is based on the algorithm described in [1]_.

    References
    ----------
    .. [1] Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers,
           Fieremans E, 2016. Denoising of Diffusion MRI using random matrix
           theory. Neuroimage 142:394-406.
           doi: 10.1016/j.neuroimage.2016.08.016
    """
    var = np.mean(L)
    c = L.size - 1
    r = L[c] - L[0] - 4 * np.sqrt((c + 1.0) / nvoxels) * var
    while r > 0:
        var = np.mean(L[:c])
        c = c - 1
        r = L[c] - L[0] - 4 * np.sqrt((c + 1.0) / nvoxels) * var
    ncomps = c + 1
    return var, ncomps

def pca_patchloop(patch_radius, arr, arr_shape, mask, jx1, tau_factor, dim, is_svd, var, calc_sigma=True, verbose=False):

    patch_size = 2 * patch_radius + 1
    sizei = arr_shape[0] - 2 * patch_radius
    Xesti = np.zeros((sizei, patch_size, patch_size, patch_size, dim))
    this_thetai = np.zeros(sizei)
    this_vari = np.zeros(sizei)
    for ix1 in range(arr_shape[0] - 2 * patch_radius):

        if mask[ix1+patch_radius]:

            #X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
            #                patch_size ** 3, dim)
            ix2 = ix1 + 2 * patch_radius + 1

            X = arr[ix1:ix2, :, :].reshape(
                patch_size ** 3, dim)
            # compute the mean and normalize
            M = np.mean(X, axis=0)
            # Upcast the dtype for precision in the SVD
            X = X - M
            patch_size = 2 * patch_radius + 1

            if is_svd:
                # PCA using an SVD
                U, S, Vt = svd(X, *svd_args)[:3]
                # Items in S are the eigenvalues, but in ascending order
                # We invert the order (=> descending), square and normalize
                # \lambda_i = s_i^2 / n
                d = S[::-1] ** 2 / X.shape[0]
                # Rows of Vt are eigenvectors, but also in ascending
                # eigenvalue order:
                W = Vt[::-1].T

            else:
                # PCA using an Eigenvalue decomposition
                C = np.transpose(X).dot(X)
                C = C / X.shape[0]
                [d, W] = eigh(C, turbo=True)

            if calc_sigma:
                # Random matrix theory
                this_vari[ix1], ncomps = _pca_classifier(d, patch_size ** 3)
            else:
                # Predefined variance
                this_vari[ix1] = var[ix1 + patch_radius]

            # Threshold by tau:
            tau = tau_factor ** 2 * this_vari[ix1]

            # Update ncomps according to tau_factor
            ncomps = np.sum(d < tau)
            W[:, :ncomps] = 0

            # This is equations 1 and 2 in Manjon 2013:
            Xest = X.dot(W).dot(W.T) + M
            Xesti[ix1,:,:,:,:] = Xest.reshape(patch_size,
                                patch_size,
                                patch_size, dim)
            # This is equation 3 in Manjon 2013:
            this_thetai[ix1] = 1.0 / (1.0 + dim - ncomps)
            """""
            theta[ix1:ix2, jx1:jx2, kx1:kx2] = this_theta
            thetax[ix1:ix2, jx1:jx2, kx1:kx2] = Xest * this_theta
            if calc_sigma:
                var[ix1:ix2, jx1:jx2, kx1:kx2] = this_var * this_theta
                thetavar[ix1:ix2, jx1:jx2, kx1:kx2] = this_theta
            else:
                var = 0
                thetavar = 0
            """

            #return theta,thetax,var,thetavar
        else:
            Xesti[ix1] = 0
            this_thetai[ix1] = 0
            this_vari[ix1] = 0
    return [Xesti, this_thetai, this_vari, jx1]



def genpca_parallel(arr, sigma=None, mask=None, patch_radius=2, pca_method='eig',
           tau_factor=None, return_sigma=False, out_dtype=None, processes=1, verbose=False):
    """General function to perform PCA-based denoising of diffusion datasets.

    Parameters
    ----------
    arr : 4D array
        Array of data to be denoised. The dimensions are (X, Y, Z, N), where N
        are the diffusion gradient directions.
    sigma : float or 3D array (optional)
        Standard deviation of the noise estimated from the data. If no sigma
        is given, this will be estimated based on random matrix theory
        [1]_,[2]_
    mask : 3D boolean array (optional)
        A mask with voxels that are true inside the brain and false outside of
        it. The function denoises within the true part and returns zeros
        outside of those voxels.
    patch_radius : int (optional)
        The radius of the local patch to be taken around each voxel (in
        voxels). Default: 2 (denoise in blocks of 5x5x5 voxels).
    pca_method : 'eig' or 'svd' (optional)
        Use either eigenvalue decomposition (eig) or singular value
        decomposition (svd) for principal component analysis. The default
        method is 'eig' which is faster. However, occasionally 'svd' might be
        more accurate.
    tau_factor : float (optional)
        Thresholding of PCA eigenvalues is done by nulling out eigenvalues that
        are smaller than:

        .. math ::

                \tau = (\tau_{factor} \sigma)^2

        \tau_{factor} can be set to a predefined values (e.g. \tau_{factor} =
        2.3 [3]_), or automatically calculated using random matrix theory
        (in case that \tau_{factor} is set to None).
        Default: None.
    return_sigma : bool (optional)
        If true, the Standard deviation of the noise will be returned.
        Default: False.
    out_dtype : str or dtype (optional)
        The dtype for the output array. Default: output has the same dtype as
        the input.

    Returns
    -------
    denoised_arr : 4D array
        This is the denoised array of the same size as that of the input data,
        clipped to non-negative values

    References
    ----------
    .. [1] Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers,
           Fieremans E, 2016. Denoising of Diffusion MRI using random matrix
           theory. Neuroimage 142:394-406.
           doi: 10.1016/j.neuroimage.2016.08.016
    .. [2] Veraart J, Fieremans E, Novikov DS. 2016. Diffusion MRI noise
           mapping using random matrix theory. Magnetic Resonance in Medicine.
           doi: 10.1002/mrm.26059.
    .. [3] Manjon JV, Coupe P, Concha L, Buades A, Collins DL (2013)
           Diffusion Weighted Image Denoising Using Overcomplete Local
           PCA. PLoS ONE 8(9): e73021.
           https://doi.org/10.1371/journal.pone.0073021
    """
    if mask is None:
        # If mask is not specified, use the whole volume
        mask = np.ones_like(arr, dtype=bool)[..., 0]

    if out_dtype is None:
        out_dtype = arr.dtype

    # We retain float64 precision, iff the input is in this precision:
    if arr.dtype == np.float64:
        calc_dtype = np.float64
    # Otherwise, we'll calculate things in float32 (saving memory)
    else:
        calc_dtype = np.float32

    if not arr.ndim == 4:
        raise ValueError("PCA denoising can only be performed on 4D arrays.",
                         arr.shape)

    if pca_method.lower() == 'svd':
        is_svd = True
    elif pca_method.lower() == 'eig':
        is_svd = False
    else:
        raise ValueError("pca_method should be either 'eig' or 'svd'")

    patch_size = 2 * patch_radius + 1

    if return_sigma is True or sigma is None:
        calc_sigma = True
    else:
        calc_sigma = False

    if patch_size ** 3 < arr.shape[-1]:
        e_s = "You asked for PCA denoising with a "
        e_s += "patch_radius of {0} ".format(patch_radius)
        e_s += "for data with {0} directions. ".format(arr.shape[-1])
        e_s += "This would result in an ill-conditioned PCA matrix. "
        e_s += "Please increase the patch_radius."
        raise ValueError(e_s)

    if isinstance(sigma, np.ndarray):
        var = sigma ** 2
        if not sigma.shape == arr.shape[:-1]:
            e_s = "You provided a sigma array with a shape"
            e_s += "{0} for data with".format(sigma.shape)
            e_s += "shape {0}. Please provide a sigma array".format(arr.shape)
            e_s += " that matches the spatial dimensions of the data."
            raise ValueError(e_s)
    elif isinstance(sigma, (int, float)):
        var = sigma ** 2 * np.ones(arr.shape[:-1])

    dim = arr.shape[-1]
    if tau_factor is None:
        tau_factor = 1 + np.sqrt(dim / (patch_size ** 3))

    pool = mp.Pool(processes)

    #for k in range(patch_radius, arr.shape[2] - patch_radius):
    #    mykrange[i] = k
    #    i += 1

    theta = np.zeros(arr.shape, dtype=calc_dtype)
    thetax = np.zeros(arr.shape, dtype=calc_dtype)

    if calc_sigma:
        var = np.zeros(arr.shape[:-1], dtype=calc_dtype)
        thetavar = np.zeros(arr.shape[:-1], dtype=calc_dtype)

    allnum = (arr.shape[2] - 2 * patch_radius) * (arr.shape[1] - 2 * patch_radius) * (arr.shape[0] - 2 * patch_radius)

    duration1=time()
    if verbose:
        print("Begin mpca denoising")
    arr_shape = arr.shape
    for kx1 in range(0, arr_shape[2] - 2 * patch_radius):
        """
        kx1 = np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))
        jx1 = np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius))
        ix1 = np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))%(arr_shape[0] - 2 * patch_radius))
        jx2 = jx1 + 2 * patch_radius + 1
        ix2 = ix1 + 2 * patch_radius + 1
        """
        kx2 = kx1 + 2*patch_radius + 1

        jlist= list(range(arr_shape[1] - 2 * patch_radius))

        # Shorthand for indexing variables:

        #if not mask[ix1 + patch_radius, jx1 + patch_radius, kx1 + patch_radius]:

        #    X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
        #                    patch_size ** 3, dim)
        #minimask=mask[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
        #                    patch_size ** 3, dim)
        resultslist=[]
        resultslist = pool.starmap_async(pca_patchloop,
                                         [(patch_radius, arr[:, jx1:jx1 + 2*patch_radius + 1, kx1:kx2],
                                           arr_shape,mask[:, jx1+patch_radius, kx1+patch_radius],
                                           jx1, tau_factor, dim, is_svd, var[:, jx1+patch_radius, kx1+patch_radius],
                                           calc_sigma, verbose) for jx1 in jlist]).get()

        #Xest=resultslist[0][0]
        #this_theta=resultslist[0][1]
        #this_var=resultslist[0][2]
        #resultslist[0][3]
        #Xest, this_theta, this_var, numlist = resultslist

        for jj in range(len(jlist)):
            jx1 = resultslist[jj][3]
            jx2 = jx1 + 2 * patch_radius + 1
            for ix1 in range(arr_shape[0] - 2 * patch_radius):
                ix2 = ix1 + 2 * patch_radius + 1
                theta[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[jj][1][ix1]
                thetax[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[jj][0][ix1] * resultslist[jj][1][ix1]
                if calc_sigma:
                    var[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[jj][2][ix1] * resultslist[jj][1][ix1]
                    thetavar[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[jj][1][ix1]

        if verbose:
            print("finished " + str(kx1) + " of " + str(arr_shape[2] - 2 * patch_radius) )
            print("Process has been running for "+ str(time()-duration1) + "s")
            #pool.starmap_async(create_tracts, [(mypath, outpath, subject, step_size, function_processes,
            #                                                                        saved_streamlines, denoise, savefa, verbose) for subject in l]).get()
    #theta, thetax, var, thetavar = pca_patchloop(patch_radius, arr, k)


    #for k in range(patch_radius, arr.shape[2] - patch_radius):
    #    for j in range(patch_radius, arr.shape[1] - patch_radius):
    #        for i in range(patch_radius, arr.shape[0] - patch_radius):
    #mynumrange = range((arr.shape[2] - 2 * patch_radius) * (arr.shape[1] - 2 * patch_radius) *
    #                   (arr.shape[0] - 2 * patch_radius))
    #Xest, this_theta, this_var, ix1l, ix2l, jx1l, jx2l, kx1l, kx2l = pool.starmap_async(pca_patchloop_jusk, [(patch_radius,arr,num, is_svd) for num in mynumrange]).get()
    #print("hi")
    #print(mykrange)
    #arr=np.zeros(5)
    #thetap, thetaxp, varp, thetavarp = pool.starmap_async(pca_patchloop_justk, [(patch_radius, arr, k, is_svd,
    #                                                                             calc_sigma, calc_dtype, verbose) for k in mykrange]).get()
    if verbose:
        print("finished main computations, preparing matrix")

    denoised_arr = thetax / theta
    denoised_arr.clip(min=0, out=denoised_arr)
    denoised_arr[mask == 0] = 0
    if verbose:
        print("finished calculating denoised matrix")
    if return_sigma is True:
        if sigma is None:
            var = var / thetavar
            var[mask == 0] = 0
            return denoised_arr.astype(out_dtype), np.sqrt(var)
        else:
            return denoised_arr.astype(out_dtype), sigma
    else:
        return denoised_arr.astype(out_dtype)



def genpca(arr, sigma=None, mask=None, patch_radius=2, pca_method='eig',
           tau_factor=None, return_sigma=False, out_dtype=None, verbose = False):
    r"""General function to perform PCA-based denoising of diffusion datasets.

    Parameters
    ----------
    arr : 4D array
        Array of data to be denoised. The dimensions are (X, Y, Z, N), where N
        are the diffusion gradient directions.
    sigma : float or 3D array (optional)
        Standard deviation of the noise estimated from the data. If no sigma
        is given, this will be estimated based on random matrix theory
        [1]_,[2]_
    mask : 3D boolean array (optional)
        A mask with voxels that are true inside the brain and false outside of
        it. The function denoises within the true part and returns zeros
        outside of those voxels.
    patch_radius : int (optional)
        The radius of the local patch to be taken around each voxel (in
        voxels). Default: 2 (denoise in blocks of 5x5x5 voxels).
    pca_method : 'eig' or 'svd' (optional)
        Use either eigenvalue decomposition (eig) or singular value
        decomposition (svd) for principal component analysis. The default
        method is 'eig' which is faster. However, occasionally 'svd' might be
        more accurate.
    tau_factor : float (optional)
        Thresholding of PCA eigenvalues is done by nulling out eigenvalues that
        are smaller than:

        .. math ::

                \tau = (\tau_{factor} \sigma)^2

        \tau_{factor} can be set to a predefined values (e.g. \tau_{factor} =
        2.3 [3]_), or automatically calculated using random matrix theory
        (in case that \tau_{factor} is set to None).
        Default: None.
    return_sigma : bool (optional)
        If true, the Standard deviation of the noise will be returned.
        Default: False.
    out_dtype : str or dtype (optional)
        The dtype for the output array. Default: output has the same dtype as
        the input.

    Returns
    -------
    denoised_arr : 4D array
        This is the denoised array of the same size as that of the input data,
        clipped to non-negative values

    References
    ----------
    .. [1] Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers,
           Fieremans E, 2016. Denoising of Diffusion MRI using random matrix
           theory. Neuroimage 142:394-406.
           doi: 10.1016/j.neuroimage.2016.08.016
    .. [2] Veraart J, Fieremans E, Novikov DS. 2016. Diffusion MRI noise
           mapping using random matrix theory. Magnetic Resonance in Medicine.
           doi: 10.1002/mrm.26059.
    .. [3] Manjon JV, Coupe P, Concha L, Buades A, Collins DL (2013)
           Diffusion Weighted Image Denoising Using Overcomplete Local
           PCA. PLoS ONE 8(9): e73021.
           https://doi.org/10.1371/journal.pone.0073021
    """
    if mask is None:
        # If mask is not specified, use the whole volume
        mask = np.ones_like(arr, dtype=bool)[..., 0]

    if out_dtype is None:
        out_dtype = arr.dtype

    # We retain float64 precision, iff the input is in this precision:
    if arr.dtype == np.float64:
        calc_dtype = np.float64
    # Otherwise, we'll calculate things in float32 (saving memory)
    else:
        calc_dtype = np.float32

    if not arr.ndim == 4:
        raise ValueError("PCA denoising can only be performed on 4D arrays.",
                         arr.shape)

    if pca_method.lower() == 'svd':
        is_svd = True
    elif pca_method.lower() == 'eig':
        is_svd = False
    else:
        raise ValueError("pca_method should be either 'eig' or 'svd'")

    patch_size = 2 * patch_radius + 1

    if patch_size ** 3 < arr.shape[-1]:
        e_s = "You asked for PCA denoising with a "
        e_s += "patch_radius of {0} ".format(patch_radius)
        e_s += "for data with {0} directions. ".format(arr.shape[-1])
        e_s += "This would result in an ill-conditioned PCA matrix. "
        e_s += "Please increase the patch_radius."
        raise ValueError(e_s)

    if isinstance(sigma, np.ndarray):
        var = sigma ** 2
        if not sigma.shape == arr.shape[:-1]:
            e_s = "You provided a sigma array with a shape"
            e_s += "{0} for data with".format(sigma.shape)
            e_s += "shape {0}. Please provide a sigma array".format(arr.shape)
            e_s += " that matches the spatial dimensions of the data."
            raise ValueError(e_s)
    elif isinstance(sigma, (int, float)):
        var = sigma ** 2 * np.ones(arr.shape[:-1])

    dim = arr.shape[-1]
    if tau_factor is None:
        tau_factor = 1 + np.sqrt(dim / (patch_size ** 3))

    theta = np.zeros(arr.shape, dtype=calc_dtype)
    thetax = np.zeros(arr.shape, dtype=calc_dtype)

    if return_sigma is True and sigma is None:
        var = np.zeros(arr.shape[:-1], dtype=calc_dtype)
        thetavar = np.zeros(arr.shape[:-1], dtype=calc_dtype)

    # loop around and find the 3D patch for each direction at each pixel
    duration1=time()
    if verbose:
        print("Start of mpca process")
    for kk in range(patch_radius, arr.shape[2] - patch_radius):
        for jj in range(patch_radius, arr.shape[1] - patch_radius):
            for ii in range(patch_radius, arr.shape[0] - patch_radius):
                # Shorthand for indexing variables:
                if not mask[ii, jj, kk]:
                    continue
                ix1 = ii - patch_radius
                ix2 = ii + patch_radius + 1
                jx1 = jj - patch_radius
                jx2 = jj + patch_radius + 1
                kx1 = kk - patch_radius
                kx2 = kk + patch_radius + 1

                X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
                                patch_size ** 3, dim)
                # compute the mean and normalize
                M = np.mean(X, axis=0)
                # Upcast the dtype for precision in the SVD
                X = X - M

                if is_svd:
                    # PCA using an SVD
                    U, S, Vt = svd(X, *svd_args)[:3]
                    # Items in S are the eigenvalues, but in ascending order
                    # We invert the order (=> descending), square and normalize
                    # \lambda_i = s_i^2 / n
                    d = S[::-1] ** 2 / X.shape[0]
                    # Rows of Vt are eigenvectors, but also in ascending
                    # eigenvalue order:
                    W = Vt[::-1].T

                else:
                    # PCA using an Eigenvalue decomposition
                    C = np.transpose(X).dot(X)
                    C = C / X.shape[0]
                    [d, W] = eigh(C, turbo=True)

                if sigma is None:
                    # Random matrix theory
                    this_var, ncomps = _pca_classifier(d, patch_size ** 3)
                else:
                    # Predefined variance
                    this_var = var[ii, jj, kk]

                # Threshold by tau:
                tau = tau_factor ** 2 * this_var

                # Update ncomps according to tau_factor
                ncomps = np.sum(d < tau)
                W[:, :ncomps] = 0

                # This is equations 1 and 2 in Manjon 2013:
                Xest = X.dot(W).dot(W.T) + M
                Xest = Xest.reshape(patch_size,
                                    patch_size,
                                    patch_size, dim)
                # This is equation 3 in Manjon 2013:
                this_theta = 1.0 / (1.0 + dim - ncomps)
                theta[ix1:ix2, jx1:jx2, kx1:kx2] += this_theta
                thetax[ix1:ix2, jx1:jx2, kx1:kx2] += Xest * this_theta
                if return_sigma is True and sigma is None:
                    var[ix1:ix2, jx1:jx2, kx1:kx2] += this_var * this_theta
                    thetavar[ix1:ix2, jx1:jx2, kx1:kx2] += this_theta

        if verbose:
            print("Ran loop on line ", str(kk))
            print("Process has been running for "+ str(time()-duration1) + "s")

    denoised_arr = thetax / theta
    denoised_arr.clip(min=0, out=denoised_arr)
    denoised_arr[mask == 0] = 0
    if return_sigma is True:
        if sigma is None:
            var = var / thetavar
            var[mask == 0] = 0
            return denoised_arr.astype(out_dtype), np.sqrt(var)
        else:
            return denoised_arr.astype(out_dtype), sigma
    else:
        return denoised_arr.astype(out_dtype)


def localpca(arr, sigma, mask=None, patch_radius=2, pca_method='eig',
             tau_factor=2.3, processes = 1, out_dtype=None, verbose=False):
    r""" Performs local PCA denoising according to Manjon et al. [1]_.

    Parameters
    ----------
    arr : 4D array
        Array of data to be denoised. The dimensions are (X, Y, Z, N), where N
        are the diffusion gradient directions.
    sigma : float or 3D array
        Standard deviation of the noise estimated from the data.
    mask : 3D boolean array (optional)
        A mask with voxels that are true inside the brain and false outside of
        it. The function denoises within the true part and returns zeros
        outside of those voxels.
    patch_radius : int (optional)
        The radius of the local patch to be taken around each voxel (in
        voxels). Default: 2 (denoise in blocks of 5x5x5 voxels).
    pca_method : 'eig' or 'svd' (optional)
        Use either eigenvalue decomposition (eig) or singular value
        decomposition (svd) for principal component analysis. The default
        method is 'eig' which is faster. However, occasionally 'svd' might be
        more accurate.
    tau_factor : float (optional)
        Thresholding of PCA eigenvalues is done by nulling out eigenvalues that
        are smaller than:

        .. math ::

                \tau = (\tau_{factor} \sigma)^2

        \tau_{factor} can be change to adjust the relationship between the
        noise standard deviation and the threshold \tau. If \tau_{factor} is
        set to None, it will be automatically calculated using the
        Marcenko-Pastur distribution [2]_.
        Default: 2.3 (according to [1]_)
    out_dtype : str or dtype (optional)
        The dtype for the output array. Default: output has the same dtype as
        the input.

    Returns
    -------
    denoised_arr : 4D array
        This is the denoised array of the same size as that of the input data,
        clipped to non-negative values

    References
    ----------
    .. [1] Manjon JV, Coupe P, Concha L, Buades A, Collins DL (2013)
           Diffusion Weighted Image Denoising Using Overcomplete Local
           PCA. PLoS ONE 8(9): e73021.
           https://doi.org/10.1371/journal.pone.0073021
    .. [2] Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers,
           Fieremans E, 2016. Denoising of Diffusion MRI using random matrix
           theory. Neuroimage 142:394-406.
           doi: 10.1016/j.neuroimage.2016.08.016
    """

    if processes == 1:
        return genpca(arr, sigma=sigma, mask=mask, patch_radius=patch_radius,
                  pca_method=pca_method,  tau_factor=tau_factor,
                  return_sigma=False, out_dtype=out_dtype, verbose=verbose)
    elif processes > 1:
        return genpca_parallel(arr, sigma=sigma, mask=mask, patch_radius=patch_radius,
                               pca_method=pca_method,  tau_factor=tau_factor,
                               return_sigma=False, out_dtype=out_dtype, processes=processes, verbose=verbose)
    else:
        print("unrecognized number of processes, run as standard genpca")
        return genpca(arr, sigma=sigma, mask=mask, patch_radius=patch_radius,
                      pca_method=pca_method, tau_factor=tau_factor,
                      return_sigma=False, out_dtype=out_dtype, verbose=verbose)


def mppca(arr, mask=None, patch_radius=2, pca_method='eig',
          return_sigma=False, out_dtype=None, processes=1, verbose = False):
    r"""Performs PCA-based denoising using the Marcenko-Pastur
    distribution [1]_.

    Parameters
    ----------
    arr : 4D array
        Array of data to be denoised. The dimensions are (X, Y, Z, N), where N
        are the diffusion gradient directions.
    mask : 3D boolean array (optional)
        A mask with voxels that are true inside the brain and false outside of
        it. The function denoises within the true part and returns zeros
        outside of those voxels.
    patch_radius : int (optional)
        The radius of the local patch to be taken around each voxel (in
        voxels). Default: 2 (denoise in blocks of 5x5x5 voxels).
    pca_method : 'eig' or 'svd' (optional)
        Use either eigenvalue decomposition (eig) or singular value
        decomposition (svd) for principal component analysis. The default
        method is 'eig' which is faster. However, occasionally 'svd' might be
        more accurate.
    return_sigma : bool (optional)
        If true, a noise standard deviation estimate based on the
        Marcenko-Pastur distribution is returned [2]_.
        Default: False.
    out_dtype : str or dtype (optional)
        The dtype for the output array. Default: output has the same dtype as
        the input.

    Returns
    -------
    denoised_arr : 4D array
        This is the denoised array of the same size as that of the input data,
        clipped to non-negative values
    sigma : 3D array (when return_sigma=True)
        Estimate of the spatial varying standard deviation of the noise

    References
    ----------
    .. [1] Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers,
           Fieremans E, 2016. Denoising of Diffusion MRI using random matrix
           theory. Neuroimage 142:394-406.
           doi: 10.1016/j.neuroimage.2016.08.016
    .. [2] Veraart J, Fieremans E, Novikov DS. 2016. Diffusion MRI noise
           mapping using random matrix theory. Magnetic Resonance in Medicine.
           doi: 10.1002/mrm.26059.
    """

    if processes == 1:
        return genpca(arr, sigma=None, mask=mask, patch_radius=patch_radius,
                  pca_method=pca_method, tau_factor=None,
                  return_sigma=return_sigma, out_dtype=out_dtype, verbose=verbose)
    elif processes > 1:
        return genpca_parallel(arr, sigma=None, mask=mask, patch_radius=patch_radius,
                               pca_method=pca_method, tau_factor=None,
                               return_sigma=return_sigma, out_dtype=out_dtype, processes=processes, verbose=verbose)
    else:
        print("unrecognized number of processes, run as standard genpca")
        return genpca(arr, sigma=None, mask=mask, patch_radius=patch_radius,
                      pca_method=pca_method, tau_factor=None,
                      return_sigma=return_sigma, out_dtype=out_dtype, verbose=verbose)

# End

"""
earlier version
def genpca_parallel(arr, sigma=None, mask=None, patch_radius=2, pca_method='eig',
           tau_factor=None, return_sigma=False, out_dtype=None, processes=1, verbose=False):

    if mask is None:
        # If mask is not specified, use the whole volume
        mask = np.ones_like(arr, dtype=bool)[..., 0]

    if out_dtype is None:
        out_dtype = arr.dtype

    # We retain float64 precision, iff the input is in this precision:
    if arr.dtype == np.float64:
        calc_dtype = np.float64
    # Otherwise, we'll calculate things in float32 (saving memory)
    else:
        calc_dtype = np.float32

    if not arr.ndim == 4:
        raise ValueError("PCA denoising can only be performed on 4D arrays.",
                         arr.shape)

    if pca_method.lower() == 'svd':
        is_svd = True
    elif pca_method.lower() == 'eig':
        is_svd = False
    else:
        raise ValueError("pca_method should be either 'eig' or 'svd'")

    patch_size = 2 * patch_radius + 1

    if return_sigma is True and sigma is None:
        calc_sigma = True
    else:
        calc_sigma = False

    if patch_size ** 3 < arr.shape[-1]:
        e_s = "You asked for PCA denoising with a "
        e_s += "patch_radius of {0} ".format(patch_radius)
        e_s += "for data with {0} directions. ".format(arr.shape[-1])
        e_s += "This would result in an ill-conditioned PCA matrix. "
        e_s += "Please increase the patch_radius."
        raise ValueError(e_s)

    if isinstance(sigma, np.ndarray):
        var = sigma ** 2
        if not sigma.shape == arr.shape[:-1]:
            e_s = "You provided a sigma array with a shape"
            e_s += "{0} for data with".format(sigma.shape)
            e_s += "shape {0}. Please provide a sigma array".format(arr.shape)
            e_s += " that matches the spatial dimensions of the data."
            raise ValueError(e_s)
    elif isinstance(sigma, (int, float)):
        var = sigma ** 2 * np.ones(arr.shape[:-1])

    dim = arr.shape[-1]
    if tau_factor is None:
        tau_factor = 1 + np.sqrt(dim / (patch_size ** 3))

    pool = mp.Pool(processes)

    #for k in range(patch_radius, arr.shape[2] - patch_radius):
    #    mykrange[i] = k
    #    i += 1

    theta = np.zeros(arr.shape, dtype=calc_dtype)
    thetax = np.zeros(arr.shape, dtype=calc_dtype)

    if calc_sigma:
        var = np.zeros(arr.shape[:-1], dtype=calc_dtype)
        thetavar = np.zeros(arr.shape[:-1], dtype=calc_dtype)

    allnum=(arr.shape[2] - 2 * patch_radius) * (arr.shape[1] - 2 * patch_radius) * (arr.shape[0] - 2 * patch_radius)

    duration1=time()
    if verbose:
        print("Begin mpca denoising")
    arr_shape = arr.shape
    for num_start in range(0,allnum,processes):

        #print(num_start)
        #kx1 = np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))
        #jx1 = np.floor((num - kx1 * ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius))
        #ix1 = num - kx1 * ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)) - jx1 * (arr_shape[0] - 2 * patch_radius)
        #kx2 = kx1 + 2 * patch_radius + 1
        #jx2 = jx1 + 2 * patch_radius + 1
        #ix2 = ix1 + 2 * patch_radius + 1


        if (num_start+processes)<allnum:
            minilist= list(range(num_start, num_start+processes))
        else:
            minilist= list(range(num_start, allnum))

        # Shorthand for indexing variables:

        #if not mask[ix1 + patch_radius, jx1 + patch_radius, kx1 + patch_radius]:

        #    X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
        #                    patch_size ** 3, dim)
        #minimask=mask[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
        #                    patch_size ** 3, dim)
        resultslist=[]
        resultslist = pool.starmap_async(pca_patchloop, [(patch_radius, arr[int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) % (arr_shape[0] - 2 * patch_radius))):int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))%(arr_shape[0] - 2 * patch_radius)) + 2 * patch_radius + 1),
                                                                        int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius))):int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius)) + 2*patch_radius + 1),
                                                                        int(np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))):int(np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))+2*patch_radius + 1)].reshape(patch_size ** 3, dim),
                                                          arr_shape,
                                                          mask[int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))%(arr_shape[0] - 2 * patch_radius)) + patch_radius),
                                                                        int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius)) + patch_radius),
                                                                        int(np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))+patch_radius)],
                                                          num, tau_factor, dim, is_svd, calc_sigma, verbose) for num in minilist]).get()

        #Xest=resultslist[0][0]
        #this_theta=resultslist[0][1]
        #this_var=resultslist[0][2]
        #resultslist[0][3]
        #Xest, this_theta, this_var, numlist = resultslist

        for i in range(len(minilist)):
            num = resultslist[i][3]
            kx1 = int(np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))))
            jx1 = int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (
                        arr_shape[0] - 2 * patch_radius)))
            ix1 = int(np.floor((num % ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) % (
                        arr_shape[0] - 2 * patch_radius)))
            kx2 = kx1 + 2 * patch_radius + 1
            jx2 = jx1 + 2 * patch_radius + 1
            ix2 = ix1 + 2 * patch_radius + 1

            theta[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[i][1]
            thetax[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[i][0] * resultslist[i][1]
            if calc_sigma:
                var[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[i][2] * resultslist[i][1]
                thetavar[ix1:ix2, jx1:jx2, kx1:kx2] += resultslist[i][1]
            else:
                var = 0
                thetavar = 0

        if verbose & ((((num_start/processes)+1) % 10000) == 0):
            print("finished " + str(int(np.floor(num_start/processes))+1) + " of " + str(int(np.floor(allnum/processes)+1)) )
            print("Process has been running for "+ str(time()-duration1) + "s")
            #pool.starmap_async(create_tracts, [(mypath, outpath, subject, step_size, function_processes,
            #                                                                        saved_streamlines, denoise, savefa, verbose) for subject in l]).get()
    #theta, thetax, var, thetavar = pca_patchloop(patch_radius, arr, k)


    #for k in range(patch_radius, arr.shape[2] - patch_radius):
    #    for j in range(patch_radius, arr.shape[1] - patch_radius):
    #        for i in range(patch_radius, arr.shape[0] - patch_radius):
    #mynumrange = range((arr.shape[2] - 2 * patch_radius) * (arr.shape[1] - 2 * patch_radius) *
    #                   (arr.shape[0] - 2 * patch_radius))
    #Xest, this_theta, this_var, ix1l, ix2l, jx1l, jx2l, kx1l, kx2l = pool.starmap_async(pca_patchloop_jusk, [(patch_radius,arr,num, is_svd) for num in mynumrange]).get()
    #print("hi")
    #print(mykrange)
    #arr=np.zeros(5)
    #thetap, thetaxp, varp, thetavarp = pool.starmap_async(pca_patchloop_justk, [(patch_radius, arr, k, is_svd,
    #                                                                             calc_sigma, calc_dtype, verbose) for k in mykrange]).get()
    if verbose:
        print("finished main computations, preparing matrix")


def pca_patchloop(patch_radius, X, arr_shape, mask, num, tau_factor, dim, is_svd, calc_sigma=True, verbose=False):

    # loop around and find the 3D patch for each direction at each pixel
    # for k in range(patch_radius, arr.shape[2] - patch_radius):

    #for k in range(patch_radius, arr.shape[2] - patch_radius):
    #    for j in range(patch_radius, arr.shape[1] - patch_radius):
    #        for i in range(patch_radius, arr.shape[0] - patch_radius):

    #arr_shape = np.shape(arr)
    #kx1 = np.floor(num / ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)))
    #jx1 = np.floor((num - k * ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius))) / (arr_shape[0] - 2 * patch_radius))
    #ix1 = num - k * ((arr_shape[1] - 2 * patch_radius) * (arr_shape[0] - 2 * patch_radius)) - j * (arr_shape[0] - 2 * patch_radius)
    #kx2 = kx1 + 2 * patch_radius + 1
    #jx2 = jx1 + 2 * patch_radius + 1
    #ix2 = ix1 + 2 * patch_radius + 1

    # Shorthand for indexing variables:
    #print("Start of process with num " + str(num))
    if mask:

        #X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
        #                patch_size ** 3, dim)
        # compute the mean and normalize
        M = np.mean(X, axis=0)
        # Upcast the dtype for precision in the SVD
        X = X - M
        patch_size = 2 * patch_radius + 1

        if is_svd:
            # PCA using an SVD
            U, S, Vt = svd(X, *svd_args)[:3]
            # Items in S are the eigenvalues, but in ascending order
            # We invert the order (=> descending), square and normalize
            # \lambda_i = s_i^2 / n
            d = S[::-1] ** 2 / X.shape[0]
            # Rows of Vt are eigenvectors, but also in ascending
            # eigenvalue order:
            W = Vt[::-1].T

        else:
            # PCA using an Eigenvalue decomposition
            C = np.transpose(X).dot(X)
            C = C / X.shape[0]
            [d, W] = eigh(C, turbo=True)

        if calc_sigma:
            # Random matrix theory
            this_var, ncomps = _pca_classifier(d, patch_size ** 3)
        else:
            # Predefined variance
            this_var = var[i, j, k]

        # Threshold by tau:
        tau = tau_factor ** 2 * this_var

        # Update ncomps according to tau_factor
        ncomps = np.sum(d < tau)
        W[:, :ncomps] = 0

        # This is equations 1 and 2 in Manjon 2013:
        Xest = X.dot(W).dot(W.T) + M
        Xest = Xest.reshape(patch_size,
                            patch_size,
                            patch_size, dim)
        # This is equation 3 in Manjon 2013:
        this_theta = 1.0 / (1.0 + dim - ncomps)

        #return theta,thetax,var,thetavar
    else:
        Xest = 0
        this_theta = 0
        this_var = 0
    return [Xest, this_theta, this_var, num]
    
    
def pca_patchloop_justk(patch_radius, arr, k, is_svd, calc_sigma, calc_dtype, verbose= False):

    # loop around and find the 3D patch for each direction at each pixel
    # for k in range(patch_radius, arr.shape[2] - patch_radius):
    theta = np.zeros(arr.shape, dtype=calc_dtype)
    thetax = np.zeros(arr.shape, dtype=calc_dtype)

    if verbose:
        print("beginning analysis of line " + str(k))

    if calc_sigma:
        var = np.zeros(arr.shape[:-1], dtype=calc_dtype)
        thetavar = np.zeros(arr.shape[:-1], dtype=calc_dtype)

    for j in range(patch_radius, arr.shape[1] - patch_radius):
        for i in range(patch_radius, arr.shape[0] - patch_radius):
            # Shorthand for indexing variables:
            if not mask[i, j, k]:
                continue

            ix1 = i - patch_radius
            ix2 = i + patch_radius + 1
            jx1 = j - patch_radius
            jx2 = j + patch_radius + 1
            kx1 = k - patch_radius
            kx2 = k + patch_radius + 1

            X = arr[ix1:ix2, jx1:jx2, kx1:kx2].reshape(
                            patch_size ** 3, dim)
            # compute the mean and normalize
            M = np.mean(X, axis=0)
            # Upcast the dtype for precision in the SVD
            X = X - M

            if is_svd:
                # PCA using an SVD
                U, S, Vt = svd(X, *svd_args)[:3]
                # Items in S are the eigenvalues, but in ascending order
                # We invert the order (=> descending), square and normalize
                # \lambda_i = s_i^2 / n
                d = S[::-1] ** 2 / X.shape[0]
                # Rows of Vt are eigenvectors, but also in ascending
                # eigenvalue order:
                W = Vt[::-1].T

            else:
                # PCA using an Eigenvalue decomposition
                C = np.transpose(X).dot(X)
                C = C / X.shape[0]
                [d, W] = eigh(C, turbo=True)

            if sigma is None:
                # Random matrix theory
                this_var, ncomps = _pca_classifier(d, patch_size ** 3)
            else:
                # Predefined variance
                this_var = var[i, j, k]

            # Threshold by tau:
            tau = tau_factor ** 2 * this_var

            # Update ncomps according to tau_factor
            ncomps = np.sum(d < tau)
            W[:, :ncomps] = 0

            # This is equations 1 and 2 in Manjon 2013:
            Xest = X.dot(W).dot(W.T) + M
            Xest = Xest.reshape(patch_size,
                                patch_size,
                                patch_size, dim)
            # This is equation 3 in Manjon 2013:
            this_theta = 1.0 / (1.0 + dim - ncomps)

            theta[ix1:ix2, jx1:jx2, kx1:kx2] += this_theta
            thetax[ix1:ix2, jx1:jx2, kx1:kx2] += Xest * this_theta
            if calc_sigma:
                var[ix1:ix2, jx1:jx2, kx1:kx2] += this_var * this_theta
                thetavar[ix1:ix2, jx1:jx2, kx1:kx2] += this_theta
            else:
                var = 0
                thetavar = 0

    if verbose:
        print("analysis of line " + str(k))
    return theta,thetax,var,thetavar
    #return Xest, this_theta, this_var, i, j, k
"""
