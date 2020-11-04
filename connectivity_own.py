from __future__ import division, print_function, absolute_import

import numpy as np
from collections import defaultdict
from tract_handler import prune_streamlines
from dipy.tracking.streamline import Streamlines
"""This is a helper module for dipy.tracking.utils"""


def _mapping_to_voxel(affine):
    """Inverts affine and returns a mapping so voxel coordinates. This
    function is an implementation detail and only meant to be used with
    ``_to_voxel_coordinates``.

    Parameters
    ----------
    affine : array_like (4, 4)
        The mapping from voxel indices, [i, j, k], to real world coordinates.
        The inverse of this mapping is used unless `affine` is None.

    Returns
    --------
    lin_T : array (3, 3)
        Transpose of the linear part of the mapping to voxel space, (ie
        ``inv(affine)[:3, :3].T``)
    offset : array or scalar
        Offset part of the mapping (ie, ``inv(affine)[:3, 3]``) + ``.5``. The
        half voxel shift is so that truncating the result of this mapping
        will give the correct integer voxel coordinate.

    Raises
    ------
    ValueError
        If both affine and voxel_size are None.

    """
    if affine is None:
        raise ValueError("no affine specified")

    affine = np.array(affine, dtype=float)
    inv_affine = np.linalg.inv(affine)
    lin_T = inv_affine[:3, :3].T.copy()
    offset = inv_affine[:3, 3] + .5

    return lin_T, offset


def _to_voxel_coordinates(streamline, lin_T, offset):
    """Applies a mapping from streamline coordinates to voxel_coordinates,
    raises an error for negative voxel values."""
    inds = np.dot(streamline, lin_T)
    inds += offset
    if inds.min().round(decimals=6) < 0:
        raise IndexError('streamline has points that map to negative voxel'
                         ' indices')
    return inds.astype(int)


def ndbincount(x, weights=None, shape=None):
    """Like bincount, but for nd-indices.

    Parameters
    ----------
    x : array_like (N, M)
        M indices to a an Nd-array
    weights : array_like (M,), optional
        Weights associated with indices
    shape : optional
        the shape of the output
    """
    x = np.asarray(x)
    if shape is None:
        shape = x.max(1) + 1

    x = np.ravel_multi_index(x, shape)
    out = np.bincount(x, weights, minlength=np.prod(shape))
    out.shape = shape

    return out


def connectivity_matrix_special(streamlines, affine, label_volume, symmetric=True, return_mapping=False,
                                mapping_as_streamlines=False, cutoff=4):
    """Counts the streamlines that start and end at each label pair.

    Parameters
    ----------
    streamlines : sequence
        A sequence of streamlines.
    affine : array_like (4, 4)
        The mapping from voxel coordinates to streamline coordinates.
        The voxel_to_rasmm matrix, typically from a NIFTI file.
    label_volume : ndarray
        An image volume with an integer data type, where the intensities in the
        volume map to anatomical structures.
    symmetric : bool, True by default
        Symmetric means we don't distinguish between start and end points. If
        symmetric is True, ``matrix[i, j] == matrix[j, i]``.
    return_mapping : bool, False by default
        If True, a mapping is returned which maps matrix indices to
        streamlines.
    mapping_as_streamlines : bool, False by default
        If True voxel indices map to lists of streamline objects. Otherwise
        voxel indices map to lists of integers.

    Returns
    -------
    matrix : ndarray
        The number of connection between each pair of regions in
        `label_volume`.
    mapping : defaultdict(list)
        ``mapping[i, j]`` returns all the streamlines that connect region `i`
        to region `j`. If `symmetric` is True mapping will only have one key
        for each start end pair such that if ``i < j`` mapping will have key
        ``(i, j)`` but not key ``(j, i)``.

    """
    # Error checking on label_volume
    kind = label_volume.dtype.kind
    labels_positive = ((kind == 'u') or
                       ((kind == 'i') and (label_volume.min() >= 0)))
    valid_label_volume = (labels_positive and label_volume.ndim == 3)
    if not valid_label_volume:
        raise ValueError("label_volume must be a 3d integer array with"
                         "non-negative label values")

    # If streamlines is an iterators
    if return_mapping and mapping_as_streamlines:
        streamlines_orig = streamlines
        streamlines = list(streamlines)
    # take the first and last point of each streamline
    endpoints = [sl[0::len(sl)-1] for sl in streamlines]

    # Map the streamlines coordinates to voxel coordinates
    lin_T, offset = _mapping_to_voxel(affine)
    endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)

    # get labels for label_volume
    i, j, k = endpoints.T
    try:
        endlabels = label_volume[i, j, k]
    except IndexError:
        pruned_streamlines = prune_streamlines(list(streamlines_orig), label_volume, cutoff=cutoff, verbose=False)
        streamlines = Streamlines(pruned_streamlines)
        if return_mapping and mapping_as_streamlines:
            streamlines = list(streamlines)
        endpoints = [sl[0::len(sl) - 1] for sl in streamlines]
        endpoints = _to_voxel_coordinates(endpoints, lin_T, offset)
        i, j, k = endpoints.T
        endlabels = label_volume[i, j, k]

    mx = label_volume.max() + 1
    matrix = ndbincount(endlabels, shape=(mx, mx))
    if symmetric:
        matrix = np.maximum(matrix, matrix.T)

    if return_mapping:
        mapping = defaultdict(list)
        for i, (a, b) in enumerate(endlabels.T):
            mapping[a, b].append(i)

        # Replace each list of indices with the streamlines they index
        if mapping_as_streamlines:
            for key in mapping:
                mapping[key] = [streamlines[i] for i in mapping[key]]

        # Return the mapping matrix and the mapping
        return matrix, mapping
    else:
        return matrix