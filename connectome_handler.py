import numpy as np
from numpy import ravel_multi_index
#from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)
from collections import defaultdict, OrderedDict
from dipy.tracking.streamline import Streamlines
import warnings

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


def _to_voxel_coordinates_warning(streamline, lin_T, offset):
    """Applies a mapping from streamline coordinates to voxel_coordinates,
    raises an error for negative voxel values."""
    inds = np.dot(streamline, lin_T)
    inds += offset
    if inds.min().round(decimals=6) < 0:
        warnings.warn('streamline has points that map to negative voxel'
                         ' indices')
    return inds.astype(int)


def convert_label_3dto4d(label_volume):
    value_set = np.unique(label_volume)
    value_set.remove(0)
    value_set = np.delete(value_set,0)
    num_vals = len(value_set)
    label_shape = np.shape(label_volume)
    label_volume_new = np.zeros([label_shape,num_vals])
    for i in label_shape[0]:
        for j in label_shape[1]:
            for k in label_shape[2]:
                if label_volume[i,j,k] != 0:
                    l = value_set.index(label_volume[i,j,k])
                    label_volume_new[i,j,k,l] = label_volume[i,j,k]
    return label_volume_new


def dilate_label_volume(label_volume_4d, dilation = 1, volumes_todilate = None):
    if volumes_todilate is None:
        #volumes_todilate = np.shape(label_volume_4d)[3]
        volumes_todilate = np.unique(label_volume_4d).remove(0)


def label_weights_matrix(label_volume):

    value_set = np.unique(label_volume)
    value_set = np.delete(value_set,0)
    size_matrix = len(value_set)
    matrix = np.zeros([size_matrix,size_matrix])
    volume_sizes = np.zeros(len(value_set))

    for i, value in enumerate(value_set):
        volume_sizes[i] = np.sum(label_volume == value)
    for i in np.arange(np.shape(volume_sizes)[0]):
        for j in np.arange(np.shape(volume_sizes)[0]):
            matrix[i,j] = (volume_sizes[i] + volume_sizes[j])/2

    return(matrix)



def connectivity_matrix_test(streamlines, affine, label_volume, inclusive=False,
                        symmetric=True, return_mapping=False,
                        mapping_as_streamlines=False):
    print('success')
    print(np.shape(streamlines))

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

    x = x.astype(int)
    x = ravel_multi_index(x, shape)
    out = np.bincount(x, weights, minlength=np.prod(shape))
    out.shape = shape

    return out


def connectivity_matrix_func(pruned_streamlines_SL, affine_streams, labelmask, inclusive=False, symmetric = True,return_mapping = False,mapping_as_streamlines = False, reference_weighting = None, volume_weighting = False, function_processes = 1, verbose=False):

    n = function_processes
    size_SL = np.size(pruned_streamlines_SL)
    listcut = []
    listcut.append(0)
    for i in np.arange(n - 1):
        listcut.append(np.int(((i + 1) * size_SL) / n))
        if verbose:
            print(size_SL, i + 1, n)
    listcut.append(size_SL)
    if verbose:
        print(listcut)
    pruned_cut = []
    for i in np.arange(n):
        pruned_cut.append(pruned_streamlines_SL[listcut[i]:listcut[i+1]])
    pool = Pool()
    if verbose:
        print("The streamline is split into "+str(function_processes)+" of size "+str(np.int(size_SL / n)))

    return_mapping = True
    connectomic_results = pool.starmap_async(connectivity_matrix_custom, [(Streamlines(pruned_streamlines_SL[listcut[i]:listcut[i+1]]), affine_streams, labelmask,
                                                                          inclusive, symmetric, return_mapping,
                                                                          mapping_as_streamlines,reference_weighting, volume_weighting) for i in np.arange(n)]).get()

    M = np.zeros(np.shape(connectomic_results[0][0]))
    grouping = {}
    i=0
    for connectome_results in connectomic_results:
        M += connectome_results[0]
        for key, val in connectome_results[1].items():
            if key in grouping:
                grouping[key].extend([j+listcut[i] for j in val])
            else:
                grouping[key] = val
        i = i + 1

    return M, grouping

def retweak_points(points, shape):
    pointsT = points.T
    for axis in [0,1,2]:
        if np.min(pointsT[axis])<0:
            print(f'There are {np.sum(pointsT[axis]<0)} points that are negative in axis {axis}')
            pointsT[axis][pointsT[axis]<0] = 0
        if np.max(pointsT[axis]>=shape[axis]):
            print(f'There are {np.sum(pointsT[axis]>=shape[axis])} points that are above maximum in axis {axis}')
            pointsT[axis][pointsT[axis]>0] = shape[axis] - 1
    pointsnew = pointsT.T

    return pointsnew

def connectivity_matrix_custom(streamlines, affine, label_volume,
                        inclusive=False, symmetric=True, return_mapping=False,
                        mapping_as_streamlines=False,  reference_weighting = None, volume_weighting=False):

    # Error checking on label_volume
    from itertools import combinations, groupby

    if volume_weighting:
        volume_weights = label_weights_matrix(label_volume)

    kind = label_volume.dtype.kind
    labels_positive = ((kind == 'u') or
                       ((kind == 'i') and (label_volume.min() >= 0)))
    valid_label_volume = (labels_positive and label_volume.ndim == 3)
    if not valid_label_volume:
        raise ValueError("label_volume must be a 3d integer array with"
                         "non-negative label values")

    # If streamlines is an iterator
    if return_mapping and mapping_as_streamlines:
        streamlines = list(streamlines)

    if inclusive:
        # Create ndarray to store streamline connections
        edges = np.ndarray(shape=(4, 0), dtype=float)
        lin_T, offset = _mapping_to_voxel(affine)
        for sl, _ in enumerate(streamlines):
            # Convert streamline to voxel coordinates

            entire = _to_voxel_coordinates_warning(streamlines[sl], lin_T, offset)
            entire = retweak_points(entire, np.shape(label_volume))

            i, j, k = entire.T

            ## streamlinep3794 endpoints for 02363 is [[ 99.23122215  97.81464577  44.45293291] [167.02122498  74.64361191  46.30043031]]

            if reference_weighting is not None:
                ref_values = reference_weighting[i, j, k]
                ref_mean = np.mean(ref_values)
            else:
                ref_mean = 1
            if symmetric:
                # Create list of all labels streamline passes through
                entirelabels = list(OrderedDict.fromkeys(label_volume[i, j, k]))
                # Append all connection combinations with streamline number

                for comb in combinations(entirelabels, 2):
                    edges = np.append(edges, [[comb[0]], [comb[1]], [sl], [ref_mean]],
                                      axis=1)
            else:
                # Create list of all labels streamline passes through, keeping
                # order and whether a label was entered multiple times
                entirelabels = list(groupby(label_volume[i, j, k]))
                # Append connection combinations along with streamline number,
                # removing duplicates and connections from a label to itself
                combs = set(combinations([z[0] for z in entirelabels], 2))
                for comb in combs:
                    if comb[0] == comb[1]:
                        pass
                    else:
                        edges = np.append(edges, [[comb[0]], [comb[1]], [sl], [ref_mean]],
                                          axis=1)
        if symmetric:
            edges[0:2].sort(0)
        mx = label_volume.max() + 1
        matrix = ndbincount(edges[0:2], weights=None, shape=(mx, mx))
        if reference_weighting is not None:
            matrix_refweighted = ndbincount(edges[0:2], weights = edges[3], shape=(mx, mx))
        else:
            matrix_refweighted = np.copy(matrix)

        matrix_refweighted = matrix_refweighted.astype('float64')

        if symmetric:
            matrix = np.maximum(matrix, matrix.T)
            matrix_refweighted = np.maximum(matrix_refweighted,matrix_refweighted.T)

        matrix_vol = np.copy(matrix)
        matrix_vol = matrix_vol.astype('float64')
        matrix_vol_refweighted = np.copy(matrix_refweighted)

        if volume_weighting:
            matrix_vol[1:,1:] = matrix_vol[1:,1:] / volume_weights
            if reference_weighting is not None:
                matrix_vol_refweighted[1:, 1:] = matrix_vol_refweighted[1:, 1:] / volume_weights
            else:
                matrix_vol_refweighted = matrix_vol

        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b, c, d) in enumerate(edges.T):
                mapping[int(a), int(b)].append(int(c))
            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            return matrix, matrix_vol, matrix_refweighted, matrix_vol_refweighted, mapping

        return matrix, matrix_vol, matrix_refweighted, matrix_vol_refweighted
    else:
        # take the first and last point of each streamline
        endpoints = np.array([sl[0::len(sl)-1] for sl in streamlines])

        # Map the streamlines coordinates to voxel coordinates
        lin_T, offset = _mapping_to_voxel(affine)

        endpoints = _to_voxel_coordinates_warning(endpoints, lin_T, offset)
        endpoints = retweak_points(endpoints, np.shape(label_volume))

        # get labels for label_volume

        i, j, k = endpoints.T
        endlabels = label_volume[i, j, k]  ## CREATE CONSTRAIN FUNCTION SO WE STOP HAVING THE out of bounds error!!!!
        if symmetric:
            endlabels.sort(0)
        mx = label_volume.max() + 1

        if reference_weighting is not None:
            weights = np.zeros(np.size(streamlines),dtype = float)
            l=0
            for streamline in streamlines:
                coordinates = _to_voxel_coordinates(streamline, lin_T, offset)
                i, j, k = coordinates.T
                ref_values = reference_weighting[i, j, k]
                ref_mean = np.mean(ref_values)
                weights[l] = ref_mean
                l+=1

        matrix = ndbincount(endlabels, shape=(mx, mx))
        if symmetric:
            matrix = np.maximum(matrix, matrix.T)

        if reference_weighting is not None:
            matrix_refweighted = ndbincount(endlabels, weights = weights, shape=(mx, mx))
            if symmetric:
                matrix_refweighted = np.maximum(matrix_refweighted, matrix_refweighted.T)
        else:
            matrix_refweighted = np.copy(matrix)

        matrix_refweighted = matrix_refweighted.astype('float64')
        matrix_vol = np.copy(matrix)
        matrix_vol = matrix_vol.astype('float64')
        matrix_vol_refweighted = np.copy(matrix_refweighted)

        if volume_weighting:
            matrix_vol[1:,1:] = matrix_vol[1:,1:] / volume_weights
            if reference_weighting is not None:
                matrix_vol_refweighted[1:, 1:] = matrix_vol_refweighted[1:, 1:] / volume_weights
            else:
                matrix_vol_refweighted = matrix_vol

        if return_mapping:
            mapping = defaultdict(list)
            for i, (a, b) in enumerate(endlabels.T):
                mapping[int(a), int(b)].append(i)

            # Replace each list of indices with the streamlines they index
            if mapping_as_streamlines:
                for key in mapping:
                    mapping[key] = [streamlines[i] for i in mapping[key]]

            # Return the mapping matrix and the mapping
            return matrix, matrix_vol, matrix_refweighted, matrix_vol_refweighted, mapping

        return matrix, matrix_vol, matrix_refweighted, matrix_vol_refweighted