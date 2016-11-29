#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

# def merge_duplicates_rows(arr, decimals=20, return_index=False):
#     # TODO: ecrire docstring
#     # FIXME: fonction non robuste !!!
#     # The technique relies on np.unique and has been proposed in the following post:
#     # https://stackoverflow.com/questions/17273022/python-numpy-build-2d-array-without-adding-duplicate-rows-for-triangular-mesh
#
#     # Rounding array to the specified number of decimals for fair comparison in np.unique
#     rounded_arr = np.round(arr, decimals=decimals)
#
#     # Defining a row_dtype so that the rows are compared as a whole
#     row_dtype = np.dtype((np.void, (3 * rounded_arr.dtype.itemsize)))
#     _, index, inv = np.unique(rounded_arr.view(row_dtype), return_index=True, return_inverse=True)
#
#     # Re-introducing initial unrounded values
#     uniq = arr[index]
#
#     print "%u rows merged" % (arr.shape[0] - uniq.shape[0])
#
#     if return_index:
#         return uniq, inv
#     else:
#         return uniq


def merge_duplicate_rows(arr, atol=1e-8, return_index=False):
    """merge_duplicates(_vertices, _faces, verbose=False, atol=1e-8)

    Returns a new node array where close nodes have been merged into one node (following atol). It also returns
    the connectivity array _faces with the new node IDs.

    Parameters:
        arr: ndarray
            numpy array of the coordinates of the mesh's nodes
        verbose[optional]: bool
            if set to True, displays information on the merge procedure
        atol[optional]: float
            the tolerance used to define nodes that are coincident and
            that have to be merged

    Returns:
        _vertices: ndarray
            numpy array of the coordinates of the mesh's nodes where
            every node is different
        _faces: ndarray
            numpy array of the _faces' nodes connectivities, accordingly
            to the new node list that has been merged
    """
    # TODO: Refaire la documentation --> les entrees sorties ont change !!

    # TODO : Set a tolerance option in command line arguments

    # This function is a bottleneck in the clipping routines
    # TODO: use np.unique to cluster groups --> acceleration !!

    # atol = pow(10, -decimals)

    nv, nbdim = arr.shape

    levels = [0, nv]
    iperm = np.arange(nv)

    for dim in range(nbdim):
        # Sorting the first dimension
        values = arr[:, dim].copy()
        if dim > 0:
            values = values[iperm]
        levels_tmp = []
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            if istop-istart > 1:
                level_values = values[istart:istop]
                iperm_view = iperm[istart:istop]

                iperm_tmp = level_values.argsort()

                level_values[:] = level_values[iperm_tmp]
                iperm_view[:] = iperm_view[iperm_tmp]

                levels_tmp.append(istart)
                vref = values[istart]

                for idx in xrange(istart, istop):
                    cur_val = values[idx]
                    if np.abs(cur_val - vref) > atol:
                        levels_tmp.append(idx)
                        vref = cur_val

            else:
                levels_tmp.append(levels[ilevel])
        if len(levels_tmp) == nv:
            # No duplicate rows
            # if verbose:
            # print "\t -> No duplicate _vertices detected :)"
            if return_index:
                newID = np.arange(nv)
            break

        levels_tmp.append(nv)
        levels = levels_tmp

    else:
        # Building the new merged node list
        arr_tmp = []
        newID = np.arange(nv)
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            arr_tmp.append(arr[iperm[istart]])
            newID[iperm[range(istart, istop)]] = ilevel
        arr = np.array(arr_tmp, dtype=float)
        # Applying renumbering to cells
        # if F is not None:
        #     for cell in F:
        #         cell[:] = newID[cell]

        # if verbose:
        # nv_new = arr.shape[0]
        # print "\t -> Initial number of nodes : {:d}".format(nv)
        # print "\t -> New number of nodes     : {:d}".format(nv_new)
        # print "\t -> {:d} nodes have been merged".format(nv-nv_new)

    # if F is not None:
    #     if return_index:
    #         return arr, F, newID
    #     else:
    #         return arr, F
    # else:
    if return_index:
        return arr, newID
    else:
        return arr











