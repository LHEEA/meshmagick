#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

def merge_duplicate_rows(arr, atol=1e-8, return_index=False):
    """Returns a new node array where close nodes have been merged into one node (following atol).

    Parameters
    ----------
    arr : array_like
        array of the coordinates of the mesh's nodes
    atol[optional] : float
        the tolerance used to define nodes that are coincident and
        that have to be merged
    return_index : bool
        If true, it also returns the array for new indices of vertices

    Returns
    -------
    arr : ndarray
        array of the coordinates of the mesh's nodes where
        every node is different
    newID : ndarray, optional
        array of the new new vertices IDs
    """
    # TODO: Refaire la documentation --> les entrees sorties ont change !!

    # TODO : Set a tolerance option in command line arguments

    # This function is a bottleneck in the clipping routines
    # TODO: use np.unique to cluster groups --> acceleration !!

    # atol = pow(10, -decimals)
    
    arr = np.asarray(arr)
    
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











