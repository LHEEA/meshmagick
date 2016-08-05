#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

def merge_duplicates_rows(arr, decimals=8, return_index=False):
    # The technique relies on np.unique and has been proposed in the following post:
    # https://stackoverflow.com/questions/17273022/python-numpy-build-2d-array-without-adding-duplicate-rows-for-triangular-mesh

    # Rounding array to the specified number of decimals for fair comparison in np.unique
    rounded_arr = np.round(arr, decimals=decimals)

    # Defining a row_dtype so that the rows are compared as a whole
    row_dtype = np.dtype((np.void, (3 * rounded_arr.dtype.itemsize)))
    _, index, inv = np.unique(rounded_arr.view(row_dtype), return_index=True, return_inverse=True)

    # Re-introducing initial unrounded values
    uniq = arr[index]

    print "%u rows merged" % (arr.shape[0] - uniq.shape[0])

    if return_index:
        return uniq, inv
    else:
        return uniq
