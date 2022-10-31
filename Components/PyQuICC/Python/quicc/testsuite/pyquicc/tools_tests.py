"""Module provides tests for PyQuICC tools."""

import numpy as np


def check_tuple(t):
    return (t == (1,2,3,5))

def arr():
    return np.array([0.1,0.2,0.3,0.5])

def check_arr2list(a):
    result = isinstance(a, list)
    ref = [0.1,0.2,0.3,0.5]
    result = result and isinstance(ref, list)
    result = result and (a == ref)
    if(not result):
        print(a)
        print(ref)
    return result

def mat():
    return np.array([[0.1,0.2,0.3,0.5],[1.1,1.2,1.3,1.5],[2.1,2.2,2.3,2.5]])

def check_mat2list(a):
    result = isinstance(a, list)
    ref = [0.1,1.1,2.1,0.2,1.2,2.2,0.3,1.3,2.3,0.5,1.5,2.5]
    result = result and isinstance(ref, list)
    result = result and (a == ref)
    if(not result):
        print(a)
        print(ref)
    return result

def spMat():
    import scipy.sparse as spsp

    d0 = np.arange(1,6)
    d1 = -np.arange(1,5)
    mat = spsp.diags([d0,d1], [0,1]).tocoo()
    return mat

def spMat2():
    mat = spMat().tocsr()
    return mat

def check_sparse2triplets(s):
    result = isinstance(s, list)
    import quicc.base.utils as utils
    mat = spMat()
    ref = utils.triplets(mat)
    result = result and (sorted(s) == sorted(ref))
    if(not result):
        print(sorted(s))
        print(sorted(ref))
    return result

def spMatZ():
    import scipy.sparse as spsp

    d0 = np.arange(1,6) + (np.arange(1,6) + 10)*1j
    d1 = -np.arange(1,5) - (np.arange(1,5) + 6)*1j
    mat = spsp.diags([d0,d1], [0,1]).tocoo()
    return mat
