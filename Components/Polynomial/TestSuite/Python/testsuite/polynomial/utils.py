"""Module provides generic functions for the testsuite."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

def read_complex(fname):
    """Read complex data in split storage"""

    tmp = np.genfromtxt(fname + "z")
    mid = tmp.shape[0]/2
    data = tmp[0:mid,:] + tmp[mid:,:]*1j
    return data

def read_real(fname):
    """Read real data"""

    return np.genfromtxt(fname)

def quadrature_error(data, ref, ignore = False):
    """Compute error between C++ data and reference"""

    err_g = np.max(np.abs(ref[0,:] - data[0,:]))
    err_w = np.max(np.abs(ref[1,:] - data[1,:]))
    if err_g < 1e-14 and err_w < 1e-14:
        code_msg = "\t\t" + "OK"
        code = 0
    else:
        code_msg = "\t\t" + str((err_g,err_w))
        code = 1
    if ignore:
        code_msg += " (ignored) "
    print(code_msg)
    return 0 if ignore else code

def test_summary(code):
    """Print summary of tests"""

    if code == 0:
        print("\t"+"-----------------------")
        print("\t"+"-- All tests passed! --")
        print("\t"+"-----------------------")
    else:
        print("\t"+"!!!!!!!!!!!!!!!!!!!!!!")
        print("\t"+"!! " + str(code) + " tests failed! !!")
        print("\t"+"!!!!!!!!!!!!!!!!!!!!!!")
