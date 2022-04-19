#! /usr/bin/env python

import numpy as np
import testsuite.polynomial.utils as utils

base_dir = 'data/Polynomial/Quadrature/'

def check_worlandchebyshev(n):
    """Check D^1 projector"""

    print("\t" + "Checking Worland-Chebyshev rule...")
    data = utils.read_real(base_dir + "worlandchebyshev.dat")
    import numpy.polynomial.chebyshev as cheby
    ref = np.vstack(cheby.chebgauss(n))[:,::-1]
    # Map to x = 2r^2-1
    ref[0,:] = ((ref[0,:] + 1.0)/2.0)**0.5
    # Map for integral r=[0,1]
    ref[1,:] *= 0.5
    return utils.quadrature_error(data, ref)

print("Quadrature rules")
code = 0
code += check_worlandchebyshev(256)

utils.test_summary(code)

import sys
sys.exit(code)
