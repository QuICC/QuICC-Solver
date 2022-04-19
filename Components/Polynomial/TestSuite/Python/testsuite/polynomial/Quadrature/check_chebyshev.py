#! /usr/bin/env python

import numpy as np
import testsuite.polynomial.utils as utils

base_dir = 'data/Polynomial/Quadrature/'

def check_chebyshev(n):
    """Check Chebyshev rule"""

    print("\t" + "Checking Chebyshev rule...")
    data = utils.read_real(base_dir + "chebyshev.dat")
    import numpy.polynomial.chebyshev as cheby
    ref = np.vstack(cheby.chebgauss(n))
    return utils.quadrature_error(data, ref)

print("Quadrature rules")
code = 0
code += check_chebyshev(256)

utils.test_summary(code)

import sys
sys.exit(code)
