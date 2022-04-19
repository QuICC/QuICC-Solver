#! /usr/bin/env python

import numpy as np
import testsuite.polynomial.utils as utils

base_dir = 'data/Polynomial/Quadrature/'

def check_legendre(n):
    """Check Legendre rule"""

    print("\t" + "Checking Legendre rule...")
    data = utils.read_real(base_dir + "legendre.dat")
    import numpy.polynomial.legendre as leg
    ref = np.vstack(leg.leggauss(n))
    return utils.quadrature_error(data, ref)

print("Quadrature rules")
code = 0
code += check_legendre(256)

utils.test_summary(code)

import sys
sys.exit(code)
