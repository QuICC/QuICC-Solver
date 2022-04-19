"""Generate test problem matrices for sparse solver tests"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.io as io
import scipy.sparse as spsp
if True:
    import matplotlib.pylab as pl
    has_error_plot = True
else:
    has_error_plot = False

import quicc.geometry.cylindrical.cylinder_radius_worland as cyl
import quicc.transform.cylinder_worland as transf
import quicc.base.utils as utils

save_path = ""

def dclaplh(n, m, restriction = None):
    """Accuracy test for i2lapl operator"""

    B = cyl.i2laplh(n, m, cyl.radbc.no_bc()).tocsr()
    A = cyl.i2(n, m, {0:10}).tocsr()
    print(B.todense())

    return (A, B)

if __name__ == "__main__":
    # Set test parameters
    print("Doing nothing!")
