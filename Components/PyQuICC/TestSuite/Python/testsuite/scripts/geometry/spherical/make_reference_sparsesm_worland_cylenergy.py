"""Script to check Python full sphere worland of Chebyshev type operators 
"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import scipy.io as io

import quicc.geometry.worland.setup as wsetup
wsetup.type = "CylEnergy"
import quicc.geometry.spherical.sphere_radius_worland as sph

# Global resolution setup
nrs = list([16])
ls = list([0, 1, 2, 5, 20, 128])
ls_0 = [l for l in ls if l > 0]
datadir = "ref/SparseSM/Worland/CylEnergy/"

def writeOperator(nrs, ls, fctOp, op):
    print("-"*80)
    print(op + ":")
    for nr in nrs:
        for l in ls:
            A = fctOp(nr, l, {0:0})
            fname = op + "_l" + str(l) + "_r" + str(nr) + "_c" + str(nr) + "_t2" + ".mtx"
            io.mmwrite(datadir+fname, A)

# Generate r2 block
writeOperator(nrs, ls, sph.r2, "R2")

# Generate i2 block
writeOperator(nrs, ls, sph.i2, "I2")

# Generate i4 block
writeOperator(nrs, ls, sph.i4, "I4")

# Generate i4d1r1 block
#writeOperator(nrs, ls, sph.i4d1r1, "I4D1R1")

# Generate i6 block
#writeOperator(nrs, ls, sph.i6, "I6")

# Generate i6cyllaplh block
#writeOperator(nrs, ls, sph.i6cyllaplh, "I6CylLaplh")
