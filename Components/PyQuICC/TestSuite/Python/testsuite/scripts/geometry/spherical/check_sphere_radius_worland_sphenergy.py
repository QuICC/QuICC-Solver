"""Script to check Python full sphere worland of Chebyshev type operators 
"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import scipy.io as io

import quicc.geometry.worland.setup as wsetup
wsetup.type = "SphEnergy"
import quicc.geometry.spherical.sphere_radius_worland as sph

# Global resolution setup
keep_data = True
nrs = list([10, 15, 23, 32])
ls = list([0, 1, 2, 3, 16, 32])
ls_0 = [l for l in ls if l > 0]
allbcs = list([{0:0},{0:10},{0:11},{0:12},{0:13},{0:14},{0:20},{0:21}])
nobcs = list([{0:0}])
datadir = "data/geometry/spherical/sphere_worland_sphenergy/"
refdir = "ref/geometry/spherical/sphere_worland_sphenergy/"

def checkOperator(nrs, ls, bcs, fctOp, op, keep_data):
    has_ref = True
    print("-"*80)
    print(op + ":")
    for bc in bcs:
        for nr in nrs:
            for l in ls:
                A = fctOp(nr, l, bc)
                if keep_data:
                    fname = op + "_n" + str(nr) + "_l" + str(l) + "_bc" + str(bc[0]) + ".mtx"
                    io.mmwrite(datadir+fname, A)
                try:
                    R = io.mmread(refdir+fname)
                except FileNotFoundError:
                    has_ref = False
                    msg = "!?! Missing reference !?!"

                if has_ref:
                    # Make sure sparsity is identical
                    assert(A.nnz == R.nnz)
                    if A.nnz == 0:
                        err = 1e-16
                    else:
                        idx = np.argmax(np.abs(A.data-R.data)/(1.0 + np.abs(A.data)))
                        err = np.abs(A.data[idx]-R.data[idx])/(1.0 + np.abs(A.data[idx]))
                        relerr = np.abs(A.data[idx]-R.data[idx])/(np.abs(A.data[idx]))
                        err_val = A.data[idx]
                    if err > 1e-13:
                        msg = "!!! ERROR !!! max = {:.1e}".format(err) + ", value = {:.3e}".format(err_val)
                    elif err > 1e-14:
                        msg = "OK ({:.1e}".format(err) + ")"
                    else:
                        msg = "Perfect!"
                print("N = " + str(nr) + ", l = " + str(l) + ": " + msg)

# Generate zero block
checkOperator(nrs, ls, allbcs, sph.zblk, "zblk", keep_data)

# Generate r2 block
checkOperator(nrs, ls, nobcs, sph.r2, "r2", keep_data)

# Generate i1 block
checkOperator(nrs, ls, nobcs, sph.i1, "i1", keep_data)

# Generate i1qm block
checkOperator(nrs, ls_0, nobcs, sph.i1qm, "i1qm", keep_data)

# Generate i1qp block
checkOperator(nrs, ls, nobcs, sph.i1qp, "i1qp", keep_data)

# Generate i2 block
checkOperator(nrs, ls, nobcs, sph.i2, "i2", keep_data)

# Generate i2lapl block
checkOperator(nrs, ls, nobcs, sph.i2lapl, "i2lapl", keep_data)

# Generate i2qm block
checkOperator(nrs, ls_0, nobcs, sph.i2qm, "i2qm", keep_data)

# Generate i2qp block
checkOperator(nrs, ls, nobcs, sph.i2qp, "i2qp", keep_data)

# Generate i4 block
checkOperator(nrs, ls, nobcs, sph.i4, "i4", keep_data)

# Generate i4lapl block
checkOperator(nrs, ls, nobcs, sph.i4lapl, "i4lapl", keep_data)

# Generate i4lapl2 block
checkOperator(nrs, ls, nobcs, sph.i4lapl2, "i4lapl2", keep_data)

# Generate i4qm block
checkOperator(nrs, ls_0, nobcs, sph.i4qm, "i4qm", keep_data)

# Generate i4qp block
checkOperator(nrs, ls, nobcs, sph.i4qp, "i4qp", keep_data)

# Generate qid block
def qid_q1(nr, l, bc, coeff = 1.0):
    return sph.qid(nr, l, 1, bc, coeff)

def qid_q2(nr, l, bc, coeff = 1.0):
    return sph.qid(nr, l, 2, bc, coeff)

def qid_q4(nr, l, bc, coeff = 1.0):
    return sph.qid(nr, l, 4, bc, coeff)
checkOperator(nrs, ls, nobcs, qid_q1, "qid_q1", keep_data)
checkOperator(nrs, ls, nobcs, qid_q2, "qid_q2", keep_data)
checkOperator(nrs, ls, nobcs, qid_q4, "qid_q4", keep_data)

# Generate stencil block
def stencil_sq0(nr, l, bc):
    return sph.stencil(nr, l, bc, False)
def stencil_sq1(nr, l, bc):
    return sph.stencil(nr, l, bc, True)
#checkOperator(nrs, ls, nobcs, stencil_sq0, "stencil_sq0", keep_data)
checkOperator(nrs, ls, nobcs, stencil_sq1, "stencil_sq1", keep_data)
