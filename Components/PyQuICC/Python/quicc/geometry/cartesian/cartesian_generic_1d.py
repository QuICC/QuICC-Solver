"""Module provides functions to generate generic sparse operators for a single cartesian direction"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.linalg as splin
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.transform.cartesian as transf
import quicc.geometry.cartesian.cartesian_boundary_1d as c1dbc


def generic(nx, xi, xo, q, expr, var, bc, coeff = 1.0, ntrunc = -1):
    """Compute the spectral operator for a generic expression"""

    # Convert expression to function
    func = sy.utilities.lambdify(var, expr.simplify(), 'numpy')

    # Convert to physical space values
    grid = transf.grid(2*nx)
    phys = func(grid)
    spec = transf.tocheb(phys)[0:nx]
    spec *= (np.abs(spec)/np.max(np.abs(spec)) > np.spacing(1))
    #spec *= (np.abs(spec) > np.spacing(1))
    spec[1:] *= 2

    # Truncate non constant coefficient expansion
    if ntrunc > -1 and ntrunc+1 < nx:
        spec[ntrunc+1:] = 0

    # Form operator matrix
    h = splin.hankel(spec)
    h[0,:] = 0
    mat = splin.toeplitz(np.append(2*spec[0], spec[1:]))
    mat = 0.5*(mat + h)
    mat[0,:] *= 2.0
    mat[:,0] *= 0.5

    # Clear q top rows
    mat[0:q,:] = 0
    mat = coeff*spsp.coo_matrix(mat)
    return c1dbc.constrain(mat, xi, xo, bc)

def mult_generic(op, nx, xi, xo, q, expr, var, bc, ntrunc = -1):
    """Compute the spectral operator for a generic expression"""

    gen = generic(2*nx, xi, xo, q, expr, var, c1dbc.no_bc(), ntrunc = nx)
    mat = op(2*nx)*gen
    mat = mat[0:nx,0:nx]
    return c1dbc.constrain(mat, xi, xo, bc)
