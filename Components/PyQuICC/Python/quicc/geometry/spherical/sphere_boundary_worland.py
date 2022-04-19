"""Module provides functions to generate the boundary conditions in a sphere with Worland expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.geometry.spherical.sphere_radius_boundary_worland as radbc
from quicc.geometry.spherical.sphere_radius_boundary_worland import no_bc

def standard_truncation(nr, l):
    """Impose standard truncation in radial direction"""

    return nr

def triangular_truncation(nr, l):
    """Impose triangular truncation in radial direction"""

    return nr - l

def radial_truncation(nr, l):
    """Generic radial truncation modification"""

    #return triangular_truncation(nr,l)
    return standard_truncation(nr,l)

def no_bc():
    """Get a no boundary condition flag"""

    return radbc.no_bc()

def constrain(mat, nr, maxnl, m, bc, zero_l_zero = False, restriction = None):
    """Contrain the matrix with the tau boundary condition"""

    lnr = radial_truncation(nr, m)
    if bc[0] > 0:
        if m == 0 and zero_l_zero:
            bc_mat = spsp.coo_matrix((lnr,lnr))
        else:
            if restriction is None or m in restriction:
                bcMat = spsp.coo_matrix((lnr,lnr))
                bc_mat = radbc.constrain(bcMat, m, bc)
            else:
                bc_mat = spsp.coo_matrix((lnr,lnr))
        for l in range(m+1, maxnl):
            lnr = radial_truncation(nr, l)
            if restriction is None or l in restriction:
                bcMat = spsp.coo_matrix((lnr,lnr))
                bcMat = radbc.constrain(bcMat, l, bc)
            else:
                bcMat = spsp.coo_matrix((lnr,lnr))
            bc_mat = spsp.block_diag((bc_mat,bcMat), format = 'coo')

        bc_mat = mat + bc_mat

        return bc_mat
    else:
        return mat
