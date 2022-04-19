"""Module provides functions to generate the boundary conditions in a cartesian 2D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.polynomial.chebyshev as cheby
import scipy.sparse as spsp

import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import quicc.transform.cartesian as phys
import quicc.base.utils as utils


def no_bc():
    """Get a no boundary condition flag"""

    return {'x':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

def bid(nx, xi, xo, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(nx-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = c1dbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return c1dbc.constrain(mat, xi, xo, tbc)

def constrain(mat, nx, xi, xo, nz, zi, zo, qx, qz, bc, location = 't', restriction = None):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'x')
    sx = 0
    dx = 0
    sz = 0
    dz = 0
    if priority == 'x':
        sx = qx
    elif priority == 'z':
        sz = qz
    elif priority == 'n':
        sx = qx
        sz = qz
    elif priority == 'sx':
        sx = qx
        dx = -qx
    elif priority == 'sz':
        sz = qz
        dz = -qz
    else:
        raise RuntimeError("Unknown boundary condition priority!")

    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, xi, xo, bc['x'], location = location)
        if bc['x'].get('kron',0) == 0 or bc['x']['kron'] == "id":
            bc_mat = bc_mat + utils.restricted_kron_2d(bid(nz, zi, zo, sz, dz, bc['z'], location = location), bcMat, restriction = restriction)
        elif bc['x']['kron'] == "d1":
            bc_mat = bc_mat + utils.restricted_kron_2d(bid(nz, zi, zo, sz, dz, bc['z'], location = location)*c1d.d1(nz, zi, zo, c1dbc.no_bc()), bcMat, restriction = restriction)
        elif bc['x']['kron'] == "d2":
            bc_mat = bc_mat + utils.restricted_kron_2d(bid(nz, zi, zo, sz, dz, bc['z'], location = location)*c1d.d2(nz, zi, zo, c1dbc.no_bc()), bcMat, restriction = restriction)
        elif bc['x']['kron'] == "i1":
            bc_mat = bc_mat + utils.restricted_kron_2d(bid(nz, zi, zo, sz, dz, bc['z'], location = location)*c1d.i1(nz, zi, zo, c1dbc.no_bc()), bcMat, restriction = restriction)
        elif bc['x']['kron'] == "q1":
            bc_mat = bc_mat + utils.restricted_kron_2d(bid(nz, zi, zo, sz, dz, bc['z'], location = location)*c1d.qid(nz, zi, zo, 1, c1dbc.no_bc()), bcMat, restriction = restriction)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, zi, zo, bc['z'], location = location)
        if bc['x'][0] >= 0:
            if bc['z'].get('kron',0) == 0 or bc['z']['kron'] == "id":
                bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, bid(nx, xi, xo, sx, dx, bc['x'], location = location), restriction = restriction)
            elif bc['z']['kron'] == "d1":
                bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, bid(nx, xi, xo, sx, dx, bc['x'], location = location)*c1d.d1(nx, xi, xo, c1dbc.no_bc()), restriction = restriction)
            elif bc['z']['kron'] == "d2":
                bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, bid(nx, xi, xo, sx, dx, bc['x'], location = location)*c1d.d2(nx, xi, xo, c1dbc.no_bc()), restriction = restriction)
            elif bc['z']['kron'] == "i1":
                bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, bid(nx, xi, xo, sx, dx, bc['x'], location = location)*c1d.i1(nx, xi, xo, c1dbc.no_bc()), restriction = restriction)
            elif bc['z']['kron'] == "q1":
                bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, bid(nx, xi, xo, sx, dx, bc['x'], location = location)*c1d.qid(nx, xi, xo, 1, c1dbc.no_bc()), restriction = restriction)

        else:
            tmpB = c1dbc.constrain(bid(nx, xi, xo, 0, 0, c1dbc.no_bc(), location = location),bc['x'], location = location)
            bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, tmpB, restriction = restriction)

    return bc_mat
