"""Module provides functions to generate the boundary conditions in a cartesian 3D geometry"""

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

    return {'x':c1dbc.no_bc(), 'y':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

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

def constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, qx, qy, qz, bc, location = 't', restriction = None):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'xz')
    sx = 2*[0]
    dx = 2*[0]
    sy = 2*[0]
    dy = 2*[0]
    sz = 2*[0]
    dz = 2*[0]
    if priority == 'xy':
        sx = [qx,qx]
        sy = [0,qy]
    elif priority == 'xz':
        sx = [qx,qx]
        sz = [0,qz]
    elif priority == 'yx':
        sx = [0,qx]
        sy = [qy,qy]
    elif priority == 'yz':
        sy = [qy,qy]
        sz = [qz,0]
    elif priority == 'zx':
        sx = [qx,0]
        sz = [qz,qz]
    elif priority == 'zy':
        sy = [qy,0]
        sz = [qz,qz]
    elif priority == 'xsy':
        sx = [qx,qx]
        sy = [0,qy]
        dy = [0,-qy]
    elif priority == 'xsz':
        sx = [qx,qx]
        sz = [0,qz]
        dz = [0,-qz]
    elif priority == 'ysx':
        sx = [0,qx]
        dx = [0,-qx]
        sy = [qy,qy]
    elif priority == 'ysz':
        sy = [qy,qy]
        sz = [qz,0]
        dz = [-qz,0]
    elif priority == 'zsx':
        sx = [qx,0]
        dx = [-qx,0]
        sz = [qz,qz]
    elif priority == 'zsy':
        sy = [qy,0]
        dy = [-qy,0]
        sz = [qz,qz]
    elif priority == 'n':
        sx = [qx,qx]
        sy = [qy,qy]
        sz = [qz,qz]
    else:
        raise RuntimeError("Unknown boundary condition priority!")

    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, xi, xo, bc['x'], location = location)
        bc_mat = bc_mat + utils.restricted_kron_3d(bid(ny, yi, yo, sy[0], dy[0], bc['y'], location = location), bid(nz, zi, zo, sz[0], dz[0], bc['z']), bcMat, restriction = restriction)

    if bc['y'][0] > 0:
        bcMat = spsp.lil_matrix((ny,ny))
        bcMat = c1dbc.constrain(bcMat, yi, yo, bc['y'], location = location)
        bc_mat = bc_mat + utils.restricted_kron_3d(bcMat, bid(nz, zi, zo, sz[1], dz[1], bc['z']), bid(nx, xi, xo, sx[0], dx[0], bc['x'], location = location), restriction = restriction)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, zi, zo, bc['z'], location = location)
        bc_mat = bc_mat + utils.restricted_kron_3d(bid(ny, yi, yo, sy[1], dy[1], bc['y'], location = location), bcMat, bid(nx, xi, xo, sx[1], dx[1], bc['x']), restriction = restriction)

    return bc_mat
