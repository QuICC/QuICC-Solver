"""Module provides functions to generate sparse operators in a cartesian box."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_2d as c2d
import quicc.geometry.cartesian.cartesian_boundary_3d as c3dbc
import quicc.base.utils as utils


def convert_bc(bc):
    """Convert boundary dictionary into x, y and z kronecker product boundaries"""

    if bc['x'][0] < 0:
        bcx = bc['x']
    else:
        bcx = c1d.c1dbc.no_bc()
        for key, val in bc['x'].items():
            if key != 0:
                bcx[key] = val

    if bc['y'][0] < 0:
        bcy = bc['y']
    else:
        bcy = c1d.c1dbc.no_bc()
        for key, val in bc['y'].items():
            if key != 0:
                bcy[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcx, bcy, bcz)

def d1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, sx = 1, sy = 0, sz = 0, restriction = None):
    """Create operator for 1st X derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.sid(ny, yi, yo, sy, bcy), c1d.sid(nz, zi, zo, sz, bcz), c1d.d1(nx, xi, xo, bcx, zr = sx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, sx, sy, sz, bc, location = 'b', restriction = restriction)

def e1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, sx = 0, sy = 1, sz = 0, restriction = None):
    """Create operator for 1st Y derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.d1(ny, yi, yo, bcy, zr = sy), c1d.sid(nz, zi, zo, sz, bcz), c1d.sid(nx, xi, xo, sx, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, sx, sy, sz, bc, location = 'b', restriction = restriction)

def f1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, sx = 0, sy = 0, sz = 1, restriction = None):
    """Create operator for 1st Z derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.sid(ny, yi, yo, sy, bcy), c1d.d1(nz, zi, zo, bcz, zr = sz), c1d.sid(nx, xi, xo, sx, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, sx, sy, sz, bc, location = 'b', restriction = restriction)

def zblk(nx, xi, xo, ny, yi, yo, nz, zi, zo, qx, qy, qz, bc, location = 't', restriction = None):
    """Create a block of zeros"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = utils.restricted_kron_3d(c1d.zblk(ny, yi, yo, bcy), c1d.zblk(nz, zi, zo, bcz), c1d.zblk(nx, xi, xo, bcx), None)
    return c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, qx, qy, qz, bc, location = location, restriction = restriction)

def i1j1k1d1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral in x,y,z of T'_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i1(ny, yi, yo, bcy), c1d.i1(nz, zi, zo, bcz), c1d.i1d1(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 1, 1, 1, bc, restriction = restriction)

def i1j1k1e1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral in x,y,z of T_n(x)T'_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i1d1(ny, yi, yo, bcy), c1d.i1(nz, zi, zo, bcz), c1d.i1(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 1, 1, 1, bc, restriction = restriction)

def i1j1k1f1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral in x,y,z of T_n(x)T_n(y)T'_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i1(ny, yi, yo, bcy), c1d.i1d1(nz, zi, zo, bcz), c1d.i1(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 1, 1, 1, bc, restriction = restriction)

def i2j2k2d2e2f2(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block of order 2,2,2"""

    return qid(nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, coeff, restriction = restriction)

def i2j2k2(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i2j2k2d1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of T'_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2d1(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i2j2k2e1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of T_n(x)T'_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i2d1(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i2j2k2f1(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T'_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2d1(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i2j2k2laplh(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of horizontal Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2d2(nx, xi, xo, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_3d(c1d.i2d2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i2j2k2lapl(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2d2(nx, xi, xo, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_3d(c1d.i2d2(ny, yi, yo, bcy), c1d.i2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + utils.restricted_kron_3d(c1d.i2(ny, yi, yo, bcy), c1d.i2d2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 2, 2, 2, bc, restriction = restriction)

def i4j4k4(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 4, 4, 4, bc, restriction = restriction)

def i4j4k4lapl(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz), c1d.i4d2(nx, xi, xo, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_3d(c1d.i4d2(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4d2(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 4, 4, 4, bc, restriction = restriction)

def i4j4k4lapl2(nx, xi, xo, ny, yi, yo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,y,z of Laplacian^2 T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz), c1d.i4d4(nx, xi, xo, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_3d(c1d.i4d4(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4d4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + 2*utils.restricted_kron_3d(c1d.i4d2(ny, yi, yo, bcy), c1d.i4(nz, zi, zo, bcz),c1d.i4d2(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + 2*utils.restricted_kron_3d(c1d.i4d2(ny, yi, yo, bcy), c1d.i4d2(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = mat + 2*utils.restricted_kron_3d(c1d.i4(ny, yi, yo, bcy), c1d.i4d2(nz, zi, zo, bcz), c1d.i4d2(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, 4, 4, 4, bc, restriction = restriction)

def qid(nx, xi, xo, ny, yi, yo, nz, zi, zo, qx, qy, qz, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block order qx,qy,qz"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_3d(c1d.qid(ny, yi, yo, qy, bcy), c1d.qid(nz, zi, zo, qz, bcz), c1d.qid(nx, xi, xo, qx, bcx), restriction = restriction)
    return c3dbc.constrain(mat, nx, xi, xo, ny, yi, yo, nz, zi, zo, qx, qy, qz, bc, restriction = restriction)

def volumeAvg(nx, xi, xo, ny, yi, yo, nz, zi, zo):
    """Compute a volume average"""

    mat = c1d.avg(ny, yi, yo)*spsp.kron(c1d.qid(ny, yi, yo, 0, c1d.c1dbc.no_bc()), c2d.surfaceAvg(nx, xi, xo, nz, zi, zo))
    return mat

def avgFlux_z(nx, xi, xo, ny, yi, yo, nz, zi, zo):
    """Compute average flux through Z surface"""

    mat = spsp.kron(c1d.avg(ny, yi, yo), spsp.kron(c1d.surfaceFlux(nz, zi, zo), c1d.avg(nx, xi, xo)))
    return mat
