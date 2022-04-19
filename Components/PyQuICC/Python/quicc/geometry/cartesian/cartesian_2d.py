"""Module provides functions to generate sparse operators in a cartesian box with a single periodic dimension."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_boundary_2d as c2dbc
import quicc.base.utils as utils


def convert_bc(bc):
    """Convert boundary dictionary into x and z kronecker product boundaries"""

    if bc['x'][0] < 0:
        bcx = bc['x']
    else:
        bcx = c1d.c1dbc.no_bc()
        for key, val in bc['x'].items():
            if key != 0:
                bcx[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcx, bcz)

def d1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, sx = 1, sz = 0, restriction = None):
    """Create operator for the 1st Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.sid(nz, zi, zo, sz, bcz), c1d.d1(nx, xi, xo, bcx, zr = sx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, sx, sz, bc, location = 'b', restriction = restriction)

def e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, sx = 0, sz = 1, restriction = None):
    """Create operator for the 1st Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.d1(nz, zi, zo, bcz, zr = sz), c1d.sid(nx, xi, xo, sx, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, sx, sz, bc, location = 'b', restriction = restriction)

def lapl(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for the 2nd X derivative and 2nd Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.sid(nz, zi, zo, 2, bcz), c1d.laplh(nx, xi, xo, k, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.d2(nz, zi, zo, bcz), c1d.sid(nx, xi, xo, 2, bcx), restriction = restriction)
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, location = 'b', restriction = restriction)

def laplh(nx, xi, xo, nz, zi, zo, k, sz, bc, coeff = 1.0, restriction = None):
    """Create operator for the horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.sid(nz, zi, zo, sz, bcz), c1d.laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, sz, bc, location = 'b', restriction = restriction)

def lapl2h(nx, xi, xo, nz, zi, zo, k, sz, bc, coeff = 1.0, restriction = None):
    """Create operator for the horizontal bilaplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.sid(nz, zi, zo, sz, bcz), c1d.lapl2h(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, sz, bc, location = 'b', restriction = restriction)

def zblk(nx, xi, xo, nz, zi, zo, qx, qz, bc, location = 't', restriction = None):
    """Create a block of zeros"""

    bcx, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.zblk(nz, zi, zo, bcz),c1d.zblk(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, qx, qz, bc, location = location, restriction = restriction)

def i1j1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral in X and 1st integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 1, 1, bc, restriction = restriction)

def i1j1d1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral of 1st derivative in X and 1st integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i1d1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 1, 1, bc, restriction = restriction)

def i1j1e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 1st integral in X and 1st integral of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1d1(nz, zi, zo, bcz), c1d.i1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 1, 1, bc, restriction = restriction)

def i2j1e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X and 1st integral of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1d1(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 1, bc, restriction = restriction)

def i2j2d1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral of 1st derivative in X and 2nd integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i2d1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2j2e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X and 2nd integrazl of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d1(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2j2d1e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral of 1st derivative in X and 2nd integrazl of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d1(nz, zi, zo, bcz), c1d.i2d1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2j2d2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X and 2nd integrazl of 2nd derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i2d2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2j2e2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X and 2nd integrazl of 2nd derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2j2d2e2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block of order 2,2"""

    return qid(nx, xi, xo, nz, zi, zo, 2, 2, bc, coeff, restriction = restriction)

def i2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, 0, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 0, bc, restriction = restriction)

def i2j1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x and 1st integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 1, bc, restriction = restriction)

def i2j2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in X,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i2laplh(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x and 1st integral in z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, 0, bcz), c1d.i2laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 0, bc, restriction = restriction)

def i2j1laplh(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x and 1st integral in z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i2laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 1, bc, restriction = restriction)

def i2j2laplh(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i2laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc)

def i2j2lapl(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i2laplh(nx, xi, xo, k, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i2d2(nz, zi, zo, bcz), c1d.i2(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i4j1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x and 1st integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 1, bc, restriction = restriction)

def i4j2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x and 2nd integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 2, bc, restriction = restriction)

def i4j1e1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x and 1st integral of 1st derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1d1(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 1, bc, restriction = restriction)

def i4j2e2(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x, 2nd integral of 2nd derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d2(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 2, bc, restriction = restriction)

def i4j4(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 4, bc, restriction = restriction)

def i4j4d1(nx, xi, xo, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral of 1st derivative in X and 4th integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), c1d.i4d1(nx, xi, xo, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 4, bc, restriction = restriction)

def i4j1laplh(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x 1st integral of z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i4laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 1, bc, restriction = restriction)

def i4j2laplh(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x 2nd integral of z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i4laplh(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 2, bc, restriction = restriction)

def i4j4lapl(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), c1d.i4laplh(nx, xi, xo, k, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 4, bc, restriction = restriction)

def i4j1lapl2h(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x and 1st integral of z of horizontal Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i1(nz, zi, zo, bcz), c1d.i4lapl2h(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 1, bc, restriction = restriction)

def i4j2lapl2h(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x and 2nd integral in z of horizontal Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), c1d.i4lapl2h(nx, xi, xo, k, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 2, bc, restriction = restriction)

def i4j4lapl2(nx, xi, xo, nz, zi, zo, k, bc, coeff = 1.0, restriction = None):
    """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), c1d.i4lapl2h(nx, xi, xo, k, bcx), restriction = restriction)
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + 2.0*utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz), c1d.i4laplh(nx, xi, xo, k, bcx), restriction = restriction)
    mat = mat + utils.restricted_kron_2d(c1d.i4d4(nz, zi, zo, bcz), c1d.i4(nx, xi, xo, bcx), restriction = restriction)
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 4, 4, bc, restriction = restriction)

def qid(nx, xi, xo, nz, zi, zo, qx, qz, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block order qx in x in qz in z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, qz, bcz), c1d.qid(nx, xi, xo, qx, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, qx, qz, bc, restriction = restriction)

def sid(nx, xi, xo, nz, zi, zo, sx, sz, bc, coeff = 1.0, restriction = None):
    """Create a identity block order with last sx, sz rows zeroed"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.sid(nz, zi, zo, sz, bcz), c1d.sid(nx, xi, xo, sx, bcx), restriction = restriction)
    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, sx, sz, bc, restriction = restriction)

def surfaceAvg(nx, xi, xo, nz, zi, zo):
    """Compute a surface average"""

    mat = c1d.avg(nz)*spsp.kron(c1d.qid(nz, zi, zo, 0, c1d.c1dbc.no_bc()), c1d.avg(nx, xi, xo))
    return mat

def stencil(nx, xi, xo, nz, zi, zo, bc, make_square, restriction = None):
    """Create a galerkin stencil matrix"""

    bcx, bcz = convert_bc(bc)
    mat_x = c1d.stencil(nx, xi, xo, bcx, make_square)
    mat_z = c1d.stencil(nz, zi, zo, bcz, make_square)
    mat = spsp.kron(mat_z, mat_x)

    return c2dbc.constrain(mat, nx, xi, xo, nz, zi, zo, 0, 0, bc, restriction = restriction)
