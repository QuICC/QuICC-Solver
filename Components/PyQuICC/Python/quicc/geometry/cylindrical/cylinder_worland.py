"""Module provides functions to generate sparse operators in a cylinder with Worland expansion in radius and Chebyshev in the vertical."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cylindrical.cylinder_radius_worland as rad
import quicc.geometry.cylindrical.cylinder_boundary_worland as cylbc


def convert_bc(bc):
    """Convert boundary dictionary into r and z kronecker product boundaries"""

    if bc['r'][0] < 0:
        bcr = bc['r']
    else:
        bcr = rad.radbc.no_bc()

        for key, val in bc['r'].items():
            if key != 0:
                bcr[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcr, bcz)

def zblk(nr, m, nz, zi, zo, qr, qz, bc, restriction = None):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.zblk(nz, zi, zo, bcz),rad.zblk(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, qr, qz, bc, restriction = restriction)

def i2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with I in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, 0, bcz), rad.i2(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 0, bc, restriction = restriction)

def i2laplh(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2laph in R kronecker with I in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, 0, bcz), rad.i2laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 0, bc, restriction = restriction)

def i2j2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i2(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 2, bc, restriction = restriction)

def i2r_1drj2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2r_1dr in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i2r_1dr(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m+1, nz, zi, zo, 1, 2, bc, restriction = restriction)

def i2laplhj2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2laph in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i2laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 2, bc, restriction = restriction)

def i2j2e1(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z of the laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d1(nz, zi, zo, bcz), rad.i2(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 2, bc, restriction = restriction)

def i2j2lapl(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z of the laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i2laplh(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i2d2(nz, zi, zo, bcz), rad.i2(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 1, 2, bc, restriction = restriction)

def i4j2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i4laplhj2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4laplh in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i4laplhj2e1(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4laplh in R kronecker with i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d1(nz, zi, zo, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i4j2lapllaplh(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i2 in Z of the laplacian of horizontal laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i2(nz, zi, zo, bcz), rad.i4lapl2h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i2d2(nz, zi, zo, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 2, bc, restriction = restriction)

def i4j4(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i4j4e1(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i4e1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4d1(nz, zi, zo, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i4drj4(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4dr in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i4dr(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m-1, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i4j4lapl2(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i4 in Z of the bilaplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i4lapl2h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz, 2.0), rad.i4laplh(nr, m, bcr), restriction = restriction)
    mat = mat + utils.restricted_kron_2d(c1d.i4d4(nz, zi, zo, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i4laplhj4(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i4laplh in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i4j4lapl(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z of the laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 2, 4, bc, restriction = restriction)

def i6j4(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i6(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 3, 4, bc, restriction = restriction)

def i6laplhj4(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i6laplh in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i6laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 3, 4, bc, restriction = restriction)

def i6laplhj4e1(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i6laplh in R kronecker with i4d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4d1(nz, zi, zo, bcz), rad.i6laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 3, 4, bc, restriction = restriction)

def i6j4lapllaplh(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z of the laplacian of horizontal laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i6lapl2h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz), rad.i6laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 3, 4, bc, restriction = restriction)

def i6j4lapl2laplh(nr, m, nz, zi, zo, bc, coeff = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z of the bilaplacia of horizontal laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, zi, zo, bcz), rad.i6lapl3h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, zi, zo, bcz, 2.0), rad.i6lapl2h(nr, m, bcr), restriction = restriction)
    mat = mat + utils.restricted_kron_2d(c1d.i4d4(nz, zi, zo, bcz), rad.i6laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, m, nz, zi, zo, 3, 4, bc, restriction = restriction)

def qid(nr, m, nz, zi, zo, qr, qz, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz, zi, zo, qz, bcz), rad.qid(nr, m, qr,bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, qr, qz, bc, restriction = restriction)

def stencil(nr, m, nz, zi, zo, bc, make_square, restriction = None):
    """Create a galerkin stencil matrix"""

    bcr, bcz = convert_bc(bc)
    mat_r = rad.stencil(nr, m, bcr, make_square)
    mat_z = c1d.stencil(nz, zi, zo, bcz, make_square)
    mat = spsp.kron(mat_z, mat_r)

    return cylbc.constrain(mat, nr, m, nz, zi, zo, 0, 0, bc, restriction = restriction)

def tau_mat_r(nr, m, nz, zi, zo, tau, kron_op, qr, qz, bc, location = 't', restriction = None):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    sr, dr, sz, dz = cylbc.convert_priority(bc.get('priority', 'r'), qr, qz)

    pad = tau.get('pad',0)
    s = tau.get('kron_shift',0)
    matR = spsp.lil_matrix((nr,nr))
    matR = rad.radbc.constrain(matR, m, tau, pad_zeros = pad, location = location)
    matZ = cylbc.bzid(nz, zi, zo, sz, dz, c1d.c1dbc.no_bc(), location = location)*kron_op(nz+s, zi, zo, bc = {0:0, 'rt':s, 'cr':s})

    matR = rad.radbc.constrain(matR, m, bcr, location = location)
    matZ = c1d.c1dbc.constrain(matZ, zi, zo, bcz, location = location)
    mat = utils.restricted_kron_2d(matZ, matR, restriction = restriction)
    return cylbc.constrain(mat, nr, m, nz, zi, zo, qr, qz, bc, restriction = restriction)

def tau_mat_z(nr, m, nz, zi, zo, tau, kron_op, qr, qz, bc, location = 't', restriction = None):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    sr, dr, sz, dz = cylbc.convert_priority(bc.get('priority', 'r'), qr, qz)

    pad = tau.get('pad',0)
    s = tau.get('kron_shift',0)
    matR = cylbc.brid(nr,m,sr,dr,rad.radbc.no_bc(), location = location)*kron_op(nr+s, m, bc = {0:0, 'rt':s, 'cr':s})
    matZ = spsp.lil_matrix((nz,nz))
    matZ = c1d.c1dbc.constrain(matZ, zi, zo, tau, pad_zeros = pad, location = location)

    matR = rad.radbc.constrain(matR, m, bcr, location = location)
    matZ = c1d.c1dbc.constrain(matZ, zi, zo, bcz, location = location)
    mat = utils.restricted_kron_2d(matZ, matR, restriction = restriction)

    return cylbc.constrain(mat, nr, m, nz, zi, zo, qr, qz, bc, restriction = restriction)
