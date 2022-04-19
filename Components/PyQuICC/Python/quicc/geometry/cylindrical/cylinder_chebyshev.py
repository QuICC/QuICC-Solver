"""Module provides functions to generate sparse operators in a cylinder with Chebyshev expansion in radius."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cylindrical.cylinder_radius_chebyshev as rad
import quicc.geometry.cylindrical.cylinder_boundary_chebyshev as cylbc


def convert_bc(bc):
    """Convert boundary dictionary into x and z kronecker product boundaries"""

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

def r1d1(nr, nz, parity, bc, coeff = 1.0, sr = 1, sz = 0, restriction = None):
    """Create a r1d1 in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz, sz, bcz), rad.r1d1(nr, parity, bcr, zr = sr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, sr, sz, bc, location = 'b')

def r1div(nr, nz, parity, bc, coeff = 1.0, sr = 1, sz = 0, restriction = None):
    """Create a xdiv in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz, sz, bcz), rad.r1div(nr, parity, bcr, zr = sr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, sr, sz, bc, location = 'b')

def r1e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0, sr = 0, sz = 1, restriction = None):
    """Create operator for x in R and 1st derivative in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.d1(nz, bcz, cscale = zscale, zr = sz), rad.r1(nr, parity, bcr, zr = sr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, sr, sz, bc, location = 'b')

def r2e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0, sr = 0, sz = 1, restriction = None):
    """Create operator for r^2 in R and 1st derivative in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.d1(nz, bcz, cscale = zscale, zr = sz), rad.r2(nr, parity, bcr, zr = sr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, sr, sz, bc, location = 'b')

def zblk(nr, nz, parity, qr, qz, bc, restriction = None):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,bcz),rad.zblk(nr,parity,bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, qr, qz, bc)

def i1j1(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i1 in R kronecker with i1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz, bcz), rad.i1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 0, 1, bc)

def i1j1r1d1(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i1r1d1 in R kronecker with i1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz, bcz), rad.i1r1d1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 0, 1, bc)

def i1j1r1div(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i1r1div in R kronecker with i1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz, bcz), rad.i1r1div(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 0, 1, bc)

def i1j1r1e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i1r1 in R kronecker with i1d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1d1(nz, bcz, cscale = zscale), rad.i1r1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 0, 1, bc)

def i1j1r2e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i1r2 in R kronecker with i1d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1d1(nz, bcz, cscale = zscale), rad.i1r2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 0, 1, bc)

def i2j2(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz, bcz), rad.i2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r1(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i2r1 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz, bcz), rad.i2r1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i2r2 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r3(nr, nz, parity, bc, coeff = 1.0, restriction = None):
    """Create a i2r2 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r3(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2d1(nr, nz, parity, bc, coeff = 1.0):
    """Create a i2r2d1 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r2d1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r3d1(nr, nz, parity, bc, coeff = 1.0):
    """Create a i2r3d1 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r3d1(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r3d1r_2(nr, nz, parity, bc, coeff = 1.0):
    """Create a i2r3d1r_2 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r3d1r_2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2 in R kronecker with an i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(nz,bcz, cscale = zscale), rad.i2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2e1(nr, nz, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2r2 in R kronecker with an i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(nz,bcz, cscale = zscale), rad.i2r2(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2div(nr, nz, parity, bc, coeff = 1.0):
    """Create a i2r2div in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2r2div(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2r2laplh(nr, nz, m, parity, bc, coeff = 1.0):
    """Create a i2r2laplh in R kronecker with identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz, 0 ,bcz), rad.i2r2laplh(nr, m, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2lapl(nr, nz, m, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2r2lapl in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), rad.i2r2laplh(nr, m, parity, bcr), format ='coo')
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i2d2(nz, bcz, cscale = zscale), rad.i2r2(nr, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r2vlapl(nr, nz, m, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2r2vlapl in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), rad.i2r2vlaplh(nr, m, parity, bcr), format ='coo')
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i2d2(nz, bcz, cscale = zscale), rad.i2r2(nr, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i2j2r3vlaplr_1(nr, nz, m, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2r3vlaplr_1 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), rad.i2r3vlaplhr_1(nr, m, parity, bcr), format ='coo')
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i2d2(nz, bcz, cscale = zscale), rad.i2r2(nr, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 1, 2, bc)

def i4j4(nr, nz, parity, bc, coeff = 1.0):
    """Create a i4 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4j4r4(nr, nz, parity, bc, coeff = 1.0):
    """Create a i4r4 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4r4(nr, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4r4laplh(nr, nz, m, parity, bc, coeff = 1.0):
    """Create a i4r4lapl in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz, 0, bcz), rad.i4r4laplh(nr, m, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4j4r4lapl(nr, nz, m, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i4r4lapl in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), rad.i4r4laplh(nr, m, parity, bcr), format ='coo')
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i4d2(nz,bcz, cscale = zscale), rad.i4r4(nr, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4r4lapl2h(nr, nz, m, parity, bc, coeff = 1.0):
    """Create a i4r4lapl2 in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz, 0, bcz), rad.i4r4lapl2h(nr, m, parity, bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4j4r4lapl2h(nr, nz, m, parity, bc, coeff = 1.0):
    """Create a i4r4lapl2 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.qid(nz, 0, bcz), rad.i4r4lapl2h(nr, m, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def i4j4r4lapl2(nr, nz, m, parity, bc, coeff = 1.0, zscale = 1.0):
    """Create a i4r4lapl2 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), rad.i4r4lapl2h(nr, m, parity, bcr), format ='coo')
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + 2.0*spsp.kron(c1d.i4d2(nz, bcz, cscale = zscale), rad.i4r4laplh(nr, m, parity, bcr), format ='coo')
    mat = mat + spsp.kron(c1d.i4d4(nz, bcz, cscale = zscale), rad.i4r4(nr, parity, bcr), format ='coo')
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, parity, 2, 4, bc)

def qid(nr, nz, parity, qr, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), rad.qid(nr, parity, qr,bcr), format ='coo')
    return cylbc.constrain(mat, nr, nz, parity, qr, qz, bc)
