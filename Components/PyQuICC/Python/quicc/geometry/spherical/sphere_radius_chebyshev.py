"""Module provides functions to generate sparse operators for the radial direction in a sphere with Chebyshev expansion."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.sphere_radius_boundary_chebyshev as radbc


def zblk(nr, l, bc):
    """Create a block of zeros"""

    mat = spsp.coo_matrix((nr,nr))
    return radbc.constrain(mat,l,bc)

def r1(nr, l, bc, coeff = 1.0, zr = 0):
    """Create operator for r multiplication"""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, 1.0/2.0)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, l, bc)

def r2(nr, l, bc, coeff = 1.0, zr = 0):
    """Create operator for r^2 multiplication"""

    ns = np.arange(l%2, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 2nd subdiagonal
    def d_1(n):
        return np.full(n.shape, 1.0/4.0)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, 1.0/2.0)

    # Generate 2nd superdiagonal
    def d1(n):
        return d_1(n)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, l, bc)

def d1(nr, l, bc, coeff = 1.0):
    """Create operator for 1st derivative"""

    row = [2*j for j in range(l%2,2*nr,2)]
    mat = spsp.lil_matrix((nr,nr))
    for i in range(0,nr-1):
        mat[i,i+(l+1)%2:] = row[i+(l+1)%2:]

    mat = coeff*mat.tocoo()
    return radbc.constrain(mat, l, bc)

def i1(nr, l, bc, coeff = 1.0):
    """Create operator for 1st integral T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    ns = np.arange(l%2, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_1(n):
        return 1.0/(4.0*n*(n - 1.0))

    # Generate diagonal
    def d0(n):
        return -1.0/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d1(n):
        return d_1(n+1.0)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2r1d1r1(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r D_r r T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_2(n):
        return (n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (n + 2.0)/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_2(n)

    # Generate 3rd superdiagonal
    def d2(n):
        return -d_1(n)

    ds = [d_2, d_1, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2r2(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(16.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -1.0/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2r2lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 Laplacian T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return  -(l - n + 2.0)*(l + n - 1.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (l**2 + l + n**2 - 1.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(l - n - 1.0)*(l + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r3d1r1(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^3 D r T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-3,5)
    else:
        offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 7th subdiagonal
    def d_4(n):
        return (n - 6.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_3(n):
        return (n**2 + 7.0*n - 54.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_2(n):
        return -3.0*(n**2 - 4.0*n - 28.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*(n**3 + 9.0*n**2 - 24.0*n - 156.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*(n**3 - 9.0*n**2 - 24.0*n + 156.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d2(n):
        return 3.0*(n**2 + 4.0*n - 28.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 5th superdiagonal
    def d3(n):
        return -(n**2 - 7.0*n - 54.0)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d4(n):
        return -(n + 6.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r2(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^2 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return 1.0/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(n - 5.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(n + 17.0)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return (n**2 - 19.0)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(n - 17.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n + 5.0)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return 1.0/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r4(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return 1.0/(256.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return 3.0/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(n**2 - 19.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(3.0*(3.0*n + 11.0))/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*(n**2 - 29.0))/(128.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (3.0*(3.0*n - 11.0))/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n**2 - 19.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -d_3(n + 2.0)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 3.0)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r4lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(l - n + 6.0)*(l + n - 5.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (l**2*n - 5.0*l**2 + l*n - 5.0*l + n**3 - 3.0*n**2 - 28.0*n + 96.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (l**2*n + 17.0*l**2 + l*n + 17.0*l - n**3 + 24.0*n**2 - 5.0*n - 294.0)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -(l**2*n**2 - 19.0*l**2 + l*n**2 - 19.0*l + n**4 - 37.0*n**2 + 312.0)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (l**2*n - 17.0*l**2 + l*n - 17.0*l - n**3 - 24.0*n**2 - 5.0*n + 294.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (l**2*n + 5.0*l**2 + l*n + 5.0*l + n**3 + 3.0*n**2 - 28.0*n - 96.0)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(l - n - 5.0)*(l + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r2lapl2_l1(nr, bc, coeff = 1.0):
    """Create operator for 4th integral of r^2 Laplacian^2 T_n(x) for l = 1."""

    parity = 1
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 3

    # Generate 1st subdiagonal
    def d_1(n):
        return (n - 5.0)/(4.0*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (n**2 + 3.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (n + 5.0)/(4.0*(n + 1.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, 1, bc)

def i4r4lapl2(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian^2 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return ((l - n + 4.0)*(l - n + 6.0)*(l + n - 5.0)*(l + n - 3.0))/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -((l - n + 4.0)*(l + n - 3.0)*(l**2 + l + n**2 - 2.0*n - 9.0))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*l**4 + 6.0*l**3 + 2.0*l**2*n**2 - 47.0*l**2 + 2.0*l*n**2 - 50.0*l + 3.0*n**4 - 51.0*n**2 + 228.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -((l - n - 3.0)*(l + n + 4.0)*(l**2 + l + n**2 + 2.0*n - 9.0))/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return ((l - n - 5.0)*(l - n - 3.0)*(l + n + 4.0)*(l + n + 6.0))/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2r1(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -d_2(n + 1.0)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_2(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    ds = [d_2, d_1, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2r2d1(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 D T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return (n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (n + 3.0)/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_2(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -d_1(n)

    ds = [d_2, d_1, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r3(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^3 T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-3,5)
    else:
        offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return 1.0/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(n - 11.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -3.0*(n + 6.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 3.0*(n**2 - 5.0*n - 34.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*(n**2 + 5.0*n - 34.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -3.0*(n - 6.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(n + 11.0)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return 1.0/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4r4d1(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 D T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-3,5)
    else:
        offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return (n - 7.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (n - 5.0)*(n + 13.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -3.0*(n**2 - 5.00*n - 34.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*(n**3 + 10.0*n**2 - 29.0*n - 190.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*(n**3 - 10.0*n**2 - 29.0*n + 190.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 3.0*(n**2 + 5.0*n - 34.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(n - 13.0)*(n + 5.0)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(n + 7.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def qid(nr, l, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = spsp.coo_matrix((nr,nr))
    if coeff != 1.0:
        mat.data = coeff*np.ones((nr-q))
    else:
        mat.data = np.ones((nr-q))
    mat.row = np.arange(q,nr)
    mat.col = mat.row
    return radbc.constrain(mat, l, bc)

def stencil(nr, l, bc, make_square):
    """Create a galerkin stencil matrix"""

    mat = qid(nr, l, 0, radbc.no_bc())

    if not make_square:
        bc['rt'] = 0

    return radbc.constrain(mat, l, bc)

def integral(nr, l):
    """Compute the definite integral of the expansion"""

    mat = spsp.lil_matrix((1,nr))
    if l%2 == 0:
        mat[0,:] = [2.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,2*nr,2)]
        mat[0,0] = mat[0,0]/2.0
    else:
        mat[0,1:] = [2.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0) + (-1.0)**((n-1)//2)*n/(n**2 - 1.0)) for n in np.arange(3,2*nr,2)]
        mat[0,0] = 1

    return mat

def avg(nr, l):
    """Compute the average of the expansion"""

    mat = integral(nr,l)/2.0

    return mat
