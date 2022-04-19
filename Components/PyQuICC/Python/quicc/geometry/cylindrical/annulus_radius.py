"""Module provides functions to generate sparse operators for the radial direction in an annulus."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.cylindrical.annulus_radius_boundary as radbc


def zblk(nr, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return radbc.constrain(mat,bc)

def r1(nr, a, b, bc, coeff = 1.0, zr = 0):
    """Create operator for r multiplication"""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a/2.0)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, b)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, bc)

def r2(nr, a, b, bc, coeff = 1.0, zr = 0):
    """Create operator for r^2 multiplication."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = -1

    # Generate 2nd subdiagonal
    def d_2(n):
        return np.full(n.shape, a**2/4.0)

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a*b)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, (a**2 + 2.0*b**2)/2.0)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n)

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, bc)

def d1(nr, a, b, bc, coeff = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    row = [2*j for j in range(0,nr)]
    mat = spsp.lil_matrix((nr,nr))
    for i in range(0,nr-1):
        mat[i,i+1:nr:2] = row[i+1:nr:2]
    if zr > 0:
        mat[-zr:,:] = 0

    mat = coeff*(1.0/a)*mat
    return radbc.constrain(mat, bc)

def r1d1(nr, a, b, bc, coeff = 1.0, zr = 1):
    """Create operator for r times derivative"""

    mat = r1(nr, a, b, radbc.no_bc(), coeff, zr = zr)*d1(nr, a, b, radbc.no_bc(), coeff, zr = zr)
    return radbc.constrain(mat, bc)

def r1div(nr, a, b, bc, coeff = 1.0, zr = 1):
    """Create operator for r times radial divergence"""

    mat = sid(nr, zr, radbc.no_bc(), coeff) + r1(nr, a, b, radbc.no_bc(), coeff, zr = zr)*d1(nr, a, b, radbc.no_bc(), coeff, zr = zr)
    return radbc.constrain(mat, bc)

def i1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return a/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i1r1d1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral r 1st derivative T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(n - 1.0)/(2.0*n)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, b)

    # Generate 1st superdiagonal
    def d1(n):
        return a*(n + 1.0)/(2.0*n)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i1r1div(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral r radial divergence T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a/2.0)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, b)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i1r1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral of r T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 0

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n)

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b/(2.0*n)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, 0)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -d_2(n)

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i1r2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral of r^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 0

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3/(8.0*n)

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b/(2.0*n)

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2 + 4.0*b**2)/(8.0*n)

    # Generate main diagonal
    def d0(n):
        return np.full(n.shape, 0)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -d_2(n)

    # Generate 3rd superdiagonal
    def d3(n):
        return -d_3(n)

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*b/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a**3/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*b/(4.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2d2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 2nd derivative T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(n - 3.0)*(n - 2.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(n - 2.0)/n

    # Generate main diagonal
    def d0(n):
        return (a**2*n**2 - 3.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(n + 2.0)/n

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(n + 2.0)*(n + 3.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2d1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 1st derivative T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(n - 2.0)/(2.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*n + 3.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return a**2*b/((n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**2*n - 3.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(n + 2.0)/(2.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(n + 3.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 1

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4/(16.0*n*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (a**3*b)/(4.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (a**2*(2.0*b**2*n + a**2 + 2.0*b**2))/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(a**3*b)/(4.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -(a**2*(a**2 + 4.0*b**2))/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n - 1.0)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2 - 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n + 1.0)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 1.0)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r3(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-5,6)
    nzrow = 1

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5/(32.0*n*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return 3.0*a**4*b/(16.0*n*(n - 1))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(a**2*n + 3.0*a**2 + 12.0*b**2*n + 12.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(3.0*a**2 + 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*(a**2 + 6.0*b**2)/(16.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*b*(3.0*a**2 + 4.0*b**2)/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n - 1.0)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(3.0*a**2 - 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(a**2*n - 3.0*a**2 + 12.0*b**2*n - 12.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 1.0)

    # Generate 5th superdiagonal
    def d5(n):
        return d_5(n + 1.0)

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2div(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 radial divergence T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(2.0*n - 3.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*n + 2.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return a**2*b/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**2*n - 2.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(2.0*n + 3.0)/(4.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(n + 2.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2laplh(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(m - n + 2.0)*(m + n - 2.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(2.0*n - 3.0)/(2.0*n)

    # Generate main diagonal
    def d0(n):
        return (a**2*m**2 + a**2*n**2 - 2.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(2.0*n + 3.0)/(2.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(m - n - 2.0)*(m + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r2vlaplh(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 vector Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(m**2 - n**2 + 4.0*n - 3.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(2.0*n - 3.0)/(2.0*n)

    # Generate main diagonal
    def d0(n):
        return (a**2*m**2 + a**2*n**2 - a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(2.0*n + 3.0)/(2.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(m**2 - n**2 - 4.0*n - 3.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r3vlaplhr_1(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^3 vector Laplacian 1/x T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(m**2 - n**2 + 6.0*n - 8.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(2.0*n - 5.0)/(2.0*n)

    # Generate main diagonal
    def d0(n):
        return (a**2*m**2 + a**2*n**2 - 4.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(2.0*n + 5.0)/(2.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(m**2 - n**2 - 6.0*n - 8.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r3d1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^3 1st derivative T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 1

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(n - 4.0)/(16.0*n*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return 3.0*a**3*b*(n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(n - 2.0)*(a**2*n + 2.0*a**2 + 6.0*b**2*n + 6.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(3.0*a**2*n + 9.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return 3.0*a**2*(a**2 + 4.0*b**2)/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*b*(3.0*a**2*n - 9.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(n + 2.0)*(a**2*n - 2.0*a**2 + 6.0*b**2*n - 6.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3.0*a**3*b*(n + 3.0)/(8.0*n*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*(n + 4.0)/(16.0*n*(n + 1.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i2r3d1r_2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^3 1st derivative 1/r^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(n - 4.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b/(2.0*n)

    # Generate main diagonal
    def d0(n):
        return 3.0*a**2/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*b/(2.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(n + 4.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i4(nr, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return 3.0*a**4/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i4r4(nr, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-8,9)
    nzrow = 3

    # Generate 8th subdiagonal
    def d_8(n):
        return a**8/(256.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 7th subdiagonal
    def d_7(n):
        return (a**7*b)/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return (3.0*a**6*(2.0*b**2*n + a**2 + 2.0*b**2))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n - 11.0*a**2 - 4.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 - 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 - 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n + 6.0*a**2 + 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**4*(9.0*a**4*n + 33.0*a**4 + 6.0*a**2*b**2*n**2 + 120*a**2*b**2*n + 306.0*a**2*b**2 + 16.0*b**4*n**2 + 80.0*b**4*n + 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a**5*b*(3.0*a**2*n**2 - 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 - 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*(a**4*n**2 - 29.0*a**4 + 16.0*a**2*b**2*n**2 - 304.0*a**2*b**2 + 16.0*b**4*n**2 - 144.0*b**4))/(128.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (a**5*b*(3.0*a**2*n**2 + 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 + 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (a**4*(9.0*a**4*n - 33.0*a**4 - 6.0*a**2*b**2*n**2 + 120.0*a**2*b**2*n - 306.0*a**2*b**2 - 16.0*b**4*n**2 + 80.0*b**4*n - 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n - 6.0*a**2 - 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 + 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 + 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n + 11.0*a**2 + 4.0*b**2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -(3.0*a**6*(a**2 - 2.0*b**2*n + 2.0*b**2))/(64.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 7th superdiagonal
    def d7(n):
        return d_7(n + 3.0)

    # Generate 8th superdiagonal
    def d8(n):
        return d_8(n + 3.0)

    ds = [d_8, d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7, d8]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i4r4laplh(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-6,7)
    nzrow = 3

    # Generate 6th subdiagonal
    def d_6(n):
        return -a**6*(m - n + 6.0)*(m + n - 6.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -a**5*b*(2.0*m**2 - 4.0*n**2 + 41.0*n - 105.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(a**2*m**2*n - 5.0*a**2*m**2 + a**2*n**3 - 3.0*a**2*n**2 - 34.0*a**2*n + 120.0*a**2 - 2.0*b**2*m**2*n - 2.0*b**2*m**2 + 12.0*b**2*n**3 - 90.0*b**2*n**2 + 114.0*b**2*n + 216.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*b*(6.0*a**2*m**2 + 4.0*a**2*n**2 + 25.0*a**2*n - 183.0*a**2 + 16.0*b**2*n**2 - 44.0*b**2*n - 60.0*b**2)/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(a**4*m**2*n + 17.0*a**4*m**2 - a**4*n**3 + 27.0*a**4*n**2 - 8.0*a**4*n - 372.0*a**4 + 16.0*a**2*b**2*m**2*n + 32.0*a**2*b**2*m**2 + 216.0*a**2*b**2*n**2 - 360.0*a**2*b**2*n - 1584.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n - 96.0*b**4)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*b*(2.0*a**2*m**2*n - 16.0*a**2*m**2 + 4.0*a**2*n**3 - 33.0*a**2*n**2 - 55.0*a**2*n + 444.0*a**2 + 8.0*b**2*n**3 - 66.0*b**2*n**2 + 10.0*b**2*n + 348.0*b**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*(a**4*m**2*n**2 - 19.0*a**4*m**2 + a**4*n**4 - 43.0*a**4*n**2 + 396.0*a**4 + 6.0*a**2*b**2*m**2*n**2 - 54.0*a**2*b**2*m**2 + 12.0*a**2*b**2*n**4 - 336.0*a**2*b**2*n**2 + 2052.0*a**2*b**2 + 8.0*b**4*n**4 - 104.0*b**4*n**2 + 288.0*b**4)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a**3*b*(2.0*a**2*m**2*n + 16.0*a**2*m**2 + 4.0*a**2*n**3 + 33.0*a**2*n**2 - 55.0*a**2*n - 444.0*a**2 + 8.0*b**2*n**3 + 66.0*b**2*n**2 + 10.0*b**2*n - 348.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(a**4*m**2*n - 17.0*a**4*m**2 - a**4*n**3 - 27.0*a**4*n**2 - 8.0*a**4*n + 372.0*a**4 + 16.0*a**2*b**2*m**2*n - 32.0*a**2*b**2*m**2 - 216.0*a**2*b**2*n**2 - 360.0*a**2*b**2*n + 1584.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n + 96.0*b**4)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*b*(6.0*a**2*m**2 + 4.0*a**2*n**2 - 25.0*a**2*n - 183.0*a**2 + 16.0*b**2*n**2 + 44.0*b**2*n - 60.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*(a**2*m**2*n + 5.0*a**2*m**2 + a**2*n**3 + 3.0*a**2*n**2 - 34.0*a**2*n - 120.0*a**2 - 2.0*b**2*m**2*n + 2.0*b**2*m**2 + 12.0*b**2*n**3 + 90.0*b**2*n**2 + 114.0*b**2*n - 216.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*b*(2.0*m**2 - 4.0*n**2 - 41.0*n - 105.0)/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -a**6*(m - n - 6.0)*(m + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def i4r4lapl2h(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(m - n + 4.0)*(m - n + 6.0)*(m + n - 6.0)*(m + n - 4.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*b*(2.0*n - 9.0)*(2.0*m**2 - 2.0*n**2 + 18.0*n - 39.0)/(8.0*n*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(a**2*m**4 + 6.0*a**2*m**2*n - 28.0*a**2*m**2 - a**2*n**4 + 10.0*a**2*n**3 - 20.0*a**2*n**2 - 64.0*a**2*n + 192.0*a**2 + 2.0*b**2*m**2*n**2 - 4.0*b**2*m**2*n - 6.0*b**2*m**2 - 6.0*b**2*n**4 + 60.0*b**2*n**3 - 173.0*b**2*n**2 + 46.0*b**2*n + 285.0*b**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(4.0*a**2*m**2*n - 38.0*a**2*m**2 + 12.0*a**2*n**3 - 54.0*a**2*n**2 - 88.0*a**2*n + 461.0*a**2 + 16.0*b**2*n**3 - 72.0*b**2*n**2 + 24.0*b**2*n + 112.0*b**2)/(8.0*n*(n - 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*m**4 + 2.0*a**4*m**2*n**2 - 68.0*a**4*m**2 + 3.0*a**4*n**4 - 68.0*a**4*n**2 + 416.0*a**4 + 8.0*a**2*b**2*m**2*n**2 - 32.0*a**2*b**2*m**2 + 24.0*a**2*b**2*n**4 - 332.0*a**2*b**2*n**2 + 944.0*a**2*b**2 + 8.0*b**4*n**4 - 40.0*b**4*n**2 + 32.0*b**4)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(4.0*a**2*m**2*n + 38.0*a**2*m**2 + 12.0*a**2*n**3 + 54.0*a**2*n**2 - 88.0*a**2*n - 461.0*a**2 + 16.0*b**2*n**3 + 72.0*b**2*n**2 + 24.0*b**2*n - 112.0*b**2)/(8.0*n*(n - 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2*m**4 - 6.0*a**2*m**2*n - 28.0*a**2*m**2 - a**2*n**4 - 10.0*a**2*n**3 - 20.0*a**2*n**2 + 64.0*a**2*n + 192.0*a**2 + 2.0*b**2*m**2*n**2 + 4.0*b**2*m**2*n - 6.0*b**2*m**2 - 6.0*b**2*n**4 - 60.0*b**2*n**3 - 173.0*b**2*n**2 - 46.0*b**2*n + 285.0*b**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*b*(2.0*n + 9.0)*(2.0*m**2 - 2.0*n**2 - 18.0*n - 39.0)/(8.0*n*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*(m - n - 6.0)*(m - n - 4.0)*(m + n + 4.0)*(m + n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, bc)

def qid(nr, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = spsp.coo_matrix((nr,nr))
    if coeff != 1.0:
        mat.data = coeff*np.ones((nr-q))
    else:
        mat.data = np.ones((nr-q))
    mat.row = np.arange(q,nr)
    mat.col = mat.row
    return radbc.constrain(mat, bc)

def sid(nr, s, bc, coeff = 1.0):
    """Create a identity block with last s rows zeroed"""

    mat = spsp.coo_matrix((nr,nr))
    if coeff != 1.0:
        mat.data = coeff*np.ones((nr-s))
    else:
        mat.data = np.ones((nr-s))
    mat.row = np.arange(0,nr-s)
    mat.col = mat.row
    return radbc.constrain(mat, bc)

def linear_r2x(ro, rratio):
    """Calculat a and b for linear map r = a*x + b"""

    b = (ro*rratio + ro)/2.0;
    a = ro - b;

    return (a, b)
