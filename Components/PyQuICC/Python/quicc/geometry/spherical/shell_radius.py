"""Module provides functions to generate sparse operators for the radial direction in a spherical shell."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.shell_radius_boundary as radbc
import quicc.geometry.chebyshev.chebyshev_linear as cheb


def zblk(nr, ri, ro, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return radbc.constrain(mat, ri, ro, bc)

def r1(nr, ri, ro, bc, coeff = 1.0, zr = 0):
    """Create operator for r multiplication"""

    diags, offsets = cheb.y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, ri, ro, bc)

def r2(nr, ri, ro, bc, coeff = 1.0, zr = 0):
    """Create operator for r^2 multiplication"""

    diags, offsets = cheb.y2_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, ri, ro, bc)

def r4(nr, ri, ro, bc, coeff = 1.0, zr = 0):
    """Create operator for r^4 multiplication"""

    diags, offsets = cheb.y4_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, ri, ro, bc)

def d1(nr, ri, ro, bc, coeff = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    mat = coeff*cheb.d1(nr, ri, ro, zr)
    return radbc.constrain(mat, ri, ro, bc, location = 'b')

def d2(nr, ri, ro, bc, coeff = 1.0, zr = 2):
    """Create operator for 2nd derivative"""

    mat = coeff*cheb.d2(nr, ri, ro, zr)
    return radbc.constrain(mat, ri, ro, bc, location = 'b')

def i1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 1st integral T_n(x)."""

    diags, offsets = cheb.i1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i1r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 1st integral r T_n(x)."""

    diags, offsets = cheb.i1y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    diags, offsets = cheb.i2_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r2(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 T_n(x)."""

    diags, offsets = cheb.i2y2_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r3(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^3 T_n(x)."""

    diags, offsets = cheb.i2y3_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r2lapl(nr, ri, ro, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 Laplacian T_n(x)."""

    diags, offsets = cheb.i2y2sphlapl_diags(nr, ri, ro, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r3lapl(nr, ri, ro, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^3 Laplacian T_n(x)."""

    diags, offsets = cheb.i2y3sphlapl_diags(nr, ri, ro, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    diags, offsets = cheb.i4_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r T_n(x)."""

    diags, offsets = cheb.i4y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r1d1r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r D r T_n(x)."""

    diags, offsets = cheb.i4y1d1y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r3d2(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^3 D^2 T_n(x)."""

    diags, offsets = cheb.i4y3d2_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r4(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 T_n(x)."""

    diags, offsets = cheb.i4y4_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r4lapl(nr, ri, ro, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian T_n(x)."""

    diags, offsets = cheb.i4y4sphlapl_diags(nr, ri, ro, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r2lapl2_l1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^2 Laplacian^2 T_n(x) for l = 1."""

    diags, offsets = cheb.i4y2sphlapl2_l1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r4lapl2(nr, ri, ro, l, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 Laplacian^2 T_n(x)."""

    diags, offsets = cheb.i4y4sphlapl2_diags(nr, ri, ro, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral of r T_n(x)."""

    diags, offsets = cheb.i2y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r1d1r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral of r D_r r T_n(x)."""

    diags, offsets = cheb.i2y1d1y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i2r2d1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 2nd integral of r^2 D_r T_n(x)."""

    diags, offsets = cheb.i2y2d1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r2(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^2 T_n(x)."""

    diags, offsets = cheb.i4y2_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r3(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^3 T_n(x)."""

    diags, offsets = cheb.i4y3_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r3d1r1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^3 D r T_n(x)."""

    diags, offsets = cheb.i4y3d1y1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def i4r4d1(nr, ri, ro, bc, coeff = 1.0):
    """Create operator for 4th integral of r^4 D_r T_n(x)."""

    diags, offsets = cheb.i4y4d1_diags(nr, ri, ro)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, ri, ro, bc)

def qid(nr, ri, ro, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = cheb.qid(nr, ri, ro, q, coeff)
    return radbc.constrain(mat, ri, ro, bc)

def linear_r2x(ro, rratio):
    """Calculat a and b for linear map r = a*x + b"""

    b = (ro*rratio + ro)/2.0;
    a = ro - b;

    return (a, b)

def stencil(nr, ri, ro, bc, make_square):
    """Create a galerkin stencil matrix"""

    mat = qid(nr, ri, ro, 0, radbc.no_bc())

    if not make_square:
        bc['rt'] = 0

    return radbc.constrain(mat, ri, ro, bc)

integral = cheb.integral

avg = cheb.avg
