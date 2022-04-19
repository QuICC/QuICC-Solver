"""Module provides functions to generate sparse operators in a cartesian box with two periodic directions."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import quicc.geometry.chebyshev.chebyshev_linear as cheb


def zblk(nz, zi, zo, bc, location = 't'):
    """Create a block of zeros"""

    mat = spsp.coo_matrix((nz,nz))
    return c1dbc.constrain(mat, zi, zo, bc, location = location)

def z1(nz, zi, zo, bc, coeff = 1.0, zr = 0, location = 't'):
    """Create operator for multiplication by z = [zi, zo]"""

    diags, offsets = cheb.y1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format = 'lil')
    if zr > 0 and location == 'b':
        mat[-zr:,:] = 0
    elif zr > 0 and location == 't':
        mat[0:zr,:] = 0
    mat = mat.tocoo()

    return c1dbc.constrain(mat, zi, zo, bc)

def d1(nz, zi, zo, bc, coeff = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    mat = coeff*cheb.d1(nz, zi, zo, zr)
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def d2(nz, zi, zo, bc, coeff = 1.0, zr = 2):
    """Create operator for 2nd derivative"""

    mat = coeff*cheb.d2(nz, zi, zo, zr)
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def d4(nz, zi, zo, bc, coeff = 1.0, zr = 4):
    """Create operator for 4th derivative"""

    mat = coeff*cheb.d4(nz, zi, zo, zr)
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def lapl(nz, zi, zo, k, l, bc, coeff = 1.0):
    """Create operator for horizontal laplacian"""

    mat = d2(nz, zi, zo, bc) - k**2*sid(nz, zi, zo, 2, bc) - l**2*sid(nz, zi, zo, 2, bc)

    mat = coeff*mat.tocoo()
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def laplh(nz, zi, zo, k, bc, coeff = 1.0):
    """Create operator for horizontal laplacian"""

    mat = d2(nz, zi, zo, bc) - k**2*sid(nz, zi, zo, 2, bc)

    mat = coeff*mat.tocoo()
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def lapl2h(nz, zi, zo, k, bc, coeff = 1.0):
    """Create operator for horizontal bilaplacian"""

    mat = d4(nz, zi, zo, bc) - 2.0*k**2*sid(nz, zi, zo, 4, bc)*d2(nz, zi, zo, bc) + k**4*sid(nz, zi, zo, 4, bc)

    mat = coeff*mat.tocoo()
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def i1d1(nz, zi, zo, bc, coeff = 1.0):
    """Create a quasi identity block of order 1"""

    return qid(nz, zi, zo, 1, bc, coeff)

def i1(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 1st integral in x"""

    diags, offsets = cheb.i1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i1z1(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 1st integral of multiplication by x in x"""

    diags, offsets = cheb.i1y1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i2d2(nz, zi, zo, bc, coeff = 1.0):
    """Create a quasi identity block of order 2"""

    return qid(nz, zi, zo, 2, bc, coeff)

def i2d4(nz, zi, zo, bc, coeff = 1.0):
    """Create a quasi identity block of order 2"""

    mat = coeff*i2d2(nz, zi, zo, c1dbc.no_bc())*d2(nz, zi, zo, c1dbc.no_bc())
    return c1dbc.constrain(mat, zi, zo, bc)

def i2(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 2nd integral in x"""

    diags, offsets = cheb.i2_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i2z1(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 2nd integral in x of z"""

    diags, offsets = cheb.i2y1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i2d1(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 2nd integral in x"""

    diags, offsets = cheb.i2d1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i2lapl(nz, zi, zo, k, l, bc, coeff = 1.0):
    """Create operator for 2nd integral in x of Laplacian"""

    diags, offsets = cheb.i2lapl_3d_diags(nz, zi, zo, k, l)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i2laplh(nz, zi, zo, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x of horizontal Laplacian"""

    diags, offsets = cheb.i2lapl_2d_diags(nz, zi, zo, k)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i3(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 3rd integral in x"""

    diags, offsets = cheb.i3_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 4th integral in x"""

    diags, offsets = cheb.i4_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4d1(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 4th integral in x of D_x"""

    diags, offsets = cheb.i4d1_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4d2(nz, zi, zo, bc, coeff = 1.0):
    """Create operator for 4th integral in x of D_x^2"""

    diags, offsets = cheb.i4d2_diags(nz, zi, zo)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4d4(nz, zi, zo, bc, coeff = 1.0):
    """Create a quasi identity block of order 4"""

    return qid(nz, zi, zo, 4, bc, coeff)

def i4lapl(nz, zi, zo, k, l, bc, coeff = 1.0):
    """Create operator for 4th integral in x of Laplacian"""

    diags, offsets = cheb.i4lapl_3d_diags(nz, zi, zo, k, l)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4laplh(nz, zi, zo, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x of horizontal Laplacian"""

    diags, offsets = cheb.i4lapl_2d_diags(nz, zi, zo, k)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4lapl2(nz, zi, zo, k, l, bc, coeff = 1.0):
    """Create operator for 4th integral in x of Laplacian^2"""

    diags, offsets = cheb.i4lapl2_3d_diags(nz, zi, zo, k, l)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def i4lapl2h(nz, zi, zo, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x of horizontal Laplacian^2"""

    diags, offsets = cheb.i4lapl2_2d_diags(nz, zi, zo, k)

    mat = coeff*spsp.diags(diags, offsets, format='coo')
    return c1dbc.constrain(mat, zi, zo, bc)

def qid(nz, zi, zo, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = cheb.qid(nz, zi, zo, q, coeff)
    return c1dbc.constrain(mat, zi, zo, bc)

def sid(nz, zi, zo, s, bc, coeff = 1.0):
    """Create a identity block with last s rows zeroed"""

    mat = cheb.sid(nz, zi, zo, s, coeff)
    return c1dbc.constrain(mat, zi, zo, bc, location = 'b')

def stencil(nz, zi, zo, bc, make_square):
    """Create a galerkin stencil matrix"""

    mat = qid(nz, zi, zo, 0, c1dbc.no_bc())

    if make_square:
        bc['rb'] = bc['rt']
    bc['rt'] = 0

    return c1dbc.constrain(mat, zi, zo, bc)

integral = cheb.integral

avg = cheb.avg

def tau_mat(nz, zi, zo, tau, pad, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = spsp.coo_matrix((nz,nz))
    mat = c1dbc.constrain(mat, zi, zo, tau, pad_zeros = pad)
    return c1dbc.constrain(mat, zi, zo, bc)
