"""Module provides functions to generate sparse operators for the radial direction in a sphere with Worland expansion."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.sphere_radius_boundary_worland as radbc
import quicc.geometry.worland.wnl as wnl



def zblk(nr, l, bc):
    """Create a block of zeros"""

    # Copy BC dict as we modify it!
    bc = dict(bc)

    mat = spsp.coo_matrix((nr,nr))
    return radbc.constrain(mat,l,bc)

def r2(nr, l, bc, coeff = 1.0, zr = 0):
    """Create operator for 1st integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    diags,offsets = wnl.r2_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocoo()
    return radbc.constrain(mat, l, bc)

def i1(nr, l, bc, coeff = 1.0):
    """Create operator for 1st integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    diags,offsets = wnl.i1_diags(nr+1, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i1qm(nr, l, bc, coeff = 1.0):
    """Create operator for 1st integral of Q r^{l-1} P_n^{-1/2,l-3/2}(2r^2 -1)."""

    assert(l > 0)

    diags,offsets = wnl.i1qm_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i1qp(nr, l, bc, coeff = 1.0):
    """Create operator for 1st integral of Q r^{l+1} P_n^{-1/2,l+1/2}(2r^2 -1)."""

    diags,offsets = wnl.i1qp_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i2(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    diags,offsets = wnl.i2_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i2lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    diags,offsets = wnl.i2lapl_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i2qm(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of Q r^{l-1} P_n^{-1/2,l-3/2}(2r^2 -1)."""

    diags,offsets = wnl.i2qm_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i2qp(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of Q r^{l+1} P_n^{-1/2,l+1/2}(2r^2 -1)."""

    diags,offsets = wnl.i2qp_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, l, bc)

def i4(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    diags,offsets = wnl.i4_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, l, bc)

def i4lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    diags,offsets = wnl.i4lapl_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, l, bc)

def i4lapl2(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral bilaplacian r^l P_n^{-1/2, l-1/2}(2r^2 - 1)."""

    diags,offsets = wnl.i4lapl2_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, l, bc)

def i4qm(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of Q r^{l-1} P_n^{-1/2, l-3/2}(2r^2 - 1)."""

    diags,offsets = wnl.i4qm_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, l, bc)

def i4qp(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of Q r^{l+1} P_n^{-1/2, l+1/2}(2r^2 - 1)."""

    diags,offsets = wnl.i4qp_diags(nr, l)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
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

    if make_square:
        bc['rb'] = bc['rt']
    bc['rt'] = 0

    return radbc.constrain(mat, l, bc)
