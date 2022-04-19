"""Module provides functions to generate sparse operators for the spherical harmonic expansion in a sphere."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils


def zblk(maxnl, m):
    """Create a block of zeros"""

    nl = maxnl - m
    mat = spsp.lil_matrix((nl,nl))
    return mat

def coriolisdr(maxnl, m, coeff = 1.0):
    """Create operator for the spherical harmonic coriolis Q term with D_r radial dependence."""

    ls = np.arange(m, maxnl)
    offsets = np.arange(-1,2,2)
    nzrow = -1

    def dgen(l):
        return -(l - 1.0)*(l + 1.0)*np.sqrt(((l - m)*(l + m))/((2.0*l - 1.0)*(2.0*l + 1.0)))

    # Generate 1st subdiagonal
    def d_1(l):
        return dgen(l)

    # Generate 1st superdiagonal
    def d1(l):
        return dgen(l + 1.0)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ls, nzrow, ds, offsets, None, False)

    mat = coeff*spsp.diags(diags, offsets)
    return mat

def coriolis_r(maxnl, m, coeff = 1.0):
    """Create operator for the spherical harmonic coriolis Q term with 1/r radial dependence."""

    ls = np.arange(m, maxnl)
    offsets = np.arange(-1,2,2)
    nzrow = -1

    def dgen(l):
        return  (l - 1.0)*(l + 1.0)*np.sqrt(((l - m)*(l + m))/((2.0*l - 1.0)*(2.0*l + 1.0)))

    # Generate 1st subdiagonal
    def d_1(l):
        return  (l - 1.0)*dgen(l)

    # Generate 1st superdiagonal
    def d1(l):
        return -(l + 2.0)*dgen(l + 1.0)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ls, nzrow, ds, offsets, None, False)

    mat = coeff*spsp.diags(diags, offsets)
    return mat

def coriolism(maxnl, m, coeff = 1.0):
    """Create operator for the spherical harmonic coriolis Q (l-1) term."""

    ls = np.arange(m, maxnl)
    offsets = np.array([-1])
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(l):
        return (l - 1.0)*(l + 1.0)*np.sqrt(((l - m)*(l + m))/((2.0*l - 1.0)*(2.0*l + 1.0)))

    ds = [d_1]
    diags, offsets = utils.build_diagonals(ls, nzrow, ds, offsets, None, False)

    mat = coeff*spsp.diags(diags, offsets)
    return mat

def coriolisp(maxnl, m, coeff = 1.0):
    """Create operator for the spherical harmonic coriolis Q (l+1) term."""

    ls = np.arange(m, maxnl)
    offsets = np.array([1])
    nzrow = -1

    # Generate 1st superdiagonal
    def d1(l):
        return  -l*(l + 2.0)*np.sqrt(((l - m + 1.0)*(l + m + 1.0))/((2.0*l + 1.0)*(2.0*l + 3.0)))

    ds = [d1]
    diags, offsets = utils.build_diagonals(ls, nzrow, ds, offsets, None, False)

    mat = coeff*spsp.diags(diags, offsets)
    return mat

def qid(maxnl, m, coeff = 1.0):
    """Create a quasi identity block"""

    nl = maxnl - m
    mat = coeff*spsp.eye(nl).tocoo()
    return mat
