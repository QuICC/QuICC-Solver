"""Module provides functions to generate values for the Worland expansion of Chebyshev type"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import scipy.special as special

import quicc.base.utils as utils
from quicc.geometry.worland.base import (eval_wnl, eval_pnab_origin,
        apply_norm, compute_norm_row, compute_norm_row_l_1)

def jacobi_alpha(l):
    """Jacobi polynomial alpha parameter"""

    return -0.5

def jacobi_beta(l):
    """Jacobi polynomial beta parameter"""

    return  l - 0.5

def get_grid(nr):
    """Physical space grid"""

    x = np.cos(np.pi*(np.arange(0,nr)+0.5)/(nr))
    return np.sqrt((x + 1.0)/2.0)

def get_weights(nr):
    """Gaussian integration weight"""

    return np.full(nr, np.pi/(2.0*nr))

def get_normln(n, l):
    """Natural log of unit norm"""

    assert(l >= 0)

    if l == 0:
        normln = -np.log(2.0) + special.gammaln(n + 0.5) - special.gammaln(n + 1.0)
        if (np.ndim(n) == 0 and n == 0):
            normln += 0.5*np.log(2.0)
        elif np.ndim(n) == 1 and n[0] == 0:
            normln[0] += 0.5*np.log(2.0)
    else:
        normln = -np.log(2.0*(2.0*n + l)) + special.gammaln(n + 0.5) + special.gammaln(n + l + 0.5) - special.gammaln(n + l) - special.gammaln(n + 1.0)
        normln *= 0.5

    return normln

def get_norm(n, l):
    """Unit norm"""

    norm = np.exp(get_normln(n, l))
    return norm

def get_invnorm(n, l):
    """Inverse unit norm"""

    invnorm = np.exp(-get_normln(n, l))
    return invnorm

def eval_poly(n, l, nr):

    return eval_wnl(n, l, -0.5, -0.5, nr, get_grid, get_invnorm)

def normalize(val, l):
    """Normalize array of values"""

    apply_norm(val, l, get_invnorm)

def normalize_row(n, l, k, p = 0):
    """Normalization factor for matrix row. l is from LHS, p allows to shift RHS to l+p"""

    norm = compute_norm_row(n, l, k, p, get_normln)

    return norm

def normalize_row_l_1(n, l, k):
    """Normalization factor for matrix row from W_n^l to Wn^{l-1}"""

    norm = compute_norm_row_l_1(n, l, k, get_normln)

    return norm

def eval_jacobi_origin(nr, l, k = 0, normalized = True):
    """Compute the value at origin for Worland polynomials"""

    return eval_pnab_origin(nr, l, -0.5, k, normalized, normalize)

def eval_bc_poly(nr, l, k = 0, normalized = True):
    """Compute the endpoint value for Worland polynomials"""

    val = np.zeros(nr)
    val[0] = 1.0
    if nr > 0:
        for i in range(1,nr):
            val[i] = val[i-1]*(2.0*i-1.0 + 2*k)/(2.0*i)

    # Normalize
    if normalized:
        normalize(val, l)

    return val

def r2_diags(nr, l):
    """Create operator for 1st integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        if l == 0: # Continuity
            val = normalize_row(n,l,-1)*n/(2.0*(l + 2.0*n - 1.0))
            val[0] = normalize_row(n[0:1],l,-1)*1.0/(l + 1.0)
        else:
            val = normalize_row(n,l,-1)*n*(l + n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return val

    # Generate main diagonal
    def d0(n):
        if l == 1:
            val = normalize_row(n,l,0)*(1.0 + n)/(2.0*(n + 1.0))
        else:
            val = normalize_row(n,l,0)*(2.0*l**2 + 4.0*l*n - l + 4.0*n**2 - 1.0)/(2*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))
        if l == 0 or l == 1: # Continuity
            val[0] = normalize_row(n[0:1],l,0)*(2.0*l + 1.0)/(2.0*(l + 1.0))
        return val

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)/(4.0*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1_diags(nr, l):
    """Create operator for 1st integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        if l == 0:
            val = normalize_row(n,l,-1)/((2.0*n - 1.0))
            val[0] = normalize_row(n[0:1],l,-1)*2.0/(l + 1.0)
        else:
            val = normalize_row(n,l,-1)*2.0*(l + n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return val

    # Generate main diagonal
    def d0(n):
        return -normalize_row(n,l,0)*2.0*l/((l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1)*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)/(2.0*(l + n)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qm_diags(nr, l):
    """Create operator for 1st integral of Q r^{l-1} P_n^{-1/2,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(0,2)
    nzrow = 0

    # Generate main diagonal
    def d0(n):
        if l == 0:
            val = -normalize_row(n,l,0,-1)*4.0*(n - 1.0)/(2.0*n - 1.0)
        else:
            val = -normalize_row(n,l,0,-1)*4.0*(l + n - 1.0)/(l + 2.0*n - 1.0)
        return val

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,-1)*2.0*(2.0*n + 1.0)/(l + 2.0*n + 1.0)

    ds = [d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qp_diags(nr, l):
    """Create operator for 1st integral of Q r^{l+1} P_n^{-1/2,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1,1)*2.0*(2.0*l + 2.0*n + 1.0)/(l + 2.0*n - 1.0)

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0,1)*(2.0*n - 1.0)*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n + 1.0))

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2_diags(nr, l):
    """Create operator for 2nd integral r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        if l == 0:
            val = normalize_row(n,l,-2)/((l + 2.0*n - 3.0)*(l + 2.0*n - 1.0))
            val[0] = normalize_row(n[0:1],l,-2)*4.0*(l + 1.0)/((l + 1.0)*(l + 2.0)*(l + 3.0))
        else:
            val = normalize_row(n,l,-2)*4.0*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return val

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1)*8.0*l*(l + n - 1.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*2.0*(2.0*l**2 - 4.0*l*n - 4.0*n**2 + 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*2.0*l*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/(4.0*(l + n)*(l + n + 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2lapl_diags(nr, l):
    """Create operator for 2nd integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*8.0*(l + n - 1.0)*(2.0*l + 2.0*n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0)*8.0*(4.0*l*n + 4.0*n**2 - 1.0)/((l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*2.0*(2.0*n + 1.0)**2*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qm_diags(nr, l):
    """Create operator for 2nd integral of Q r^{l-1} P_n^{-1/2,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,3)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1,-1)*8.0*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0,-1)*4.0*(l + n - 1.0)*(2.0*l - 2.0*n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1,-1)*2.0*(2.0*n + 1.0)*(4.0*l + 2.0*n - 1.0)/((l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2,-1)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    ds = [d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qp_diags(nr, l):
    """Create operator for 2nd integral of Q r^{l+1} P_n^{-1/2,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2,1)*4.0*(l + n - 1.0)*(2.0*l + 2.0*n - 1.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1,1)*2.0*(4.0*l**2 + 4.0*l - 4.0*n**2 + 4.0*n + 3.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -normalize_row(n,l,0,1)*(2.0*l + 2.0*n + 1.0)*(8.0*l*n + 4.0*n**2 + 4.0*n - 3.0)/((l + n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,1)*(2.0*n + 1.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/(2.0*(l + n)*(l + n + 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    ds = [d_2, d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4_diags(nr, l):
    """Create operator for 4th integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        if l == 0:
            val = normalize_row(n,l,-4)/((l + 2.0*n - 7.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 1.0))
            val[0] = normalize_row(n[0:1],l,-4)*16.0/((l + 4.0)*(l + 5.0)*(l + 6.0)*(l + 7.0))
        else:
            val = normalize_row(n,l,-4)*16.0*(l + n - 4.0)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return val

    # Generate 3rd subdiagonal
    def d_3(n):
        return -normalize_row(n,l,-3)*64.0*l*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2)*16.0*(l + n - 2.0)*(l + n - 1.0)*(6.0*l**2 - 4.0*l*n + 4.0*l - 4.0*n**2 + 8.0*n + 5.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1)*16.0*l*(l + n - 1.0)*(4.0*l**2 - 12.0*l*n + 6.0*l - 12.0*n**2 + 12.0*n + 17.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*2.0*(8.0*l**4 - 96.0*l**3*n - 48.0*l**2*n**2 + 100.0*l**2 + 96.0*l*n**3 - 120.0*l*n + 48.0*n**4 - 120.0*n**2 + 27.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*4.0*l*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)*(4.0*l**2 - 12.0*l*n - 6.0*l - 12.0*n**2 - 12.0*n + 17.0)/((l + n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(6.0*l**2 - 4.0*l*n - 4.0*l - 4.0*n**2 - 8.0*n + 5.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return normalize_row(n,l,3)*l*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/((l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))

    # Generate 4th superdiagonal
    def d4(n):
        return normalize_row(n,l,4)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)/(16.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + n + 3.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl_diags(nr, l):
    """Create operator for 4th integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return normalize_row(n,l,-3)*32.0*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -normalize_row(n,l,-2)*32.0*(l + n - 2.0)*(l + n - 1.0)*(4.0*l**2 - 6.0*l - 4.0*n**2 + 8.0*n + 5.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*8.0*(l + n - 1.0)*(8.0*l**3 - 40.0*l**2*n - 4.0*l**2 - 56.0*l*n**2 + 80.0*l*n + 62.0*l - 8.0*n**3 + 36.0*n**2 - 22.0*n - 21.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0)*16.0*(16.0*l**3*n + 4.0*l**3 - 24.0*l**2*n - 24.0*l**2 - 32.0*l*n**3 - 24.0*l*n**2 + 40.0*l*n + 14.0*l - 16.0*n**4 + 40.0*n**2 - 9.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*2.0*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)*(48.0*l**2*n + 48.0*l**2 + 32.0*l*n**2 + 8.0*l*n - 84.0*l - 8.0*n**3 - 36.0*n**2 - 22.0*n + 21.0)/((l + n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*2.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(8.0*l*n + 14.0*l + 4.0*n**2 + 8.0*n - 5.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return normalize_row(n,l,3)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/(2.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl2_diags(nr, l):
    """Create operator for 4th integral bilaplacian r^l P_n^{-1/2, l-1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2)*64.0*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)*(2.0*l + 2.0*n - 3.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*128.0*(l + n - 1.0)*(2.0*l + 2.0*n - 3.0)*(4.0*l*n + 4.0*n**2 - 4.0*n - 3.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*96.0*(2.0*n + 1.0)*(2.0*l + 2.0*n - 1.0)*(4.0*l*n + 2.0*l + 4.0*n**2 - 5.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*32.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(4.0*l*n + 4.0*l + 4.0*n**2 + 4.0*n - 3.0)/((l + n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*4.0*(2.0*n + 1.0)*(2.0*n + 3.0)**2*(2.0*n + 5.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qm_diags(nr, l):
    """Create operator for 4th integral of Q r^{l-1} P_n^{-1/2, l-3/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,5)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return - normalize_row(n,l,-3,-1)*32.0*(l + n - 4.0)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return  normalize_row(n,l,-2,-1)*16.0*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)*(6.0*l - 2.0*n - 1.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return - normalize_row(n,l,-1,-1)*24.0*(l + n - 2.0)*(l + n - 1.0)*(4.0*l**2 - 8.0*l*n - 4.0*n**2 + 8.0*n + 5.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate diagonal
    def d0(n):
        return  normalize_row(n,l,0,-1)*4.0*(l + n - 1.0)*(8.0*l**3 - 72.0*l**2*n - 12.0*l**2 - 24.0*l*n**2 + 72.0*l*n + 58.0*l + 24.0*n**3 + 12.0*n**2 - 54.0*n - 27.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return  normalize_row(n,l,1,-1)*2.0*(2.0*n + 1.0)*(32.0*l**3 - 48.0*l**2*n - 72.0*l**2 - 96.0*l*n**2 - 48.0*l*n + 112.0*l - 24.0*n**3 + 12.0*n**2 + 54.0*n - 27.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return  normalize_row(n,l,2,-1)*3.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(8.0*l**2 - 8.0*l - 4.0*n**2 - 8.0*n + 5.0)/((l + n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return  normalize_row(n,l,3,-1)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(8.0*l + 2.0*n - 1.0)/(2.0*(l + n)*(l + n + 1.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    # Generate 4th superdiagonal
    def d4(n):
        return  normalize_row(n,l,4,-1)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/(4.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qp_diags(nr, l):
    """Create operator for 4th integral of Q r^{l+1} P_n^{-1/2, l+1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return normalize_row(n,l,-4,1)*16.0*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)/((l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -normalize_row(n,l,-3,1)*8.0*(l + n - 2.0)*(l + n - 1.0)*(12.0*l**2 + 8.0*l*n - 16.0*l - 4.0*n**2 + 12.0*n + 7.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2,1)*12.0*(l + n - 1.0)*(8.0*l**3 - 8.0*l**2*n + 4.0*l**2 - 24.0*l*n**2 + 40.0*l*n + 26.0*l - 8.0*n**3 + 20.0*n**2 + 2.0*n - 5.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1,1)*2.0*(16.0*l**4 - 128.0*l**3*n + 32.0*l**3 - 192.0*l**2*n**2 + 96.0*l**2*n + 272.0*l**2 + 32.0*l*n + 112.0*l + 48.0*n**4 - 48.0*n**3 - 144.0*n**2 + 108.0*n + 81.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return -normalize_row(n,l,0,1)*(2.0*l + 2.0*n + 1.0)*(64.0*l**3*n + 16.0*l**3 - 96.0*l**2*n**2 - 48.0*l**2*n - 96.0*l**2 - 192.0*l*n**3 - 144.0*l*n**2 + 320.0*l*n - 4.0*l - 48.0*n**4 - 48.0*n**3 + 144.0*n**2 + 108.0*n - 81.0)/((l + n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,1)*3.0*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(16.0*l**2*n + 16.0*l**2 - 24.0*l - 8.0*n**3 - 20.0*n**2 + 2.0*n + 5.0)/(2.0*(l + n)*(l + n + 1.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -normalize_row(n,l,2,1)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(16.0*l*n + 28.0*l + 4.0*n**2 + 12.0*n - 7.0)/(4.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -normalize_row(n,l,3,1)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)/(8.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + n + 3.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i6_diags(nr, m):
    """Create operator for 6th integral r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-6,7)
    nzrow = 5

    # Generate 6th subdiagonal
    def d_6(n):
        if m == 0:
            val = normalize_row(n, m, -6)/((m + 2.0*n - 11.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 1.0))
            val[0] = normalize_row(n[0:1], m, -6)*2.0/10395.0
        else:
            val = normalize_row(n, m, -6)*64.0*(m + n - 6.0)*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 12.0)*(m + 2.0*n - 11.0)*(m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))
        return val

    # Generate 5th subdiagonal
    def d_5(n):
        return -normalize_row(n, m, -5)*384.0*m*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 11.0)*(m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return normalize_row(n, m, -4)*96.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(10.0*m**2 - 4.0*m*n + 8.0*m - 4.0*n**2 + 16.0*n + 9.0)/((m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -normalize_row(n, m, -3)*160.0*m*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(8.0*m**2 - 12.0*m*n + 18.0*m - 12.0*n**2 + 36.0*n + 37.0)/((m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n, m, -2)*60.0*(m + n - 2.0)*(m + n - 1.0)*(16.0*m**4 - 64.0*m**3*n + 64.0*m**3 - 48.0*m**2*n**2 + 96.0*m**2*n + 236.0*m**2 + 32.0*m*n**3 - 96.0*m*n**2 - 40.0*m*n + 104.0*m + 16.0*n**4 - 64.0*n**3 - 40.0*n**2 + 208.0*n + 105.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n, m, -1)*48.0*m*(m + n - 1.0)*(8.0*m**4 - 80.0*m**3*n + 40.0*m**3 + 280.0*m**2 + 160.0*m*n**3 - 240.0*m*n**2 - 440.0*m*n + 260.0*m + 80.0*n**4 - 160.0*n**3 - 440.0*n**2 + 520.0*n + 537.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n, m, 0)*4.0*(16.0*m**6 - 480.0*m**5*n + 960.0*m**4*n**2 + 1120.0*m**4 + 2560.0*m**3*n**3 - 7840.0*m**3*n + 480.0*m**2*n**4 - 5040.0*m**2*n**2 + 5614.0*m**2 - 960.0*m*n**5 + 5600.0*m*n**3 - 5180.0*m*n - 320.0*n**6 + 2800.0*n**4 - 5180.0*n**2 + 1125.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n, m, 1)*12.0*m*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(8.0*m**4 - 80.0*m**3*n - 40.0*m**3 + 280.0*m**2 + 160.0*m*n**3 + 240.0*m*n**2 - 440.0*m*n - 260.0*m + 80.0*n**4 + 160.0*n**3 - 440.0*n**2 - 520.0*n + 537.0)/((m + n)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n, m, 2)*15.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(16.0*m**4 - 64.0*m**3*n - 64.0*m**3 - 48.0*m**2*n**2 - 96.0*m**2*n + 236.0*m**2 + 32.0*m*n**3 + 96.0*m*n**2 - 40.0*m*n - 104.0*m + 16.0*n**4 + 64.0*n**3 - 40.0*n**2 - 208.0*n + 105.0)/(4.0*(m + n)*(m + n + 1.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return normalize_row(n, m, 3)*5.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(8.0*m**2 - 12.0*m*n - 18.0*m - 12.0*n**2 - 36.0*n + 37.0)/(2.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0))

    # Generate 4th superdiagonal
    def d4(n):
        return normalize_row(n, m, 4)*3.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(10.0*m**2 - 4.0*m*n - 8.0*m - 4.0*n**2 - 16.0*n + 9.0)/(8.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0))

    # Generate 5rd superdiagonal
    def d5(n):
        return normalize_row(n, m, 5)*3.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)/(8.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0)*(m + 2.0*n + 11.0))

    # Generate 6th superdiagonal
    def d6(n):
        return normalize_row(n, m, 6)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*n + 11.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)*(2.0*m + 2.0*n + 11.0)/(64.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + n + 5.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0)*(m + 2.0*n + 11.0)*(m + 2.0*n + 12.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def stencil_value_diags(nr, l):
    """Create stencil matrix for a zero boundary value"""

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1)*2.0*n/(2.0*n - 1.0)

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_diff_diags(nr, l):
    """Create stencil matrix for a zero 1st derivative"""

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n,l,-1)*2.0*n*(4.0*(-1.0 + n)**2 + l*(-3.0 + 4.0*n))
        den = (-1.0 + 2.0*n)*(l + 4.0*l*n + 4.0*n**2)
        return num/den

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_rdiffdivr_diags(nr, l):
    """Create stencil matrix for a zero r D 1/r"""

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n,l,-1)*2.0*n*(3.0 - 3.0*l - 8.0*n + 4.0*l*n + 4.0*n**2)
        den = (-1.0 + 2.0*n)*(-1.0 + l + 4.0*l*n + 4.0*n**2)
        if l == 1:
            num[0] = 0.0
            den[0] = 1.0
        return num/den

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_insulating_sph_diags(nr, l):
    """Create stencil matrix for a insulating boundary"""

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n,l,-1)*2.0*n*(4.0*(-1.0 + n)**2 + l*(-2.0 + 4.0*n) + 1.0)
        den = (-1.0 + 2.0*n)*(2.0*l + 1.0 + 4.0*l*n + 4.0*n**2)
        return num/den

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_value_diff_diags(nr, l):
    """Create stencil matrix for a zero boundary value and zero 1st derivative"""

    ns = np.arange(0,nr)
    offsets = [-2, -1, 0]

    # Generate second subdiagonal
    def d_2(n):
        num = normalize_row(n,l,-2)*4.0*n*(n - 1.0)*(-3.0 + l + 2.0*n)
        den = (-1.0 + l + 2.0*n)*(3.0 + 4.0*(n - 2.0)*n)
        return num/den

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n,l,-1)*4.0*n*(l + 2.0*n)
        den = (-1.0 + 2.0*n)*(1.0 + l + 2.0*n)
        return num/den

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_value_diff2_diags(nr, l):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    ns = np.arange(0,nr)
    offsets = [-2, -1, 0]

    # Generate 2nd subdiagonal
    def d_2(n):
        num = normalize_row(n,l,-2)*4.0*(-1.0 + n)*n*(-3.0 + l + 2.0*n)*(19.0 + 8.0*(-3.0 + n)*n + 2.0*l*(-5.0 + 4.0*n))
        den = (-1.0 + l + 2.0*n)*(3.0 + 4.0*(-2.0 + n)*n)*(3.0 + 8.0*(-1.0 + n)*n + l*(-2.0 + 8.0*n))
        return num/den

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n,l,-1)*4.0*n*(l + 2.0*n)*(7.0 + 8.0*n**2 + l*(2.0 + 8.0*n))
        den = (-1.0 + 2.0*n)*(1.0 + l + 2.0*n)*(3.0 + 8.0*n*(1.0 + n) + l*(6.0 + 8.0*n))
        return num/den

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)
