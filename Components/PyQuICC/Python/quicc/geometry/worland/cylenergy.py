"""Module provides functions to generate values for the Worland expansion of cylindrical energy orthogonality type"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.polynomial.legendre as leg
import scipy.sparse as spsp
import scipy.special as special

import quicc.base.utils as utils
from quicc.geometry.worland.base import (eval_wnl, eval_pnab_origin,
        apply_norm, compute_norm_row, compute_norm_row_l_1)


def jacobi_alpha(l):
    """Jacobi polynomial alpha parameter"""

    return 0.0

def jacobi_beta(l):
    """Jacobi polynomial beta parameter"""

    return l

def get_grid(nr):
    """Physical space grid"""

    x, tmp = leg.leggauss(nr)
    return np.sqrt((x + 1.0)/2.0)

def get_weights(nr):
    """Gaussian integration weight"""

    t, w = leg.leggauss(nr)

    return w/4.0

def get_normln(n, l):
    """Natural log of norm"""

    normln = -0.5*np.log(2.0*(2.0*n+l+1.0))

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

    return eval_wnl(n, l, 0.0, 0.0, nr, get_grid, get_invnorm)

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

    return eval_pnab_origin(nr, l, 0.0, k, normalized, normalize)

def eval_bc_poly(nr, l, k = 0, normalized = True):
    """Compute the endpoint value for Worland polynomials"""

    if k == 0:
        val = np.ones(nr)
    else:
        val = np.zeros(nr)
        val[0] = 1.0
        if nr > 0:
            for i in range(1,nr):
                val[i] = val[i-1]*(i + k)/i

    # Normalize
    if normalized:
        normalize(val, l)

    return val

def r2_diags(nr, l):
    """Create operator for 1st integral r^l P_n^{0,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 1. subdiagonal
    def d_1(n):
        if l == 0:
            val = n/(2.0*(2.0*n - 1.0))
        else:
            val = n*(l + n)/((l + 2.0*n)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        if l == 0:
            val = (n + 1.0)/(2.0*n + 2.0)
        else:
            val = (l**2 + 2.0*l*n + l + 2.0*n**2 + 2.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = (n + 1.0)*(l + n + 1.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 1)*val

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1_diags(nr, l):
    """Create operator for 1st integral r^l P_n^{0,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1. subdiagonal
    def d_1(n):
        val = 2.0*(l + n)/((l + 2.0*n)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = -2.0*l/((l + 2.0*n)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = -2.0*(n + 1.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 1)*val

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qm_diags(nr, l):
    """Create operator for 1st integral of Q r^{l-1} P_n^{0,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(0,2)
    nzrow = 0

    # Generate main diagonal
    def d0(n):
        val = -4.0*(l + n)/(l + 2.0*n)
        return normalize_row(n, l, 0, -1)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = -4.0*(n + 1.0)/(l + 2.0*n + 2.0)
        return normalize_row(n, l, 1, -1)*val

    ds = [d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qp_diags(nr, l):
    """Create operator for 1st integral of Q r^{l+1} P_n^{0,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1. subdiagonal
    def d_1(n):
        val = 2.0*(2.0*l + 2.0*n + 1.0)/(l + 2.0*n)
        return normalize_row(n, l, -1, 1)*val

    # Generate main diagonal
    def d0(n):
        val = 2.0*(2.0*n + 1.0)/(l + 2.0*n + 2.0)
        return normalize_row(n, l, 0, 1)*val

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2_diags(nr, l):
    """Create operator for 2nd integral r^l P_n^{0,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2. subdiagonal
    def d_2(n):
        val = 4.0*(l + n)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -2)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -8.0*l*(l + n)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 4.0*(l**2 - 2.0*l*n - l - 2.0*n**2 - 2.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 8.0*l*(n + 1.0)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 4.0*(n + 1.0)*(n + 2.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 2)*val

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2lapl_diags(nr, l):
    """Create operator for 2nd integral of Laplacian r^l P_n^{0,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1. subdiagonal
    def d_1(n):
        val = 8.0*(l + n)*(2.0*l + 2.0*n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 8.0*(4.0*l*n + 3.0*l + 4.0*n**2 + 4.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 8.0*(n + 1.0)*(2.0*n + 3.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 1)*val

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qm_diags(nr, l):
    """Create operator for 2nd integral of Q r^{l-1} P_n^{0,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,3)
    nzrow = 1

    # Generate 1. subdiagonal
    def d_1(n):
        val = -8.0*(l + n)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -1, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 8.0*(l + n)*(l - n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, 0, -1)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 8.0*(2.0*l + n)*(n + 1.0)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 1, -1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 8.0*(n + 1.0)*(n + 2.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 2, -1)*val

    ds = [d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qp_diags(nr, l):
    """Create operator for 2nd integral of Q r^{l+1} P_n^{0,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2. subdiagonal
    def d_2(n):
        val = 4.0*(l + n)*(2.0*l + 2.0*n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -2, 1)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -4.0*(2.0*l**2 - 2.0*n**2 - n + 1.0)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -1, 1)*val

    # Generate main diagonal
    def d0(n):
        val = -4.0*(4.0*l*n + 3.0*l + 2.0*n**2 + 3.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 0, 1)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = -4.0*(n + 1.0)*(2.0*n + 3.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 1, 1)*val

    ds = [d_2, d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4_diags(nr, l):
    """Create operator for 4th integral r^l P_n^{0,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4. subdiagonal
    def d_4(n):
        val = 16.0*(l + n)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -4)*val

    # Generate 3. subdiagonal
    def d_3(n):
        val = -64.0*l*(l + n)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -3)*val

    # Generate 2. subdiagonal
    def d_2(n):
        val = 32.0*(l + n)*(l + n - 1.0)*(3.0*l**2 - 2.0*l*n + l - 2.0*n**2 + 2.0*n + 4.0)/((l + 2.0*n)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, -2)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -64.0*l*(l + n)*(l**2 - 3.0*l*n - 3.0*n**2 + 5.0)/((l + 2.0*n)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 16.0*(l**4 - 12.0*l**3*n - 6.0*l**3 - 6.0*l**2*n**2 - 6.0*l**2*n + 11.0*l**2 + 12.0*l*n**3 + 18.0*l*n**2 - 6.0*l*n - 6.0*l + 6.0*n**4 + 12.0*n**3 - 6.0*n**2 - 12.0*n)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 64.0*l*(n + 1.0)*(l**2 - 3.0*l*n - 3.0*l - 3.0*n**2 - 6.0*n + 2.0)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))
        return normalize_row(n, l, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 32.0*(n + 1.0)*(n + 2.0)*(3.0*l**2 - 2.0*l*n - 3.0*l - 2.0*n**2 - 6.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))
        return normalize_row(n, l, 2)*val

    # Generate 3. superdiagonal
    def d3(n):
        val = 64.0*l*(n + 1.0)*(n + 2.0)*(n + 3.0)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))
        return normalize_row(n, l, 3)*val

    # Generate 4. superdiagonal
    def d4(n):
        val = 16.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0))
        return normalize_row(n, l, 4)*val

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl_diags(nr, l):
    """Create operator for 4th integral of Laplacian r^l P_n^{0,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3. subdiagonal
    def d_3(n):
        val = 32.0*(l + n)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)/((l + 2.0*n)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -3)*val

    # Generate 2. subdiagonal
    def d_2(n):
        val = -32.0*(l + n)*(l + n - 1.0)*(4.0*l**2 - 9.0*l - 4.0*n**2 + 4.0*n + 8.0)/((l + 2.0*n)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -2)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = 32.0*(l + n)*(2.0*l**3 - 10.0*l**2*n - 9.0*l**2 - 14.0*l*n**2 + 9.0*l*n + 22.0*l - 2.0*n**3 + 9.0*n**2 + 2.0*n - 9.0)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 32.0*(8.0*l**3*n + 7.0*l**3 - 18.0*l**2*n - 21.0*l**2 - 16.0*l*n**3 - 42.0*l*n**2 - 10.0*l*n + 14.0*l - 8.0*n**4 - 16.0*n**3 + 8.0*n**2 + 16.0*n)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 32.0*(n + 1.0)*(12.0*l**2*n + 21.0*l**2 + 8.0*l*n**2 + 7.0*l*n - 21.0*l - 2.0*n**3 - 15.0*n**2 - 22.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 32.0*(n + 1.0)*(n + 2.0)*(8.0*l*n + 21.0*l + 4.0*n**2 + 12.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))
        return normalize_row(n, l, 2)*val

    # Generate 3. superdiagonal
    def d3(n):
        val = 32.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(2.0*n + 7.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))
        return normalize_row(n, l, 3)*val

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl2_diags(nr, l):
    """Create operator for 4th integral bilaplacian r^l P_n^{0, l-1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2. subdiagonal
    def d_2(n):
        val = 64.0*(l + n)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)*(2.0*l + 2.0*n - 3.0)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -2)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = 128.0*(l + n)*(2.0*l + 2.0*n - 3.0)*(4.0*l*n + 3.0*l + 4.0*n**2 - 4.0)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 64.0*(24.0*l**2*n**2 + 60.0*l**2*n + 35.0*l**2 + 48.0*l*n**3 + 108.0*l*n**2 + 26.0*l*n - 35.0*l + 24.0*n**4 + 48.0*n**3 - 10.0*n**2 - 34.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 128.0*(n + 1.0)*(2.0*n + 5.0)*(4.0*l*n + 7.0*l + 4.0*n**2 + 8.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 64.0*(n + 1.0)*(n + 2.0)*(2.0*n + 5.0)*(2.0*n + 7.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 2)*val

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qm_diags(nr, l):
    """Create operator for 4th integral of Q r^{l-1} P_n^{0, l-3/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,5)
    nzrow = 3

    # Generate 3. subdiagonal
    def d_3(n):
        val = -32.0*(l + n)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -3, -1)*val

    # Generate 2. subdiagonal
    def d_2(n):
        val = 32.0*(l + n)*(l + n - 2.0)*(l + n - 1.0)*(3.0*l - n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -2, -1)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -96.0*(l + n)*(l + n - 1.0)*(l**2 - 2.0*l*n - l - n**2 + n + 2.0)/((l + 2.0*n)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, -1, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 32.0*(l + n)*(l**3 - 9.0*l**2*n - 6.0*l**2 - 3.0*l*n**2 + 6.0*l*n + 11.0*l + 3.0*n**3 + 6.0*n**2 - 3.0*n - 6.0)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, 0, -1)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 32.0*(n + 1.0)*(4.0*l**3 - 6.0*l**2*n - 12.0*l**2 - 12.0*l*n**2 - 18.0*l*n + 8.0*l - 3.0*n**3 - 3.0*n**2 + 6.0*n)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 1, -1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 96.0*(n + 1.0)*(n + 2.0)*(2.0*l**2 - 2.0*l - n**2 - 3.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))
        return normalize_row(n, l, 2, -1)*val

    # Generate 3. superdiagonal
    def d3(n):
        val = 32.0*(4.0*l + n)*(n + 1.0)*(n + 2.0)*(n + 3.0)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))
        return normalize_row(n, l, 3, -1)*val

    # Generate 4. superdiagonal
    def d4(n):
        val = 32.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))
        return normalize_row(n, l, 4, -1)*val

    ds = [d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qp_diags(nr, l):
    """Create operator for 4th integral of Q r^{l+1} P_n^{0, l+1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4. subdiagonal
    def d_4(n):
        val = 16.0*(l + n)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)/((l + 2.0*n)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -4, 1)*val

    # Generate 3. subdiagonal
    def d_3(n):
        val = -16.0*(l + n)*(l + n - 1.0)*(6.0*l**2 + 4.0*l*n - 10.0*l - 2.0*n**2 + 3.0*n + 5.0)/((l + 2.0*n)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -3, 1)*val

    # Generate 2. subdiagonal
    def d_2(n):
        val = 48.0*(l + n)*(2.0*l**3 - 2.0*l**2*n - 2.0*l**2 - 6.0*l*n**2 + 4.0*l*n + 8.0*l - 2.0*n**3 + 3.0*n**2 + 3.0*n - 2.0)/((l + 2.0*n)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, -2, 1)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -16.0*(2.0*l**4 - 16.0*l**3*n - 8.0*l**3 - 24.0*l**2*n**2 - 6.0*l**2*n + 28.0*l**2 + 12.0*l*n**2 + 10.0*l*n + 2.0*l + 6.0*n**4 + 9.0*n**3 - 12.0*n**2 - 9.0*n + 6.0)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, -1, 1)*val

    # Generate main diagonal
    def d0(n):
        val = -16.0*(8.0*l**3*n + 7.0*l**3 - 12.0*l**2*n**2 - 27.0*l**2*n - 21.0*l**2 - 24.0*l*n**3 - 57.0*l*n**2 - 8.0*l*n + 14.0*l - 6.0*n**4 - 15.0*n**3 + 3.0*n**2 + 18.0*n)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, 0, 1)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = -48.0*(n + 1.0)*(4.0*l**2*n + 7.0*l**2 - 2.0*l*n - 7.0*l - 2.0*n**3 - 9.0*n**2 - 9.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))
        return normalize_row(n, l, 1, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = -16.0*(n + 1.0)*(n + 2.0)*(8.0*l*n + 21.0*l + 2.0*n**2 + 7.0*n)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))
        return normalize_row(n, l, 2, 1)*val

    # Generate 3. superdiagonal
    def d3(n):
        val = -16.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(2.0*n + 7.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))
        return normalize_row(n, l, 3, 1)*val

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i6_diags(nr, l):
    """Create operator for 6th integral r^l P_n^{0,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-6,7)
    nzrow = 5

    # Generate 6. subdiagonal
    def d_6(n):
        val = 64.0*(l + n)*(l + n - 5.0)*(l + n - 4.0)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 11.0)*(l + 2.0*n - 10.0)*(l + 2.0*n - 9.0)*(l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))
        return normalize_row(n, l, -6)*val

    # Generate 5. subdiagonal
    def d_5(n):
        val = -384.0*l*(l + n)*(l + n - 4.0)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n)*(l + 2.0*n - 10.0)*(l + 2.0*n - 9.0)*(l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0))
        return normalize_row(n, l, -5)*val

    # Generate 4. subdiagonal
    def d_4(n):
        val = 192.0*(l + n)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)*(5.0*l**2 - 2.0*l*n + 3.0*l - 2.0*n**2 + 6.0*n + 8.0)/((l + 2.0*n)*(l + 2.0*n - 9.0)*(l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))
        return normalize_row(n, l, -4)*val

    # Generate 3. subdiagonal
    def d_3(n):
        val = -640.0*l*(l + n)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l**2 - 3.0*l*n + 3.0*l - 3.0*n**2 + 6.0*n + 13.0)/((l + 2.0*n)*(l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))
        return normalize_row(n, l, -3)*val

    # Generate 2. subdiagonal
    def d_2(n):
        val = 960.0*(l + n)*(l + n - 1.0)*(l**4 - 4.0*l**3*n + 2.0*l**3 - 3.0*l**2*n**2 + 3.0*l**2*n + 17.0*l**2 + 2.0*l*n**3 - 3.0*l*n**2 - 7.0*l*n + 4.0*l + n**4 - 2.0*n**3 - 7.0*n**2 + 8.0*n + 12.0)/((l + 2.0*n)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))
        return normalize_row(n, l, -2)*val

    # Generate 1. subdiagonal
    def d_1(n):
        val = -384.0*l*(l + n)*(l**4 - 10.0*l**3*n + 35.0*l**2 + 20.0*l*n**3 - 70.0*l*n + 10.0*n**4 - 70.0*n**2 + 84.0)/((l + 2.0*n)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))
        return normalize_row(n, l, -1)*val

    # Generate main diagonal
    def d0(n):
        val = 64.0*(l**6 - 30.0*l**5*n - 15.0*l**5 + 60.0*l**4*n**2 + 60.0*l**4*n + 85.0*l**4 + 160.0*l**3*n**3 + 240.0*l**3*n**2 - 370.0*l**3*n - 225.0*l**3 + 30.0*l**2*n**4 + 60.0*l**2*n**3 - 270.0*l**2*n**2 - 300.0*l**2*n + 274.0*l**2 - 60.0*l*n**5 - 150.0*l*n**4 + 200.0*l*n**3 + 450.0*l*n**2 - 80.0*l*n - 120.0*l - 20.0*n**6 - 60.0*n**5 + 100.0*n**4 + 300.0*n**3 - 80.0*n**2 - 240.0*n)/((l + 2.0*n)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))
        return normalize_row(n, l, 0)*val

    # Generate 1. superdiagonal
    def d1(n):
        val = 384.0*l*(n + 1.0)*(l**4 - 10.0*l**3*n - 10.0*l**3 + 35.0*l**2 + 20.0*l*n**3 + 60.0*l*n**2 - 10.0*l*n - 50.0*l + 10.0*n**4 + 40.0*n**3 - 10.0*n**2 - 100.0*n + 24.0)/((l + 2.0*n)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))
        return normalize_row(n, l, 1)*val

    # Generate 2. superdiagonal
    def d2(n):
        val = 960.0*(n + 1.0)*(n + 2.0)*(l**4 - 4.0*l**3*n - 6.0*l**3 - 3.0*l**2*n**2 - 9.0*l**2*n + 11.0*l**2 + 2.0*l*n**3 + 9.0*l*n**2 + 5.0*l*n - 6.0*l + n**4 + 6.0*n**3 + 5.0*n**2 - 12.0*n)/((l + 2.0*n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0))
        return normalize_row(n, l, 2)*val

    # Generate 3. superdiagonal
    def d3(n):
        val = 640.0*l*(n + 1.0)*(n + 2.0)*(n + 3.0)*(2.0*l**2 - 3.0*l*n - 6.0*l - 3.0*n**2 - 12.0*n + 4.0)/((l + 2.0*n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0)*(l + 2.0*n + 10.0))
        return normalize_row(n, l, 3)*val

    # Generate 4. superdiagonal
    def d4(n):
        val = 192.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(5.0*l**2 - 2.0*l*n - 5.0*l - 2.0*n**2 - 10.0*n)/((l + 2.0*n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0)*(l + 2.0*n + 10.0)*(l + 2.0*n + 11.0))
        return normalize_row(n, l, 4)*val

    # Generate 5. superdiagonal
    def d5(n):
        val = 384.0*l*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)/((l + 2.0*n)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0)*(l + 2.0*n + 10.0)*(l + 2.0*n + 11.0)*(l + 2.0*n + 12.0))
        return normalize_row(n, l, 5)*val

    # Generate 6. superdiagonal
    def d6(n):
        val = 64.0*(n + 1.0)*(n + 2.0)*(n + 3.0)*(n + 4.0)*(n + 5.0)*(n + 6.0)/((l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0)*(l + 2.0*n + 9.0)*(l + 2.0*n + 10.0)*(l + 2.0*n + 11.0)*(l + 2.0*n + 12.0)*(l + 2.0*n + 13.0))
        return normalize_row(n, l, 6)*val

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def stencil_value_diags(nr, l):
    """Create stencil matrix for a zero boundary value"""

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -normalize_row(n, l, -1)*np.ones(n.shape)

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

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
        num = -normalize_row(n, l, -1)*(2.0*(n - 1.0)*(l + n - 1.0) + l)
        den = 2.0*n*(l + n) + l
        return num/den

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

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
        num = -normalize_row(n, l, -1.0)*(-1.0 + l - 2.0*l*n - 2.0*(-2.0 + n)*n)
        den = -1.0 + l + 2.0*l*n + 2.0*n**2
        if l == 1:
            num[0] = 1.0
            den[0] = 1.0
        return num/den

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

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
        num = -normalize_row(n, l, -1)*(2.0*n*(l + n - 2.0) + 3.0)
        den = 2.0*(l*n + l + n**2) + 1.0
        return num/den

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

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
        num = normalize_row(n, l, -2)*2.0*n*(l*(6.0*n - 7.0)+2.0*n*(5.0*n - 14.0) + 20.0)
        den = (2.0*n - 3.0)*(l*(6.0*n - 1.0) + 2.0*n*(5.0*n - 4.0) + 2.0)
        return num/den

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n, l, -1)*(l*(2.0*n*(12.0*n + 7.0) - 1.0) + 2.0*n*(n*(20.0*n + 9.0) + 2.0) + 2.0)
        den = (2.0*n - 1.0)*(l*(6.0*n + 5.0) + 2.0*n*(5.0*n + 6.0) + 4.0)
        return num/den

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)

def stencil_value_diff2_diags(nr, l):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    ns = np.arange(0,nr)
    offsets = [-2, -1, 0]

    # Generate 2. subdiagonal
    def d_2(n):
        num = normalize_row(n, l, -2)*2.0*n*(92.0 + l**2*(13.0 + 2.0*n*(-11.0 + 5.0*n)) + l*(-3.0 + 2.0*n)*(23.0 + 2.0*n*(-17.0 + 7.0*n)) + 2.0*n*(-118.0 + n*(116.0 + n*(-52.0 + 9.0*n))))
        den = (-3.0 + 2.0*n)*(2.0 - 3.0*l + l**2 - 2.0*(6.0 + (-6.0 + l)*l)*n + 2.0*(14.0 + l*(-13.0 + 5.0*l))*n**2 + 4.0*(-8.0 + 7.0*l)*n**3 + 18.0*n**4)
        return num/den

    # Generate subdiagonal
    def d_1(n):
        num = -normalize_row(n, l, -1)*(2.0 + l**2*(1.0 + 2.0*n*(9.0 + n*(21.0 + 20.0*n))) + l*(1.0 + 2.0*n)*(-3.0 + 2.0*n*(17.0 + n*(9.0 + 28.0*n))) + 2.0*n**2*(22.0 + n*(52.0 + n*(17.0 + 36.0*n))))
        den = (-1.0 + 2.0*n)*(4.0 + l**2*(9.0 + 2.0*n*(9.0 + 5.0*n)) + l*(1.0 + 2.0*n)*(11.0 + 2.0*n*(11.0 + 7.0*n)) + 2.0*n*(10.0 + n*(20.0 + n*(20.0 + 9.0*n))))
        return num/den

    # Generate main diagonal
    def d0(n):
        return normalize_row(n, l, 0)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]
    return (diags,offsets)
