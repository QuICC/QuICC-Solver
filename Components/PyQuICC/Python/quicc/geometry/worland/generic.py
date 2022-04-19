"""Module provides functions to generate values for the Worland expansion of generic type"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import scipy.special as special

import quicc.base.utils as utils
from quicc.geometry.worland.base import (eval_wnl, apply_norm, compute_norm_row, compute_norm_row_l_1)


def get_grid(nr):
    """Physical space grid"""

    raise NotImplementedError("Generic Worland quadrature grid not implemented yet")

def get_weights(nr):
    """Gaussian integration weight"""

    raise NotImplementedError("Generic Worland quadrature weights not implemented yet")

def get_normln(n, l, a):
    """Natural log of unit norm"""

    if a == -0.5 and l == 0:
        normln = -np.log(2.0) + special.gammaln(n + 0.5) - special.gammaln(n + 1.0)
        if (np.ndim(n) == 0 and n == 0):
            normln += 0.5*np.log(2.0)
        elif np.ndim(n) == 1 and n[0] == 0:
            normln[0] += 0.5*np.log(2.0)
    else:
        b = l - 0.5
        normln = -np.log(2.0*(2.0*n+a+b+1.0)) + special.gammaln(n+a+1.0) + special.gammaln(n+b+1.0) - special.gammaln(n+a+b+1.0) - special.gammaln(n+1.0)
        normln *= 0.5

    return normln

def get_norm(n, l, a):
    """Unit norm"""

    norm = np.exp(get_normln(n, l, a))
    return norm

def get_invnorm(n, l, a):
    """Inverse unit norm"""

    invnorm = np.exp(-get_normln(n, l, a))
    return invnorm

def eval_poly(n, l, a, db, nr):

    return eval_wnl(n, l, a, db, nr, get_grid, get_invnorm)

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

def eval_bc_poly(nr, l, a, k = 0, normalized = True):
    """Compute the endpoint value for Worland polynomials"""

    val = np.zeros(nr)
    val[0] = 1.0
    if nr > 0:
        for i in range(1,nr):
            val[i] = val[i-1]*(i + a + k)/i

    # Normalize
    if normalized:
        normalize(val, l)

    return val

def r2_diags(nr, l, a):
    """Create operator for 1st integral r^l P_n^{a,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        val = normalize_row(n,l,-1)*2.0*n*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))
        return val

    # Generate main diagonal
    def d0(n):
        val = normalize_row(n,l,0)*(4.0*a*l + 8.0*a*n + 2.0*a + 4.0*l**2 + 8.0*l*n + 8.0*n**2 + 4.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))
        return val

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*2.0*(a + n + 1.0)*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1_diags(nr, l, a):
    """Create operator for 1st integral r^l P_n^{a,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        val = normalize_row(n,l,-1)*4.0*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))
        return val

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0)*4.0*(2.0*a - 2.0*l + 1.0)/((2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1)*8.0*(a + n + 1.0)*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qm_diags(nr, l, a):
    """Create operator for 1st integral of Q r^{l-1} P_n^{a,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(0,2)
    nzrow = 0

    # Generate main diagonal
    def d0(n):
        return -normalize_row(n,l,0,-1)*4.0*(2.0*a + 2.0*l + 2.0*n - 1.0)/(2.0*a + 2.0*l + 4.0*n - 1.0)

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,-1)*8.0*(a + n + 1.0)/(2.0*a + 2.0*l + 4.0*n + 3.0)

    ds = [d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i1qp_diags(nr, l, a):
    """Create operator for 1st integral of Q r^{l+1} P_n^{a,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1,1)*4.0*(2.0*l + 2.0*n + 1.0)/(2.0*a + 2.0*l + 4.0*n - 1.0)

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0,1)*8.0*(a + n)*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2_diags(nr, l, a):
    """Create operator for 2nd integral r^l P_n^{a,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        val = normalize_row(n,l,-2)*16.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))
        return val

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*32.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*16.0*(4.0*a**2 - 16.0*a*l - 8.0*a*n + 4.0*a + 4.0*l**2 - 8.0*l*n - 8.0*l - 8.0*n**2 - 4.0*n + 3.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1)*64.0*(a + n + 1.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*64.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2lapl_diags(nr, l, a):
    """Create operator for 2nd integral of Laplacian r^l P_n^{a,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*16.0*(2.0*l + 2.0*n - 1.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0)*32.0*(4.0*a*l + 4.0*a*n + 4.0*l*n + 2.0*l + 4.0*n**2 + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*64.0*(a + n + 1.0)**2*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qm_diags(nr, l, a):
    """Create operator for 2nd integral of Q r^{l-1} P_n^{a,l-3/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,3)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1,-1)*16.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -normalize_row(n,l,0,-1)*16.0*(2.0*a + 2.0*l + 2.0*n - 1.0)*(4.0*a - 2.0*l + 2.0*n + 3.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,-1)*64.0*(a + n + 1.0)*(a - 2.0*l - n + 1.0)/((2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2,-1)*64.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    ds = [d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i2qp_diags(nr, l, a):
    """Create operator for 2nd integral of Q r^{l+1} P_n^{a,l+1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2,1)*16.0*(2.0*l + 2.0*n - 1.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1,1)*16.0*(8.0*a*l + 8.0*a*n - 4.0*l**2 + 4.0*n**2 - 3.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0,1)*32.0*(2.0*l + 2.0*n + 1.0)*(2.0*a**2 - 4.0*a*l - 4.0*l*n - 2.0*l - 2.0*n**2 - 2.0*n + 1.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,1)*64.0*(a + n + 1.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    ds = [d_2, d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4_diags(nr, l, a):
    """Create operator for 4th integral r^l P_n^{a,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        val = normalize_row(n,l,-4)*256.0*(2.0*a + 2.0*l + 2.0*n - 7.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))
        return val

    # Generate 3rd subdiagonal
    def d_3(n):
        return normalize_row(n,l,-3)*1024.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2)*512.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(12.0*a**2 - 32.0*a*l - 8.0*a*n + 20.0*a + 12.0*l**2 - 8.0*l*n - 8.0*l - 8.0*n**2 + 12.0*n + 17.0)/((2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*1024.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(4.0*a**2 - 20.0*a*l - 12.0*a*n + 10.0*a + 4.0*l**2 - 12.0*l*n - 4.0*l - 12.0*n**2 + 6.0*n + 21.0)/((2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*256.0*(16.0*a**4 - 256.0*a**3*l - 192.0*a**3*n + 32.0*a**3 + 576.0*a**2*l**2 + 384.0*a**2*l*n - 384.0*a**2*l - 96.0*a**2*n**2 - 288.0*a**2*n + 224.0*a**2 - 256.0*a*l**3 + 384.0*a*l**2*n + 576.0*a*l**2 + 768.0*a*l*n**2 + 384.0*a*l*n - 832.0*a*l + 192.0*a*n**3 - 96.0*a*n**2 - 384.0*a*n + 208.0*a + 16.0*l**4 - 192.0*l**3*n - 128.0*l**3 - 96.0*l**2*n**2 + 192.0*l**2*n + 344.0*l**2 + 192.0*l*n**3 + 384.0*l*n**2 - 144.0*l*n - 352.0*l + 96.0*n**4 + 96.0*n**3 - 264.0*n**2 - 144.0*n + 105.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1)*2048.0*(a + n + 1.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)*(4.0*a**2 - 20.0*a*l - 12.0*a*n - 2.0*a + 4.0*l**2 - 12.0*l*n - 16.0*l - 12.0*n**2 - 18.0*n + 15.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*2048.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(12.0*a**2 - 32.0*a*l - 8.0*a*n + 4.0*a + 12.0*l**2 - 8.0*l*n - 24.0*l - 8.0*n**2 - 20.0*n + 9.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -normalize_row(n,l,3)*8192.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0))

    # Generate 4th superdiagonal
    def d4(n):
        return normalize_row(n,l,4)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(a + n + 4.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 2.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl_diags(nr, l, a):
    """Create operator for 4th integral of Laplacian r^l P_n^{a,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return normalize_row(n,l,-3)*256.0*(2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2)*1024.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(4.0*a*l + 4.0*a*n - 7.0*a - 2.0*l**2 + 5.0*l + 2.0*n**2 - 2.0*n - 6.0)/((2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*256.0*(2.0*a + 2.0*l + 2.0*n - 1.0)*(48.0*a**2*l + 48.0*a**2*n - 48.0*a**2 - 64.0*a*l**2 - 32.0*a*l*n + 136.0*a*l + 32.0*a*n**2 + 40.0*a*n - 132.0*a + 8.0*l**3 - 40.0*l**2*n - 36.0*l**2 - 56.0*l*n**2 + 64.0*l*n + 118.0*l - 8.0*n**3 + 52.0*n**2 - 14.0*n - 75.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate main diagonal
    def d0(n):
        return normalize_row(n,l,0)*512.0*(32.0*a**3*l + 32.0*a**3*n - 8.0*a**3 - 96.0*a**2*l**2 - 96.0*a**2*l*n + 120.0*a**2*l + 96.0*a**2*n - 60.0*a**2 + 32.0*a*l**3 - 96.0*a*l**2*n - 168.0*a*l**2 - 192.0*a*l*n**2 - 96.0*a*l*n + 272.0*a*l - 64.0*a*n**3 + 48.0*a*n**2 + 152.0*a*n - 82.0*a + 32.0*l**3*n + 24.0*l**3 - 96.0*l**2*n - 108.0*l**2 - 64.0*l*n**3 - 144.0*l*n**2 + 56.0*l*n + 138.0*l - 32.0*n**4 - 32.0*n**3 + 104.0*n**2 + 56.0*n - 45.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*1024.0*(a + n + 1.0)*(2.0*l + 2.0*n + 1.0)*(4.0*a**3 - 32.0*a**2*l - 20.0*a**2*n + 8.0*a**2 + 24.0*a*l**2 - 16.0*a*l*n - 76.0*a*l - 28.0*a*n**2 - 60.0*a*n + 36.0*a + 24.0*l**2*n + 36.0*l**2 + 16.0*l*n**2 - 4.0*l*n - 72.0*l - 4.0*n**3 - 32.0*n**2 - 36.0*n + 27.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -normalize_row(n,l,2)*2048.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(4.0*a**2 - 8.0*a*l + 10.0*a - 8.0*l*n - 18.0*l - 4.0*n**2 - 8.0*n + 9.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return normalize_row(n,l,3)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4lapl2_diags(nr, l, a):
    """Create operator for 4th integral bilaplacian r^l P_n^{a, l-1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2)*256.0*(2.0*l + 2.0*n - 5.0)*(2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1)*1024.0*(2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(4.0*a*l + 4.0*a*n - 4.0*a + 4.0*l*n + 2.0*l + 4.0*n**2 - 2.0*n - 5.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0)*6144.0*(a + n + 1.0)*(2.0*l + 2.0*n - 1.0)*(2.0*a*l + 2.0*a*n - a + 2.0*l*n + 2.0*l + 2.0*n**2 + n - 3.0)/((2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 1st superdiagonal
    def d1(n):
        return normalize_row(n,l,1)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(4.0*a*l + 4.0*a*n + 4.0*l*n + 6.0*l + 4.0*n**2 + 6.0*n - 3.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2)*4096.0*(a + n + 1.0)*(a + n + 2.0)**2*(a + n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qm_diags(nr, l, a):
    """Create operator for 4th integral of Q r^{l-1} P_n^{a, l-3/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,5)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return -normalize_row(n,l,-3,-1)*256.0*(2.0*a + 2.0*l + 2.0*n - 7.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -normalize_row(n,l,-2,-1)*256.0*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(8.0*a - 6.0*l + 2.0*n + 5.0)/((2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -normalize_row(n,l,-1,-1)*768.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(8.0*a**2 - 16.0*a*l + 16.0*a + 4.0*l**2 - 8.0*l*n - 8.0*l - 4.0*n**2 + 8.0*n + 11.0)/((2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate diagonal
    def d0(n):
        return  -normalize_row(n,l,0,-1)*256.0*(2.0*a + 2.0*l + 2.0*n - 1.0)*(32.0*a**3 - 144.0*a**2*l - 48.0*a**2*n + 120.0*a**2 + 96.0*a*l**2 - 96.0*a*l*n - 240.0*a*l - 96.0*a*n**2 + 208.0*a - 8.0*l**3 + 72.0*l**2*n + 60.0*l**2 + 24.0*l*n**2 - 120.0*l*n - 142.0*l - 24.0*n**3 - 60.0*n**2 + 66.0*n + 105.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate 1st superdiagonal
    def d1(n):
        return  -normalize_row(n,l,1,-1)*2048.0*(a + n + 1.0)*(2.0*a**3 - 24.0*a**2*l - 18.0*a**2*n + 6.0*a**2 + 36.0*a*l**2 + 24.0*a*l*n - 48.0*a*l - 6.0*a*n**2 - 36.0*a*n + 19.0*a - 8.0*l**3 + 12.0*l**2*n + 36.0*l**2 + 24.0*l*n**2 + 24.0*l*n - 46.0*l + 6.0*n**3 - 6.0*n**2 - 27.0*n + 15.0)/((2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return  normalize_row(n,l,2,-1)*6144.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*a**2 - 8.0*a*l - 4.0*a*n + 2.0*a + 4.0*l**2 - 8.0*l - 2.0*n**2 - 6.0*n + 3.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return  -normalize_row(n,l,3,-1)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(3.0*a - 4.0*l - n + 2.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0))

    # Generate 4th superdiagonal
    def d4(n):
        return  normalize_row(n,l,4,-1)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(a + n + 4.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i4qp_diags(nr, l, a):
    """Create operator for 4th integral of Q r^{l+1} P_n^{a, l+1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return normalize_row(n,l,-4,1)*256.0*(2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return normalize_row(n,l,-3,1)*256.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(16.0*a*l + 16.0*a*n - 28.0*a - 12.0*l**2 - 8.0*l*n + 24.0*l + 4.0*n**2 - 4.0*n - 21.0)/((2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n,l,-2,1)*768.0*(2.0*a + 2.0*l + 2.0*n - 1.0)*(16.0*a**2*l + 16.0*a**2*n - 16.0*a**2 - 32.0*a*l**2 - 32.0*a*l*n + 48.0*a*l + 16.0*a*n - 40.0*a + 8.0*l**3 - 8.0*l**2*n - 12.0*l**2 - 24.0*l*n**2 + 24.0*l*n + 46.0*l - 8.0*n**3 + 20.0*n**2 + 6.0*n - 21.0)/((2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n,l,-1,1)*256.0*(64.0*a**3*l + 64.0*a**3*n - 16.0*a**3 - 288.0*a**2*l**2 - 384.0*a**2*l*n + 192.0*a**2*l - 96.0*a**2*n**2 + 144.0*a**2*n - 120.0*a**2 + 192.0*a*l**3 - 288.0*a*l**2 - 384.0*a*l*n**2 - 192.0*a*l*n + 656.0*a*l - 192.0*a*n**3 + 48.0*a*n**2 + 416.0*a*n - 104.0*a - 16.0*l**4 + 128.0*l**3*n + 64.0*l**3 + 192.0*l**2*n**2 - 96.0*l**2*n - 344.0*l**2 - 192.0*l*n**2 - 32.0*l*n + 176.0*l - 48.0*n**4 - 48.0*n**3 + 192.0*n**2 + 72.0*n - 105.0)/((2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n,l,0,1)*512.0*(2.0*l + 2.0*n + 1.0)*(8.0*a**4 - 96.0*a**3*l - 64.0*a**3*n + 144.0*a**2*l**2 - 144.0*a**2*l - 96.0*a**2*n**2 - 144.0*a**2*n + 124.0*a**2 - 32.0*a*l**3 + 192.0*a*l**2*n + 192.0*a*l**2 + 192.0*a*l*n**2 + 96.0*a*l*n - 328.0*a*l - 96.0*a*n**2 - 80.0*a*n + 72.0*a - 32.0*l**3*n - 24.0*l**3 + 48.0*l**2*n**2 + 120.0*l**2*n + 108.0*l**2 + 96.0*l*n**3 + 168.0*l*n**2 - 112.0*l*n - 138.0*l + 24.0*n**4 + 24.0*n**3 - 96.0*n**2 - 66.0*n + 45.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n,l,1,1)*3072.0*(a + n + 1.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(4.0*a**3 - 16.0*a**2*l - 4.0*a**2*n + 4.0*a**2 + 8.0*a*l**2 - 16.0*a*l*n - 32.0*a*l - 12.0*a*n**2 - 24.0*a*n + 14.0*a + 8.0*l**2*n + 12.0*l**2 - 8.0*l*n - 24.0*l - 4.0*n**3 - 16.0*n**2 - 10.0*n + 9.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n,l,2,1)*2048.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(6.0*a**2 - 8.0*a*l + 4.0*a*n + 14.0*a - 8.0*l*n - 18.0*l - 2.0*n**2 - 4.0*n + 9.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -normalize_row(n,l,3,1)*4096.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 2.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)

def i6_diags(nr, l, a):
    """Create operator for 6th integral r^l P_n^{a,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-6,7)
    nzrow = 5

    # Generate 6th subdiagonal
    def d_6(n):
        val = normalize_row(n, l, -6)*4096.0*(2.0*a + 2.0*l + 2.0*n - 11.0)*(2.0*a + 2.0*l + 2.0*n - 9.0)*(2.0*a + 2.0*l + 2.0*n - 7.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 23.0)*(2.0*a + 2.0*l + 4.0*n - 21.0)*(2.0*a + 2.0*l + 4.0*n - 19.0)*(2.0*a + 2.0*l + 4.0*n - 17.0)*(2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0))
        return val

    # Generate 5th subdiagonal
    def d_5(n):
        return normalize_row(n, l, -5)*24576.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 9.0)*(2.0*a + 2.0*l + 2.0*n - 7.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)/((2.0*a + 2.0*l + 4.0*n - 21.0)*(2.0*a + 2.0*l + 4.0*n - 19.0)*(2.0*a + 2.0*l + 4.0*n - 17.0)*(2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return normalize_row(n, l, -4)*12288.0*(2.0*a + 2.0*l + 2.0*n - 7.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(20.0*a**2 - 48.0*a*l - 8.0*a*n + 36.0*a + 20.0*l**2 - 8.0*l*n - 8.0*l - 8.0*n**2 + 28.0*n + 31.0)/((2.0*a + 2.0*l + 4.0*n - 19.0)*(2.0*a + 2.0*l + 4.0*n - 17.0)*(2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return normalize_row(n, l, -3)*81920.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 5.0)*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(4.0*a**2 - 14.0*a*l - 6.0*a*n + 13.0*a + 4.0*l**2 - 6.0*l*n + 2.0*l - 6.0*n**2 + 15.0*n + 24.0)/((2.0*a + 2.0*l + 4.0*n - 17.0)*(2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return normalize_row(n, l, -2)*61440.0*(2.0*a + 2.0*l + 2.0*n - 3.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(16.0*a**4 - 128.0*a**3*l - 64.0*a**3*n + 96.0*a**3 + 240.0*a**2*l**2 + 96.0*a**2*l*n - 288.0*a**2*l - 48.0*a**2*n**2 + 356.0*a**2 - 128.0*a*l**3 + 96.0*a*l**2*n + 144.0*a*l**2 + 192.0*a*l*n**2 - 288.0*a*l*n - 704.0*a*l + 32.0*a*n**3 - 144.0*a*n**2 + 8.0*a*n + 396.0*a + 16.0*l**4 - 64.0*l**3*n - 48.0*l**2*n**2 + 144.0*l**2*n + 248.0*l**2 + 32.0*l*n**3 - 208.0*l*n - 192.0*l + 16.0*n**4 - 48.0*n**3 - 100.0*n**2 + 204.0*n + 225.0)/((2.0*a + 2.0*l + 4.0*n - 15.0)*(2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return normalize_row(n, l, -1)*24576.0*(2.0*a - 2.0*l + 1.0)*(2.0*a + 2.0*l + 2.0*n - 1.0)*(16.0*a**4 - 224.0*a**3*l - 160.0*a**3*n + 112.0*a**3 + 576.0*a**2*l**2 + 480.0*a**2*l*n - 576.0*a**2*l - 240.0*a**2*n + 704.0*a**2 - 224.0*a*l**3 + 480.0*a*l**2*n + 336.0*a*l**2 + 960.0*a*l*n**2 - 480.0*a*l*n - 2408.0*a*l + 320.0*a*n**3 - 480.0*a*n**2 - 1000.0*a*n + 1148.0*a + 16.0*l**4 - 160.0*l**3*n - 32.0*l**3 + 240.0*l**2*n + 584.0*l**2 + 320.0*l*n**3 - 1240.0*l*n - 568.0*l + 160.0*n**4 - 160.0*n**3 - 1120.0*n**2 + 580.0*n + 1485.0)/((2.0*a + 2.0*l + 4.0*n - 13.0)*(2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0))

    # Generate diagonal
    def d0(n):
        return normalize_row(n, l, 0)*4096.0*(64.0*a**6 - 2304.0*a**5*l - 1920.0*a**5*n + 192.0*a**5 + 14400.0*a**4*l**2 + 17280.0*a**4*l*n - 5760.0*a**4*l + 3840.0*a**4*n**2 - 4800.0*a**4*n + 4720.0*a**4 - 25600.0*a**3*l**3 - 19200.0*a**3*l**2*n + 28800.0*a**3*l**2 + 15360.0*a**3*l*n**2 + 34560.0*a**3*l*n - 55040.0*a**3*l + 10240.0*a**3*n**3 + 7680.0*a**3*n**2 - 36160.0*a**3*n + 9120.0*a**3 + 14400.0*a**2*l**4 - 19200.0*a**2*l**3*n - 38400.0*a**2*l**3 - 57600.0*a**2*l**2*n**2 - 28800.0*a**2*l**2*n + 122400.0*a**2*l**2 - 23040.0*a**2*l*n**3 + 23040.0*a**2*l*n**2 + 79680.0*a**2*l*n - 76800.0*a**2*l + 1920.0*a**2*n**4 + 15360.0*a**2*n**3 - 14400.0*a**2*n**2 - 49440.0*a**2*n + 29236.0*a**2 - 2304.0*a*l**5 + 17280.0*a*l**4*n + 14400.0*a*l**4 + 15360.0*a*l**3*n**2 - 19200.0*a*l**3*n - 68480.0*a*l**3 - 23040.0*a*l**2*n**3 - 57600.0*a*l**2*n**2 + 39360.0*a*l**2*n + 108000.0*a*l**2 - 23040.0*a*l*n**4 - 23040.0*a*l*n**3 + 119040.0*a*l*n**2 + 62400.0*a*l*n - 103312.0*a*l - 3840.0*a*n**5 + 1920.0*a*n**4 + 30080.0*a*n**3 - 18240.0*a*n**2 - 44840.0*a*n + 24708.0*a + 64.0*l**6 - 1920.0*l**5*n - 1152.0*l**5 + 3840.0*l**4*n**2 + 8640.0*l**4*n + 8080.0*l**4 + 10240.0*l**3*n**3 + 7680.0*l**3*n**2 - 36160.0*l**3*n - 27840.0*l**3 + 1920.0*l**2*n**4 - 11520.0*l**2*n**3 - 34560.0*l**2*n**2 + 24480.0*l**2*n + 48556.0*l**2 - 3840.0*l*n**5 - 11520.0*l*n**4 + 16640.0*l*n**3 + 55680.0*l*n**2 - 6200.0*l*n - 39048.0*l - 1280.0*n**6 - 1920.0*n**5 + 11680.0*n**4 + 12480.0*n**3 - 25520.0*n**2 - 14340.0*n + 10395.0)/((2.0*a + 2.0*l + 4.0*n - 11.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -normalize_row(n, l, 1)*49152.0*(a + n + 1.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)*(16.0*a**4 - 224.0*a**3*l - 160.0*a**3*n - 48.0*a**3 + 576.0*a**2*l**2 + 480.0*a**2*l*n - 96.0*a**2*l - 240.0*a**2*n + 464.0*a**2 - 224.0*a*l**3 + 480.0*a*l**2*n + 816.0*a*l**2 + 960.0*a*l*n**2 + 1440.0*a*l*n - 1928.0*a*l + 320.0*a*n**3 + 480.0*a*n**2 - 1000.0*a*n - 12.0*a + 16.0*l**4 - 160.0*l**3*n - 192.0*l**3 + 240.0*l**2*n + 824.0*l**2 + 320.0*l*n**3 + 960.0*l*n**2 - 280.0*l*n - 1488.0*l + 160.0*n**4 + 480.0*n**3 - 640.0*n**2 - 1500.0*n + 945.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 4.0*n - 9.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return normalize_row(n, l, 2)*245760.0*(a + n + 1.0)*(a + n + 2.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(16.0*a**4 - 128.0*a**3*l - 64.0*a**3*n - 32.0*a**3 + 240.0*a**2*l**2 + 96.0*a**2*l*n - 96.0*a**2*l - 48.0*a**2*n**2 - 192.0*a**2*n + 164.0*a**2 - 128.0*a*l**3 + 96.0*a*l**2*n + 336.0*a*l**2 + 192.0*a*l*n**2 + 480.0*a*l*n - 512.0*a*l + 32.0*a*n**3 + 48.0*a*n**2 - 184.0*a*n + 92.0*a + 16.0*l**4 - 64.0*l**3*n - 128.0*l**3 - 48.0*l**2*n**2 - 48.0*l**2*n + 344.0*l**2 + 32.0*l*n**3 + 192.0*l*n**2 + 176.0*l*n - 352.0*l + 16.0*n**4 + 80.0*n**3 - 4.0*n**2 - 260.0*n + 105.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n - 7.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -normalize_row(n, l, 3)*655360.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(4.0*a**2 - 14.0*a*l - 6.0*a*n - 5.0*a + 4.0*l**2 - 6.0*l*n - 16.0*l - 6.0*n**2 - 21.0*n + 15.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n - 5.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0)*(2.0*a + 2.0*l + 4.0*n + 19.0))

    # Generate 4th superdiagonal
    def d4(n):
        return normalize_row(n, l, 4)*196608.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(a + n + 4.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)*(20.0*a**2 - 48.0*a*l - 8.0*a*n + 4.0*a + 20.0*l**2 - 8.0*l*n - 40.0*l - 8.0*n**2 - 36.0*n + 15.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 2.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n - 3.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0)*(2.0*a + 2.0*l + 4.0*n + 19.0)*(2.0*a + 2.0*l + 4.0*n + 21.0))

    # Generate 5rd superdiagonal
    def d5(n):
        return -normalize_row(n, l, 5)*786432.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(a + n + 4.0)*(a + n + 5.0)*(2.0*a - 2.0*l + 1.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)*(2.0*l + 2.0*n + 9.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 2.0*n + 7.0)*(2.0*a + 2.0*l + 2.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n - 1.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0)*(2.0*a + 2.0*l + 4.0*n + 19.0)*(2.0*a + 2.0*l + 4.0*n + 21.0)*(2.0*a + 2.0*l + 4.0*n + 23.0))

    # Generate 6th superdiagonal
    def d6(n):
        return normalize_row(n, l, 6)*262144.0*(a + n + 1.0)*(a + n + 2.0)*(a + n + 3.0)*(a + n + 4.0)*(a + n + 5.0)*(a + n + 6.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)*(2.0*l + 2.0*n + 9.0)*(2.0*l + 2.0*n + 11.0)/((2.0*a + 2.0*l + 2.0*n + 1.0)*(2.0*a + 2.0*l + 2.0*n + 3.0)*(2.0*a + 2.0*l + 2.0*n + 5.0)*(2.0*a + 2.0*l + 2.0*n + 7.0)*(2.0*a + 2.0*l + 2.0*n + 9.0)*(2.0*a + 2.0*l + 2.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 3.0)*(2.0*a + 2.0*l + 4.0*n + 5.0)*(2.0*a + 2.0*l + 4.0*n + 7.0)*(2.0*a + 2.0*l + 4.0*n + 9.0)*(2.0*a + 2.0*l + 4.0*n + 11.0)*(2.0*a + 2.0*l + 4.0*n + 13.0)*(2.0*a + 2.0*l + 4.0*n + 15.0)*(2.0*a + 2.0*l + 4.0*n + 17.0)*(2.0*a + 2.0*l + 4.0*n + 19.0)*(2.0*a + 2.0*l + 4.0*n + 21.0)*(2.0*a + 2.0*l + 4.0*n + 23.0)*(2.0*a + 2.0*l + 4.0*n + 25.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)
    return (diags,offsets)
