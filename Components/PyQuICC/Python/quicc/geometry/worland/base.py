"""Module provides generic functions required to generate values for the Worland expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.special as special


def eval_wnl(n, l, a, db, nr, func_grid, func_invnorm):

    rg = func_grid(nr)
    norm = func_invnorm(n, l)
    b = l + db

    return rg**l*special.eval_jacobi(n, a, b, 2.0*rg**2 - 1.0)*norm

def eval_pnab_origin(nr, l, db, k, normalized, func_normalize):
    """Compute the value at origin for Worland polynomials"""

    b = l + 0.5
    val = np.zeros(nr)
    val[0] = 1.0
    if nr > 0:
        for i in range(1,nr):
            val[i] = (-1.0)**i*np.exp(special.gammaln(i + b + 1.0) - special.gammaln(i + 1.0) - special.gammaln(b + 1.0))

    # Normalize
    if normalized:
        func_normalize(val, l)

    return val

def apply_norm(val, l, func_invnorm):
    """Normalize array of values"""

    for i in range(0, val.shape[0]):
         val[i] *= func_invnorm(i, l)

def compute_norm_row(n, l, k, p, func_normln):
    """Normalization factor for matrix row. l is from LHS, p allows to shift RHS to l+p"""

    if (np.ndim(n) == 0 and n + k == 0) or (np.ndim(n) > 0 and n[0] + k < 0):
        raise RuntimeError("Requested row normalization is inconsistent")

    norm = func_normln(n, l)
    norm -= func_normln(n+k, l+p)

    return np.exp(norm)

def compute_norm_row_l_1(n, l, k, func_normln):
    """Normalization factor for matrix row from W_n^l to Wn^{l-1}"""

    if (np.ndim(n) == 0 and n + k == 0) or (np.ndim(n) > 0 and n[0] + k < 0):
        raise RuntimeError("Requested row normalization is inconsistent")

    if l == 0:
        norm = func_normln(n, l+1)
    else:
        norm = func_normln(n, l-1)

    norm -= func_normln(n+k, l)

    return np.exp(norm)
