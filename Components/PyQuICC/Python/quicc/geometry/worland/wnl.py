"""Module provides functions to generate values for the Worland expansion type agnostic"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.special as special

import quicc.base.utils as utils
try:
    import quicc.geometry.worland.setup as wsetup
    worland_type = wsetup.type
except ImportError:
    worland_type = "Chebyshev"

if worland_type == "Chebyshev":
    from quicc.geometry.worland.chebyshev import (get_grid, get_weights, get_norm, get_invnorm,
            normalize, normalize_row, normalize_row_l_1,
            jacobi_alpha, jacobi_beta,
            eval_poly, eval_jacobi_origin, eval_bc_poly,
            r2_diags, i1_diags, i1qm_diags, i1qp_diags, i2_diags, i2lapl_diags, i2qm_diags, i2qp_diags,
            i4_diags, i4lapl_diags, i4lapl2_diags, i4qm_diags, i4qp_diags, i6_diags,
            stencil_value_diags, stencil_diff_diags, stencil_rdiffdivr_diags,
            stencil_insulating_sph_diags, stencil_value_diff_diags, stencil_value_diff2_diags)
elif worland_type == "Legendre":
    from quicc.geometry.worland.legendre import (get_grid, get_weights, get_norm, get_invnorm,
            normalize, normalize_row, normalize_row_l_1,
            jacobi_alpha, jacobi_beta,
            eval_poly, eval_jacobi_origin, eval_bc_poly,
            r2_diags, i1_diags, i1qm_diags, i1qp_diags, i2_diags, i2lapl_diags, i2qm_diags, i2qp_diags,
            i4_diags, i4lapl_diags, i4lapl2_diags, i4qm_diags, i4qp_diags, i6_diags,
            stencil_value_diags, stencil_diff_diags, stencil_rdiffdivr_diags,
            stencil_insulating_sph_diags, stencil_value_diff_diags, stencil_value_diff2_diags)
elif worland_type == "CylEnergy":
    from quicc.geometry.worland.cylenergy import (get_grid, get_weights, get_norm, get_invnorm,
            normalize, normalize_row, normalize_row_l_1,
            jacobi_alpha, jacobi_beta,
            eval_poly, eval_jacobi_origin, eval_bc_poly,
            r2_diags, i1_diags, i1qm_diags, i1qp_diags, i2_diags, i2lapl_diags, i2qm_diags, i2qp_diags,
            i4_diags, i4lapl_diags, i4lapl2_diags, i4qm_diags, i4qp_diags, i6_diags,
            stencil_value_diags, stencil_diff_diags, stencil_rdiffdivr_diags,
            stencil_insulating_sph_diags, stencil_value_diff_diags, stencil_value_diff2_diags)
elif worland_type == "SphEnergy":
    from quicc.geometry.worland.sphenergy import (get_grid, get_weights, get_norm, get_invnorm,
            normalize, normalize_row, normalize_row_l_1,
            jacobi_alpha, jacobi_beta,
            eval_poly, eval_jacobi_origin, eval_bc_poly,
            r2_diags, i1_diags, i1qm_diags, i1qp_diags, i2_diags, i2lapl_diags, i2qm_diags, i2qp_diags,
            i4_diags, i4lapl_diags, i4lapl2_diags, i4qm_diags, i4qp_diags, i6_diags,
            stencil_value_diags, stencil_diff_diags, stencil_rdiffdivr_diags,
            stencil_insulating_sph_diags, stencil_value_diff_diags, stencil_value_diff2_diags)
elif worland_type == "Generic":
    from quicc.geometry.worland.generic import (get_grid, get_weights, get_norm, get_invnorm,
            normalize, normalize_row, normalize_row_l_1,
            jacobi_alpha, jacobi_beta,
            eval_poly, eval_jacobi_origin, eval_bc_poly,
            r2_diags, i1_diags, i1qm_diags, i1qp_diags, i2_diags, i2lapl_diags, i2qm_diags, i2qp_diags,
            i4_diags, i4lapl_diags, i4lapl2_diags, i4qm_diags, i4qp_diags, i6_diags,
            stencil_value_diags, stencil_diff_diags, stencil_rdiffdivr_diags,
            stencil_insulating_sph_diags, stencil_value_diff_diags, stencil_value_diff2_diags)


def eval_bc_diff(nr, l):
    """Compute the first derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ns = np.arange(1,nr)
        if nr > 1:
            val[1:] = 2.0*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        if l > 0:
            val += l*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_diff2(nr, l):
    """Compute the second derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ab2 = a + b + 2.0
        ns = np.arange(2,nr)
        if nr > 2:
            val[2:] = 4.0*(ab1 + ns)*(ab2 + ns)*eval_bc_poly(nr-2, l+2, 2, False)

        ns = np.arange(1,nr)
        if l > 0:
            if nr > 1:
                val[1:] += (2.0*(ab1 + ns)*(1.0 + 2.0*l))*eval_bc_poly(nr-1, l+1, 1, False)
        else:
            if nr > 1:
                val[1:] += 2.0*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        if l > 1:
            val += l*(l-1.0)*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_diff3(nr, l):
    """Compute the third derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ab2 = a + b + 2.0
        ab3 = a + b + 3.0
        ns = np.arange(3,nr)
        if nr > 3:
            val[3:] = 8.0*(ab1 + ns)*(ab2 + ns)*(ab3 + ns)*eval_bc_poly(nr-3, l+3, 3, False)

        ns = np.arange(2,nr)
        if nr > 2:
            val[2:] += 12.0*(ab1 + ns)*(ab2 + ns)*(1.0 + l)*eval_bc_poly(nr-2, l+2, 2, False)

        if l > 0:
            if nr > 1:
                ns = np.arange(1,nr)
                val[1:] += 6.0*l**2*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        if l > 2:
            val += (l-1.0)*(l-2.0)*l*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_diff4(nr, l):
    """Compute the fourth derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ab2 = a + b + 2.0
        ab3 = a + b + 3.0
        ab4 = a + b + 4.0
        ns = np.arange(4,nr)
        if nr > 4:
            val[4:] = 16.0*(ab4 + ns)*(ab3 + ns)*(ab2 + ns)*(ab1 + ns)*eval_bc_poly(nr-4, l+4, 4, False)

        ns = np.arange(3,nr)
        if nr > 3:
            val[3:] += 16.0*(2.0*l + 3.0)*(ab3 + ns)*(ab2 + ns)*(ab1 + ns)*eval_bc_poly(nr-3, l+3, 3, False)

        ns = np.arange(2,nr)
        if nr > 2:
            val[2:] += 12.0*(2.0*l**2 + 2.0*l + 1.0)*(ab2 + ns)*(ab1 + ns)*eval_bc_poly(nr-2, l+2, 2, False)
        if l > 1:
            if nr > 1:
                ns = np.arange(1,nr)
                val[1:] += 4.0*(l-1.0)*l*(2.0*l-1.0)*(ab1 + ns)*eval_bc_poly(nr-1, l, 1, False)
        if l > 3:
            val += l*(l-1.0)*(l-2.0)*(l-3.0)*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_rdiffdivr(nr, l):
    """Compute the stress-free condition at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ns = np.arange(1,nr)
        if nr > 1:
            val[1:] = 2.0*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        val += (l-1.0)*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_divrdiffr(nr, l):
    """Compute the 1/r D r condition at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        if nr > 1:
            val[1:] = 2.0*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        val += (l+1.0)*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_insulating_sph(nr, l):
    """Compute the insulating magnetic condition for a sphere at endpoint for Worland polynomials (D + l+1/r)"""

    val = np.zeros(nr)
    if nr > 0:
        a = jacobi_alpha(l)
        b = jacobi_beta(l)
        ab1 = a + b + 1.0
        ns = np.arange(1,nr)
        if nr > 1:
            val[1:] = 2.0*(ab1 + ns)*eval_bc_poly(nr-1, l+1, 1, False)
        val += (2.0*l+1.0)*eval_bc_poly(nr, l, 0, False)

    # Normalize
    normalize(val, l)

    return val

def eval_bc_laplh_cyl(nr, m):
    """Compute the horizontal laplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        if nr > 2:
            val[2:] = 4.0*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*eval_bc_poly(nr-2, m+2, 2, False)
        if nr > 1:
            val[1:] += 4.0*(m+1.0)*(m+np.arange(1,nr))*eval_bc_poly(nr-1, m+1, 1, False)

    # Normalize
    normalize(val, m)

    return val

def eval_bc_dlaplh_cyl(nr, m):
    """Compute the radial derivative of horizontal laplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        if nr > 3:
            val[3:] = 8.0*(m+np.arange(3,nr)+2)*(m+np.arange(3,nr)+1)*(m+np.arange(3,nr))*eval_bc_poly(nr-3, m+3, 3, False)
        if nr > 2:
            val[2:] += 4.0*(3.0*m + 4.0)*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*eval_bc_poly(nr-2, m+2, 2, False)
        if m > 0:
            if nr > 1:
                val[1:] += 4.0*m*(m + 1.0)*(m+np.arange(1,nr))*eval_bc_poly(nr-1, m+1, 1, False)

    # Normalize
    normalize(val, m)

    return val

def eval_bc_lapl2h_cyl(nr, m):
    """Compute the bilaplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        if nr > 4:
            val[4:] = 16.0*(m+np.arange(4,nr)+3)*(m+np.arange(4,nr)+2)*(m+np.arange(4,nr)+1)*(m+np.arange(4,nr))*eval_bc_poly(nr-4, m+4, 4, False)
        if nr > 3:
            val[3:] += 32.0*(m + 2.0)*(m+np.arange(3,nr)+2)*(m+np.arange(3,nr)+1)*(m+np.arange(3,nr))*eval_bc_poly(nr-3, m+3, 3, False)
        if nr > 2:
            val[2:] += 16.0*(m + 1.0)*(m + 2.0)*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*eval_bc_poly(nr-2, m+2, 2, False)

    # Normalize
    normalize(val, m)

    return val
