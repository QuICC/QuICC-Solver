"""Module provides functions to generate the radial boundary conditions in a cylinder with Worland expansion in radius"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import quicc.base.utils as utils
import quicc.geometry.worland.wnl as wnl


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, m, bc, pad_zeros = 0, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, m, bc, pad_zeros = pad_zeros, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, m, bc)
    else:
        bc_mat = mat

    # top row(s) restriction if required
    if bc.get('rt', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], 'rt', bc['rt'])*bc_mat

    # bottom row(s) restriction if required
    if bc.get('rb', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], 'rb', bc['rb'])*bc_mat

    # left columns restriction if required
    if bc.get('cl', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], 'cl', bc['cl'])

    # right columns restriction if required
    if bc.get('cr', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], 'cr', bc['cr'])

    # top row(s) zeroing if required
    if bc.get('zt', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[0:bc['zt'],:] = 0
        bc_mat = bc_mat.tocoo()

    # bottom row(s) zeroing if required
    if bc.get('zb', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[-bc['zb']:,:] = 0
        bc_mat = bc_mat.tocoo()

    # left columns zeroing if required
    if bc.get('zl', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[:, 0:bc['zt']] = 0
        bc_mat = bc_mat.tocoo()

    # right columns zeroing if required
    if bc.get('zr', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[:, -bc['zr']:] = 0
        bc_mat = bc_mat.tocoo()

    return bc_mat

def apply_tau(mat, m, bc, pad_zeros = 0, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], m, bc)
    elif bc[0] == 11:
        cond = tau_diff(mat.shape[0], m, bc)
    elif bc[0] == 12:
        cond = tau_rdiffdivr(mat.shape[0], m, bc)
    elif bc[0] == 14:
        cond = tau_diff2(mat.shape[0], m, bc)
    elif bc[0] == 15:
        cond = tau_laplh(mat.shape[0], m, bc)
    elif bc[0] == 16:
        cond = tau_dlaplh(mat.shape[0], m, bc)
    elif bc[0] == 17:
        cond = tau_lapl2h(mat.shape[0], m, bc)
    elif bc[0] == 18:
        cond = tau_origin(mat.shape[0], m, bc)
    elif bc[0] == 20:
        cond = tau_value_diff(mat.shape[0], m, bc)
    elif bc[0] == 21:
        cond = tau_value_diff2(mat.shape[0], m, bc)
    elif bc[0] == 22:
        cond = tau_value_laplh(mat.shape[0], m, bc)
    elif bc[0] == 30:
        cond = tau_value_diff_laplh(mat.shape[0], m, bc)
    # Set last modes to zero
    elif bc[0] > 990 and bc[0] < 1000:
        cond = tau_last(mat.shape[1], bc[0]-990)
        nbc = bc[0]-990

    if not spsp.isspmatrix_coo(mat):
        mat = mat.tocoo()
    if location == 't':
        s = 0
    elif location == 'b':
        s = mat.shape[0]-nbc

    conc = np.concatenate
    if pad_zeros > 0:
        cond = conc((np.zeros((pad_zeros,cond.shape[1])),cond))
    elif pad_zeros < 0:
        cond = conc((cond, np.zeros((pad_zeros,cond.shape[1]))))
    for i,c in enumerate(cond):
        mat.data = conc((mat.data, c))
        mat.row = conc((mat.row, [s+i]*mat.shape[1]))
        mat.col = conc((mat.col, np.arange(0,mat.shape[1])))

    return mat

def tau_value(nr, l, bc):
    """Create the boundary value tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_poly(nr, l))

    return np.array(cond)

def tau_origin(nr, l, bc):
    """Create the value at origin tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_jacobi_origin(nr, l))

    return np.array(cond)

def tau_diff(nr, l, bc):
    """Create the 1st derivative tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_diff(nr,l))

    return np.array(cond)

def tau_diff2(nr, m, bc):
    """Create the second deriviative tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_diff2(nr,m))

    return np.array(cond)

def tau_rdiffdivr(nr, m, bc):
    """Create the r D 1/r tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_rdiffdivr(nr,m))

    return np.array(cond)

def tau_lapl2h(nr, m, bc):
    """Create the no-slip tau line(s) for toroidal component on poloidal"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_lapl2h_cyl(nr,m))

    return np.array(cond)

def tau_value_diff(nr, m, bc):
    """Create the no penetration and no-slip tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_poly(nr,m))
    c = next(it)
    cond.append(c*wnl.eval_bc_diff(nr,m))

    return np.array(cond)

def tau_laplh(nr, m, bc):
    """Create the no-slip tau line(s) for toroidal component on toroidal"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_laplh_cyl(nr,m))

    return np.array(cond)

def tau_dlaplh(nr, m, bc):
    """Create the no-slip tau line(s) for toroidal component on toroidal"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_dlaplh_cyl(nr,m))

    return np.array(cond)

def tau_value_diff2(nr, m, bc):
    """Create the no penetration and no-slip tau line(s)"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_poly(nr,m))
    c = next(it)
    cond.append(c*wnl.eval_bc_diff2(nr,m))

    return np.array(cond)

def tau_value_laplh(nr, m, bc):
    """Create the no penetration and no-slip tau line(s) for poloidal component"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_poly(nr,m))
    c = next(it)
    cond.append(c*wnl.eval_bc_laplh_cyl(nr,m))

    return np.array(cond)

def tau_value_diff_laplh(nr, m, bc):
    """Create the no penetration and no-slip tau line(s) for poloidal component"""

    it = coeff_iterator(bc.get('c',None))

    cond = []
    c = next(it)
    cond.append(c*wnl.eval_bc_poly(nr,m))
    c = next(it)
    cond.append(c*wnl.eval_bc_diff(nr,m))
    c = next(it)
    cond.append(c*wnl.eval_bc_laplh_cyl(nr,m))

    return np.array(cond)

def tau_last(nr, nrow):
    """Create the last modes to zero tau line(s)"""

    cond = np.zeros((nrow, nr))
    for j in range(0, nrow):
        cond[j,nr-nrow+j] =1.0

    return cond

def stencil(nr, m, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nr, m, bc.get('c',None))
    elif bc[0] == -11:
        mat = stencil_diff(nr, m, bc.get('c',None))
    elif bc[0] == -12:
        mat = stencil_rdiffdivr(nr, m, bc.get('c',None))
    elif bc[0] == -20:
        mat = stencil_value_diff(nr, m, bc.get('c',None))
    elif bc[0] == -21:
        mat = stencil_value_diff2(nr, m, bc.get('c',None))
    elif bc[0] == -22:
        mat = stencil_value_laplh(nr, m, bc.get('c',None))
    elif bc[0] < -1 and bc[0] > -5:
        mat = restrict_eye(nr, 'cr', -bc[0])

    return mat

def apply_galerkin(mat, m, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = mat*stencil(nr, m, bc)
    return mat

def restrict_eye(nr, t, q):
    """Create the non-square identity to restrict matrix"""

    if t == 'rt':
        offsets = [q]
        diags = [[1]*(nr-q)]
        nrows = nr - q
        ncols = nr
    elif t == 'rb':
        offsets = [0]
        diags = [[1]*(nr-q)]
        nrows = nr - q
        ncols = nr
    elif t == 'cl':
        offsets = [-q]
        diags = [[1]*(nr-q)]
        nrows = nr
        ncols = nr - q
    elif t == 'cr':
        offsets = [0]
        diags = [[1]*(nr-q)]
        nrows = nr
        ncols = nr - q

    return spsp.diags(diags, offsets, (nrows, ncols))

def stencil_value(nr, m, coeffs = None):
    """Create stencil matrix for a zero boundary value"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        return -wnl.normalize_row(n,m,-1)*c(n-1.0)*2.0*n/(2.0*n - 1.0)

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, m, coeffs = None):
    """Create stencil matrix for a zero 1st derivative"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        num = -2.0*n*(4.0*(-1.0 + n)**2 + m*(-3.0 + 4.0*n))
        den = (-1.0 + 2.0*n)*(m + 4.0*m*n + 4.0*n**2)
        return wnl.normalize_row(n,m,-1)*c(n-1.0)*num/den

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_rdiffdivr(nr, l, coeffs = None):
    """Create stencil matrix for a zero r D 1/r"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        num = -2.0*n*(3.0 - 3.0*l - 8.0*n + 4.0*l*n + 4.0*n**2)
        den = (-1.0 + 2.0*n)*(-1.0 + l + 4.0*l*n + 4.0*n**2)
        return wnl.normalize_row(n,m,-1)*c(n-1.0)*num/den

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_insulating(nr, l, coeffs = None):
    """Create stencil matrix for a insulating boundary"""

    raise NotImplementedError("Boundary condition not yet implemented!")
    assert(coeffs.get('l', None) is not None)

    l = coeffs['l']

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        val = -(l + n**2 - 4.0*n + 5.0)/(l + n**2 + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -(l + 1.0)/(l + 5.0)
                break
            if j > 2:
                break
        return wnl.normalize_row(n,m,-1)*c(n-1.0)*val

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff(nr, m, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 1st derivative"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-2, -1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate second subdiagonal
    def d_2(n):
        num = 4.0*(n - 1.0)*n*(2.0*n + m - 3.0)
        den = (2.0*n - 3.0)*(2.0*n - 1.0)*(2.0*n + m - 1.0)
        return wnl.normalize_row(n,m,-2)*c(n-2.0)*num/den

    # Generate subdiagonal
    def d_1(n):
        num = 4.0*n*(2.0*n + m)
        den = (2.0*n - 1.0)*(2.0*n + m + 1.0)
        return -wnl.normalize_row(n,m,-1)*c(n-1.0)*num/den

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff2(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    raise NotImplementedError("Boundary condition not yet implemented!")
    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-2, -1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate second subdiagonal
    def d_2(n):
        val_num = (n - 3.0)*(2.0*n**2 - 12.0*n + 19.0)
        val_den = (n - 1.0)*(2.0*n**2 - 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/38.0
                break
            if j > 4:
                break
        return wnl.normalize_row(n,m,-2)*c(n-2.0)*val

    # Generate subdiagonal
    def d_1(n):
        val_num = -2.0*n*(2.0*n**2 + 7.0)
        val_den = (n + 1.0)*(2.0*n**2 + 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -10.0/19.0
                break
            if j > 2:
                break
        return wnl.normalize_row(n,m,-1)*c(n-1.0)*val

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_laplh(nr, m, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero horizontal laplacian"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-2, -1, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate second subdiagonal
    def d_2(n):
        num = 4.0*(-1.0 + n)*n*(-3.0 + m + 2.0*n)*(11.0 + 4.0*(-3.0 + n)*n + m*(-5.0 + 4.0*n))
        den = (-1.0 + m + 2.0*n)*(3.0 + 4.0*(-2.0 + n)*n)*(3.0 + 4.0*(-1.0 + n)*n + m*(-1.0 + 4.0*n))
        return wnl.normalize_row(n,m,-2)*c(n-2.0)*num/den

    # Generate subdiagonal
    def d_1(n):
        num = -4.0*n*(m + 2.0*n)*(5.0 + m + 4.0*m*n + 4.0*n**2)
        den = (-1.0 + 2.0*n)*(1.0 + m + 2.0*n)*(3.0 + 4.0*n*(1.0 + n) + m*(3.0 + 4.0*n))
        return wnl.normalize_row(n,m,-1)*c(n-1.0)*num/den

    # Generate diagonal
    def d0(n):
        return wnl.normalize_row(n,m,0)*c(n)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def coeff_iterator(coeffs):
    """Return an iterator over the constants"""

    if coeffs is None:
        it = itertools.cycle([1.0])
    else:
        try:
            if len(coeffs) > 1:
                it = iter(coeffs)
            elif len(coeffs) == 1:
                it = itertools.cycle(coeffs)
            else:
                raise RuntimeError
        except:
            it = itertools.cycle([coeffs])

    return it
