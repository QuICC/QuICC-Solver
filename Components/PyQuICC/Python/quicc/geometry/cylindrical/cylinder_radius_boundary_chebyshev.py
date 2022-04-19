"""Module provides functions to generate the radial boundary conditions in a cylinder with Chebyshev expansion in radius"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, parity, bc, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, parity, bc, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, parity, bc)
    else:
        bc_mat = mat

    # top row(s) restriction if required
    if bc.get('rt', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], parity, 'rt', bc['rt'])*bc_mat

    # bottom row(s) restriction if required
    if bc.get('rb', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], parity, 'rb', bc['rb'])*bc_mat

    # left columns restriction if required
    if bc.get('cl', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], parity, 'cl', bc['cl'])

    # right columns restriction if required
    if bc.get('cr', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], parity, 'cr', bc['cr'])

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

def apply_tau(mat, parity, bc, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    # u = 0
    if bc[0] == 10:
        cond = tau_value(mat.shape[0], parity, bc.get('c',None))
    # D u = 0
    elif bc[0] == 11:
        cond = tau_diff(mat.shape[0], parity, bc.get('c',None))
    # 1/r D u = 0
    elif bc[0] == 13:
        cond = tau_1rdr(mat.shape[0], parity, bc.get('c',None))
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
    for i,c in enumerate(cond):
        mat.data = conc((mat.data, c))
        mat.row = conc((mat.row, [s+i]*mat.shape[1]))
        mat.col = conc((mat.col, np.arange(0,mat.shape[1])))

    return mat

def tau_value(nr, parity, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*tau_c(i) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_diff(nr, parity, coeffs = None):
    """Create the first derivative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*i**2 for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_diff2(nr, parity, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*((1.0/3.0)*(i**4 - i**2)) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_1rdr(nr, parity, coeffs = None):
    """Create the 1/r D_r tau line(s)"""

    cond = tau_value(nr, parity, coeffs) + tau_diff(nr, parity, coeffs)

    return cond

def tau_last(nr, parity):
    """Create the zero last mode tau line(s)"""

    cond = []
    cond.append([0 for i in np.arange(parity, 2*(nr-1), 2)] + [tau_c(2*nr)])

    return np.array(cond)

def stencil(nr, parity, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nr, parity)
    elif bc[0] == -11:
        mat = stencil_diff(nr, parity)
    elif bc[0] == -13:
        mat = stencil_1rdr(nr, parity)

    return mat

def apply_galerkin(mat, parity, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = restrict_eye(nr, parity, 'rt', bc['r'])*mat*stencil(nr, parity, bc)
    return mat

def restrict_eye(nr, parity, t, q):
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

def stencil_value(nr, parity):
    """Create stencil matrix for a zero boundary value"""

    ns = np.arange(parity, 2*nr, 2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -galerkin_c(n+2*offsets[0])

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, parity):
    """Create stencil matrix for a zero 1st derivative"""

    ns = np.arange(parity, 2*nr, 2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -(n+2.0*offsets[0])**2/n**2

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff2(nr, parity):
    """Create stencil matrix for a zero 2nd derivative"""

    ns = np.arange(parity, 2*nr, 2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -(n - 3.0)*(n - 2.0)**2/(n**2*(n + 1.0))

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_1rdr(nr, parity):
    """Create stencil matrix for a zero 2nd derivative"""

    ns = np.arange(parity, 2*nr, 2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -(n - 3.0)/(n + 1.0)

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def tau_c(n):
    """Compute the chebyshev normalisation c factor"""

    if n > 0:
        return 2
    else:
        return 1

def galerkin_c(n):
    """Compute the chebyshev normalisation c factor for galerkin boundary"""

    if n > 0:
        return 1
    else:
        return 0.5
