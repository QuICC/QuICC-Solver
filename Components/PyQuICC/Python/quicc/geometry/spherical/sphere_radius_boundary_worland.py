"""Module provides functions to generate the radial boundary conditions in a sphere with Worland expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.worland.wnl as wnl

def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, l, bc, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, l, bc, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, l, bc)
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

def tau_value(nr, l, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_poly(nr, l))

    return np.array(cond)

def tau_diff(nr, l, coeffs = None):
    """Create the 1st derivative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_diff(nr,l))

    return np.array(cond)

def tau_diff2(nr, l, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_diff2(nr,l))

    return np.array(cond)

def tau_diff3(nr, l, coeffs = None):
    """Create the 3. deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_diff3(nr,l))

    return np.array(cond)

def tau_diff4(nr, l, coeffs = None):
    """Create the 4. deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_diff4(nr,l))

    return np.array(cond)

def tau_rdiffdivr(nr, l, coeffs = None):
    """Create the r D 1/r tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_rdiffdivr(nr,l))

    return np.array(cond)

def tau_insulating(nr, l, coeffs = None):
    """Create the insulating boundray tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_insulating_sph(nr,l))

    return np.array(cond)

def tau_value_diff(nr, l, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_poly(nr,l))
    cond.append(c*wnl.eval_bc_diff(nr,l))

    return np.array(cond)

def tau_value_diff2(nr, l, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wnl.eval_bc_poly(nr,l))
    cond.append(c*wnl.eval_bc_diff2(nr,l))

    return np.array(cond)

def tau_last(nr, nrow):
    """Create the last modes to zero tau line(s)"""

    cond = np.zeros((nrow, nr))
    for j in range(0, nrow):
        cond[j,nr-nrow+j] =1.0

    return cond

# Tau line map
tau_map = dict()
tau_map[10] = (lambda nr, l, bc: tau_value(nr, l, bc.get('c',None)))
tau_map[11] = (lambda nr, l, bc: tau_diff(nr, l, bc.get('c',None)))
tau_map[12] = (lambda nr, l, bc: tau_rdiffdivr(nr, l, bc.get('c',None)))
tau_map[13] = (lambda nr, l, bc: tau_insulating(nr, l, bc.get('c',None)))
tau_map[14] = (lambda nr, l, bc: tau_diff2(nr, l, bc.get('c',None)))
tau_map[15] = (lambda nr, l, bc: tau_diff3(nr, l, bc.get('c',None)))
tau_map[16] = (lambda nr, l, bc: tau_diff4(nr, l, bc.get('c',None)))
tau_map[20] = (lambda nr, l, bc: tau_value_diff(nr, l, bc.get('c',None)))
tau_map[21] = (lambda nr, l, bc: tau_value_diff2(nr, l, bc.get('c',None)))

def apply_tau(mat, l, bc, location = 't'):
    """Add Tau lines to the matrix"""

    global tau_map
    nbc = bc[0]//10

    # Set last modes to zero
    if bc[0] > 990 and bc[0] < 1000:
        cond = tau_last(mat.shape[1], bc[0]-990)
        nbc = bc[0]-990
    # Apply condition from map
    elif bc[0] > 0:
        fct = tau_map.get(bc[0], None)
        if fct is None:
            raise RuntimeError("Unknown tau line " + str(bc[0]))
        cond = fct(mat.shape[0], l, bc)

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

def apply_galerkin(mat, l, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = mat*stencil(nr, l, bc)
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

def stencil_value(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_value_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, l, coeffs = None):
    """Create stencil matrix for a zero 1st derivative"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_diff_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_rdiffdivr(nr, l, coeffs = None):
    """Create stencil matrix for a zero r D 1/r"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_rdiffdivr_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_insulating(nr, l, coeffs = None):
    """Create stencil matrix for a insulating boundary"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_insulating_sph_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 1st derivative"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_value_diff_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff2(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    assert(coeffs is None)
    diags, offsets = wnl.stencil_value_diff2_diags(nr, l)

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

# Galerkin stencil map
stencil_map = dict()
stencil_map[-10] = stencil_value
stencil_map[-11] = stencil_diff
stencil_map[-12] = stencil_rdiffdivr
stencil_map[-13] = stencil_insulating
stencil_map[-20] = stencil_value_diff
stencil_map[-21] = stencil_value_diff2

def stencil(nr, l, bc):
    """Create a Galerkin stencil matrix"""

    global stencil_map

    if bc[0] < -1 and bc[0] > -5:
        mat = restrict_eye(nr, 'cr', -bc[0])
    elif bc[0] < 0:
        mat = stencil_map[bc[0]](nr, l, bc.get('c',None))

    return mat
