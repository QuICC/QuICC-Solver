"""Module provides functions to generate the radial boundary conditions in a sphere with Chebyshev expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils


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

def apply_tau(mat, l, bc, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 11:
        cond = tau_diff(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 12:
        cond = tau_rdiffdivr(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 13:
        cond = tau_insulating(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 14:
        cond = tau_diff2(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 20:
        cond = tau_value_diff(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 21:
        cond = tau_value_diff2(mat.shape[0], l%2, bc.get('c',None))
    # Set last modes to zero
    elif bc[0] > 990 and bc[0] < 1000:
        cond = tau_last(mat.shape[1], l%2, bc[0]-990)
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
    cnst = c*tau_c()
    cond.append(cnst*np.ones(nr))
    if parity == 0:
        cond[-1][0] /= tau_c()

    return np.array(cond)

def tau_diff(nr, parity, coeffs = None):
    """Create the 1st derivative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    ns = np.arange(parity, 2*nr, 2)
    cond.append(c*2.0*ns**2)

    return np.array(cond)

def tau_diff2(nr, parity, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    ns = np.arange(parity, 2*nr, 2)
    cond.append(c*(2.0/3.0)*(ns**4 - ns**2))

    return np.array(cond)

def tau_rdiffdivr(nr, parity, coeffs = None):
    """Create the r D 1/r tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    ns = np.arange(parity, 2*nr, 2)
    cond.append(c*(ns**2 - 1.0)*tau_c())
    if parity == 0:
        cond[-1][0] /= tau_c()

    return np.array(cond)

def tau_insulating(nr, parity, coeffs = None):
    """Create the insulating boundray tau line(s)"""

    assert(coeffs.get('l', None) is not None)

    c = coeffs.get('c', 1.0)
    l = coeffs['l']

    cond = []
    ns = np.arange(parity, 2*nr, 2)
    cond.append(c*(2.0*ns**2 + (l+1.0)*tau_c()))
    if parity == 0:
        cond[-1][0] /= tau_c()

    return np.array(cond)

def tau_value_diff(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    cond.append(tau_value(nr,parity,coeffs)[0])
    cond.append(tau_diff(nr,parity,coeffs)[0])

    return np.array(cond)

def tau_value_diff2(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    cond.append(tau_value(nr,parity,coeffs)[0])
    cond.append(tau_diff2(nr,parity,coeffs)[0])

    return np.array(cond)

def tau_last(nr, parity, nrow):
    """Create the last modes to zero tau line(s)"""

    # convert nrow to parity aware number
    prow = nrow//2 + parity*(nrow%2)

    cond = np.zeros((prow, nr))
    for j in range(0, prow):
        cond[j,nr-prow+j] = 1.0

    return cond

def stencil(nr, l, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nr, l%2, bc.get('c',None))
    elif bc[0] == -11:
        mat = stencil_diff(nr, l%2, bc.get('c',None))
    elif bc[0] == -12:
        mat = stencil_rdiffdivr(nr, l%2, bc.get('c',None))
    elif bc[0] == -13:
        mat = stencil_insulating(nr, l%2, bc.get('c',None))
    elif bc[0] == -14:
        mat = stencil_zeroreg0(nr, l%2, bc.get('c',None))
    elif bc[0] == -15:
        mat = stencil_zeroreg1(nr, l%2, bc.get('c',None))
    elif bc[0] == -20:
        mat = stencil_value_diff(nr, l%2, bc.get('c',None))
    elif bc[0] == -21:
        mat = stencil_value_diff2(nr, l%2, bc.get('c',None))
    elif bc[0] == -22:
        mat = stencil_zeroreg2(nr, l%2, bc.get('c',None))
    elif bc[0] == -23:
        mat = stencil_zeroreg3(nr, l%2, bc.get('c',None))
    elif bc[0] == -30:
        mat = stencil_zeroreg4(nr, l%2, bc.get('c',None))
    elif bc[0] == -31:
        mat = stencil_zeroreg5(nr, l%2, bc.get('c',None))
    elif bc[0] == -40:
        mat = stencil_zeroreg6(nr, l%2, bc.get('c',None))
    elif bc[0] == -41:
        mat = stencil_zeroreg7(nr, l%2, bc.get('c',None))
    elif bc[0] == -50:
        mat = stencil_zeroreg8(nr, l%2, bc.get('c',None))
    elif bc[0] < -1 and bc[0] > -5:
        mat = restrict_eye(nr, 'cr', -bc[0])

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

def stencil_value(nr, parity, coeffs = None):
    """Create stencil matrix for a zero boundary value"""

    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        val = -np.ones(n.shape)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -1.0/2.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, parity, coeffs = None):
    """Create stencil matrix for a zero 1st derivative"""

    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -(n - 2.0)**2/n**2

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_rdiffdivr(nr, parity, coeffs = None):
    """Create stencil matrix for a zero r D 1/r"""

    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        val = -(n - 3.0)/(n + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = 1.0/6.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_insulating(nr, parity, coeffs = None):
    """Create stencil matrix for a insulating boundary"""

    assert(coeffs.get('l', None) is not None)

    l = coeffs['l']

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        val = -(l + n**2 - 4.0*n + 5.0)/(l + n**2 + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -(l + 1.0)/(l + 5.0)
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_zeroreg0(nr, parity, coeffs = None):
    """Create stencil matrix for a 0th order regularity at the origin"""

    assert(coeffs is None)
    assert(parity == 0)

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        val = np.ones(n.shape)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = 1.0/2.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

    mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
    mat = mat[:,0:nr+offsets[0]]
    return mat.tocoo()

def stencil_zeroreg1(nr, parity, coeffs = None):
    """Create stencil matrix for a 1st order regularity at the origin"""

    assert(coeffs is None)

    if parity == 0:
        return stencil_zeroreg0(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-1, 0]

        # Generate subdiagonal
        def d_1(n):
            val = (n - 2.0)/n
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg2(nr, parity, coeffs = None):
    """Create stencil matrix for a 2nd order regularity at the origin"""

    assert(coeffs is None)

    if parity == 1:
        return stencil_zeroreg1(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-2, -1, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            val = (n - 3.0)/(n - 1.0)
            for i,j in enumerate(n):
                if j == 4:
                    val[i] = 1.0/6.0
                    break
                if j > 4:
                    break
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 2.0*n/(n + 1.0)
            for i,j in enumerate(n):
                if j == 2:
                    val[i] = 2.0/3.0
                    break
                if j > 2:
                    break
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg3(nr, parity, coeffs = None):
    """Create stencil matrix for a 3rd order regularity at the origin"""

    assert(coeffs is None)

    if parity == 0:
        return stencil_zeroreg2(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-2, -1, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            val = (n - 4.0)*(n - 3.0)/((n - 1.0)*n)
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 2.0*(n - 2.0)/(n + 1.0)
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg4(nr, parity, coeffs = None):
    """Create stencil matrix for a 4th order regularity at the origin"""

    assert(coeffs is None)

    if parity == 1:
        return stencil_zeroreg3(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-3, -2, -1, 0]

        # Generate 3rd subdiagonal
        def d_3(n):
            val = (n - 5.0)*(n - 4.0)/((n - 2.0)*(n - 1.0))
            for i,j in enumerate(n):
                if j == 6:
                    val[i] = 1.0/20.0
                    break
                if j > 6:
                    break
            return val

        # Generate 2nd subdiagonal
        def d_2(n):
            val = 3.0*(n - 3.0)/(n + 1.0)
            for i,j in enumerate(n):
                if j == 4:
                    val[i] = 3.0/10.0
                    break
                if j > 4:
                    break
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 3.0*n/(n + 2.0)
            for i,j in enumerate(n):
                if j == 2:
                    val[i] = 3.0/4.0
                    break
                if j > 2:
                    break
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_3, d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg5(nr, parity, coeffs = None):
    """Create stencil matrix for a 5th order regularity at the origin"""

    assert(coeffs is None)

    if parity == 0:
        return stencil_zeroreg4(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-3, -2, -1, 0]

        # Generate 3rd subdiagonal
        def d_3(n):
            val = (n - 6.0)*(n - 5.0)*(n - 4.0)/((n - 2.0)*(n - 1.0)*n)
            return val

        # Generate 2nd subdiagonal
        def d_2(n):
            val = 3.0*(n - 4.0)*(n - 3.0)/(n*(n + 1.0))
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 3.0*(n - 2.0)/(n + 2.0)
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_3, d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg6(nr, parity, coeffs = None):
    """Create stencil matrix for a 6th order regularity at the origin"""

    assert(coeffs is None)

    if parity == 1:
        return stencil_zeroreg5(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-4, -3, -2, -1, 0]

        # Generate 4th subdiagonal
        def d_4(n):
            val = (n - 7.0)*(n - 6.0)*(n - 5.0)/((n - 3.0)*(n - 2.0)*(n - 1.0))
            for i,j in enumerate(n):
                if j == 8:
                    val[i] = 1.0/70.0
                    break
                if j > 8:
                    break
            return val

        # Generate 3rd subdiagonal
        def d_3(n):
            val = 4.0*(n - 5.0)*(n - 4.0)/((n - 1.0)*(n + 1.0))
            for i,j in enumerate(n):
                if j == 6:
                    val[i] = 4.0/35.0
                    break
                if j > 6:
                    break
            return val

        # Generate 2nd subdiagonal
        def d_2(n):
            val = 6.0*(n - 3.0)*n/((n + 1.0)*(n + 2.0))
            for i,j in enumerate(n):
                if j == 4:
                    val[i] = 2.0/5.0
                    break
                if j > 4:
                    break
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 4.0*n/(n + 3.0)
            for i,j in enumerate(n):
                if j == 2:
                    val[i] = 4.0/5.0
                    break
                if j > 2:
                    break
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_4, d_3, d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg7(nr, parity, coeffs = None):
    """Create stencil matrix for a 7th order regularity at the origin"""

    assert(coeffs is None)

    if parity == 0:
        return stencil_zeroreg6(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-4, -3, -2, -1, 0]

        # Generate 4th subdiagonal
        def d_4(n):
            val = (n - 8.0)*(n - 7.0)*(n - 6.0)*(n - 5.0)/((n - 3.0)*(n - 2.0)*(n - 1.0)*n)
            return val

        # Generate 3rd subdiagonal
        def d_3(n):
            val = 4.0*(n - 6.0)*(n - 5.0)*(n - 4.0)/((n - 1.0)*n*(n + 1.0))
            return val

        # Generate 2nd subdiagonal
        def d_2(n):
            val = 6.0*(n - 4.0)*(n - 3.0)/((n + 1.0)*(n + 2.0))
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 4.0*(n - 2.0)/(n + 3.0)
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_4, d_3, d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_zeroreg8(nr, parity, coeffs = None):
    """Create stencil matrix for a 8th order regularity at the origin"""

    assert(coeffs is None)

    if parity == 1:
        return stencil_zeroreg7(nr, parity, coeffs)

    else:
        ns = np.arange(parity,2*nr,2)
        offsets = [-5, -4, -3, -2, -1, 0]

        # Generate 5th subdiagonal
        def d_5(n):
            val = (n - 9.0)*(n - 8.0)*(n - 7.0)*(n - 6.0)/((n - 4.0)*(n - 3.0)*(n - 2.0)*(n - 1.0))
            for i,j in enumerate(n):
                if j == 10:
                    val[i] = 1.0/252.0
                    break
                if j > 10:
                    break
            return val

        # Generate 4th subdiagonal
        def d_4(n):
            val = 5.0*(n - 7.0)*(n - 6.0)*(n - 5.0)/((n - 2.0)*(n - 1.0)*(n + 1.0))
            for i,j in enumerate(n):
                if j == 8:
                    val[i] = 5.0/126.0
                    break
                if j > 8:
                    break
            return val

        # Generate 3rd subdiagonal
        def d_3(n):
            val = 10.0*(n - 5.0)*(n - 4.0)/((n + 1.0)*(n + 2.0))
            for i,j in enumerate(n):
                if j == 6:
                    val[i] = 5.0/28.0
                    break
                if j > 6:
                    break
            return val

        # Generate 2nd subdiagonal
        def d_2(n):
            val = 10.0*(n - 3.0)*n/((n + 2.0)*(n + 3.0))
            for i,j in enumerate(n):
                if j == 4:
                    val[i] = 10.0/21.0
                    break
                if j > 4:
                    break
            return val

        # Generate 1st subdiagonal
        def d_1(n):
            val = 5.0*n/(n + 4.0)
            for i,j in enumerate(n):
                if j == 2:
                    val[i] = 5.0/6.0
                    break
                if j > 2:
                    break
            return val

        # Generate diagonal
        def d0(n):
            return np.ones(n.shape)

        ds = [d_5, d_4, d_3, d_2, d_1, d0]
        diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)

        mat = spsp.diags(diags, offsets, (nr,nr)).tolil()
        mat = mat[:,0:nr+offsets[0]]
        return mat.tocoo()

def stencil_value_diff(nr, parity, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 1st derivative"""

    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-2, -1, 0]

    # Generate second subdiagonal
    def d_2(n):
        val = (n - 3.0)/(n - 1.0)
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/6.0
                break
            if j > 4:
                break
        return val

    # Generate subdiagonal
    def d_1(n):
        val = -2.0*n/(n + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -2.0/3.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff2(nr, parity, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-2, -1, 0]

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
        return val

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
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def tau_c():
    """Compute the chebyshev normalisation c factor"""

    return 2.0

def galerkin_c(n):
    """Compute the chebyshev normalisation c factor for galerkin boundary"""

    val = np.ones(n.shape)

    for i, j in enumerate(n):
        if j == 0:
            val[i] = 0.5
        if j > 0:
            break

    return val
