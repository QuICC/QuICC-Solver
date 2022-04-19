"""Module provides functions to generate the boundary conditions in a cartesian 1D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import quicc.base.utils as utils


default_use_parity = False

def linear_map(xi, xo):
    """Calculat a and b for linear map x = a*y + b"""

    a = (xo - xi)/2.0;
    b = (xo + xi)/2.0;

    return (a, b)

def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, xi, xo, bc, pad_zeros = 0, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, xi, xo, bc, pad_zeros = pad_zeros, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, xi, xo, bc)
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

def apply_tau(mat, xi, xo, bc, pad_zeros = 0, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[1], xi, xo, 1, bc)
    elif bc[0] == 11:
        cond = tau_value(mat.shape[1], xi, xo, -1, bc)
    elif bc[0] == 12:
        cond = tau_diff(mat.shape[1], xi, xo, 1, bc)
    elif bc[0] == 13:
        cond = tau_diff(mat.shape[1], xi, xo, -1, bc)
    elif bc[0] == 14:
        cond = tau_diff2(mat.shape[1], xi, xo, 1, bc)
    elif bc[0] == 15:
        cond = tau_diff2(mat.shape[1], xi, xo, -1, bc)
    elif bc[0] == 16:
        cond = tau_integral(mat.shape[1], xi, xo, 1, bc)
    elif bc[0] == 20:
        cond = tau_value(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 21:
        cond = tau_diff(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 22:
        cond = tau_valuediff(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 23:
        cond = tau_diff2(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 24:
        cond = tau_insulating(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 25:
        cond = tau_lapl(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 40:
        cond = tau_value_diff(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 41:
        cond = tau_value_diff2(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 42:
        cond = tau_diff_diff2(mat.shape[1], xi, xo, 0, bc)
    elif bc[0] == 43:
        cond = tau_value_lapl(mat.shape[1], xi, xo, 0, bc)
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

def tau_value(nx, xi, xo, pos, bc):
    """Create the tau line(s) for a zero boundary value"""

    it = coeff_iterator(bc.get('c',None), pos)

    cond = []
    c = next(it)
    if pos >= 0:
        cnst = c*tau_c()
        cond.append(cnst*np.ones(nx))
        cond[-1][0] /= tau_c()
        c = next(it)

    if pos <= 0:
        cnst = c*tau_c()
        cond.append(cnst*alt_ones(nx, 1))
        cond[-1][0] /= tau_c()

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_diff(nx, xi, xo, pos, bc):
    """Create the tau line(s) for a zero 1st derivative"""

    it = coeff_iterator(bc.get('c',None), pos)
    a,b = linear_map(xi, xo)

    cond = []
    c = next(it)
    ns = np.arange(0,nx)
    if pos >= 0:
        cnst = c*(tau_c()/a)
        cond.append(cnst*ns**2)
        c = next(it)

    if pos <= 0:
        cnst = c*(tau_c()/a)
        cond.append(cnst*ns**2*alt_ones(nx, 0))

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_diff2(nx, xi, xo, pos, bc):
    """Create the tau line(s) for a zero 2nd derivative"""

    it = coeff_iterator(bc.get('c',None), pos)
    a,b = linear_map(xi, xo)

    cond = []
    c = next(it)
    ns = np.arange(0,nx)
    if pos >= 0:
        cnst = c*(tau_c()/a**2)
        cond.append((cnst/3.0)*(ns**4 - ns**2))
        c = next(it)

    if pos <= 0:
        cnst = c*(tau_c()/a**2)
        cond.append((cnst/3.0)*(ns**4 - ns**2)*alt_ones(nx, 1))

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_insulating(nx, xi, xo, pos, bc):
    """Create the tau line(s) for a insulating boundary"""

    assert(bc.get('k_perp', None) is not None)

    k_perp = bc['k_perp']
    it = coeff_iterator(bc.get('c',None), pos)
    a,b = linear_map(xi, xo)

    cond = []
    c = next(it)
    ns = np.arange(0,nx)
    if pos >= 0:
        # D + a P
        c_P = k_perp*c
        cond.append(c_P*np.ones(nx))
        cond[-1][0] /= tau_c()
        c_D = c/a
        cond[-1] += c_D*ns**2
        c = next(it)*tau_c()

    if pos <= 0:
        # D - a P
        c_P = -k_perp*c
        cond.append(c_P*alt_ones(nx, 1))
        cond[-1][0] /= tau_c()
        c_D = c/a
        cond[-1] += c_D*ns**2*alt_ones(nx, 0)

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_lapl(nx, xi, xo, pos, bc):
    """Create the tau line(s) for a zero laplacian"""

    assert(bc.get('kx', None) is not None)
    assert(bc.get('ky', None) is not None)

    kx = bc.get('kx')
    ky = bc.get('ky')

    cond = tau_diff2(nx, xi, xo, pos, bc) - (kx**2 + ky**2)*tau_value(nx, xi, xo, pos, bc)

    return np.array(cond)

def tau_integral(nx, xi, xo, pos, bc):
    """Create the boundary integral tau line(s)"""

    it = coeff_iterator(bc.get('c',None), pos)

    cond = []
    c = next(it)
    ns = np.arange(0,nx,2)
    if pos >= 0:
        tmp = np.zeros(nx)
        tmp[::2] = 2.0*(ns/(ns**2 - 1.0) - 1.0/(ns - 1.0))
        tmp[0] = tmp[0]/2.0
        cond.append(tmp)
        c = next(it)

    if pos <= 0:
        tmp = np.zeros(nx)
        tmp[::2] = 2.0*(ns/(ns**2 - 1.0) - 1.0/(ns - 1.0))
        tmp[0] = tmp[0]/2.0
        cond.append(tmp)

    return np.array(cond)

def tau_valuediff(nx, xi, xo, pos, bc):
    """Create the tau lines for a zero boundary value at top and a zero 1st derivative at bottom"""

    assert(pos == 0)
    itv = coeff_iterator(bc.get('cv', None), pos)
    itd = coeff_iterator(bc.get('cd', None), pos)

    cond = []
    bcv = dict()
    bcd = dict()
    bcv['c'] = next(itv)
    bcd['c'] = next(itd)
    cond.append(tau_value(nx, xi, xo, 1, bcv)[0])
    cond.append(tau_diff(nx, xi, xo, -1, bcd)[0])

    return np.array(cond)

def tau_value_diff(nx, xi, xo, pos, bc):
    """Create the tau lines for a zero boundary value and a zero 1st derivative"""

    itv = coeff_iterator(bc.get('cv', None), pos)
    itd = coeff_iterator(bc.get('cd', None), pos)

    cond = []
    bcv = dict()
    bcd = dict()
    bcv['c'] = next(itv)
    bcd['c'] = next(itd)
    if pos >= 0:
        cond.append(tau_value(nx, xi, xo, 1, bcv)[0])
        cond.append(tau_diff(nx, xi, xo, 1, bcd)[0])
        bcv['c'] = next(itv)
        bcd['c'] = next(itd)

    if pos <= 0:
        cond.append(tau_value(nx, xi, xo, -1, bcv)[0])
        cond.append(tau_diff(nx, xi, xo, -1, bcd)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    if pos == 0:
        order = [0, 2, 1, 3]
        cond = [cond[i] for i in order]

    return np.array(cond)

def tau_diff_diff2(nx, xi, xo, pos, bc):
    """Create tau lines for a zero 1st derivative and a zero 2nd derivative """

    itd = coeff_iterator(bc.get('cd', None), pos)
    itd2 = coeff_iterator(bc.get('cd2', None), pos)

    cond = []
    bcd = dict()
    bcd2 = dict()
    bcd['c'] = next(itd)
    bcd2['c'] = next(itd2)
    if pos >= 0:
        cond.append(tau_diff(nx, xi, xo, 1, bcd)[0])
        cond.append(tau_diff2(nx, xi, xo, 1, bcd2)[0])
        bcd['c'] = next(itd)
        bcd2['c'] = next(itd2)

    if pos <= 0:
        cond.append(tau_diff(nx, xi, xo, -1, bcd)[0])
        cond.append(tau_diff2(nx, xi, xo, -1, bcd2)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    if pos == 0:
        order = [0, 2, 1, 3]
        cond = [cond[i] for i in order]

    return np.array(cond)

def tau_value_diff2(nx, xi, xo, pos, bc):
    """Create tau lines for a zero boundary value and a zero 2nd derivative """

    itv = coeff_iterator(bc.get('cv', None), pos)
    itd2 = coeff_iterator(bc.get('cd2', None), pos)

    cond = []
    bcv = dict()
    bcd2 = dict()
    bcv['c'] = next(itv)
    bcd2['c'] = next(itd2)
    if pos >= 0:
        cond.append(tau_value(nx, xi, xo, 1, bcv)[0])
        cond.append(tau_diff2(nx, xi, xo, 1, bcd2)[0])
        bcv['c'] = next(itv)
        bcd2['c'] = next(itd2)

    if pos <= 0:
        cond.append(tau_value(nx, xi, xo, -1, bcv)[0])
        cond.append(tau_diff2(nx, xi, xo, -1, bcd2)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    if pos == 0:
        order = [0, 2, 1, 3]
        cond = [cond[i] for i in order]

    return np.array(cond)

def tau_value_lapl(nx, xi, xo, pos, bc):
    """Create the tau lines for a zero boundary value and a zero 1st derivative"""

    itv = coeff_iterator(bc.get('cv', None), pos)
    itl = coeff_iterator(bc.get('cl', None), pos)

    cond = []
    bcv = dict()
    bcl = dict()
    bcv['c'] = next(itv)
    bcl['c'] = next(itl)
    if pos >= 0:
        cond.append(tau_value(nx, xi, xo, 1, bcv)[0])
        cond.append(tau_lapl(nx, xi, xo, 1, bcl)[0])
        bcv['c'] = next(itv)
        bcl['c'] = next(itl)

    if pos <= 0:
        cond.append(tau_value(nx, xi, xo, -1, bcv)[0])
        cond.append(tau_lapl(nx, xi, xo, -1, bcl)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    if pos == 0:
        order = [0, 2, 1, 3]
        cond = [cond[i] for i in order]

    return np.array(cond)

def tau_last(nx, nrow):
    """Create the last modes to zero value tau line(s)"""

    cond = np.zeros((nrow, nx))
    for j in range(0, nrow):
        cond[j,nx-nrow+j] = tau_c()

    return cond

def stencil(nx, xi, xo, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] < 0 and bc[0] > -10:
        mat = restrict_eye(nx, 'cr', abs(bc[0]))
    elif bc[0] == -10:
        mat = stencil_value(nx, xi, xo, 1)
    elif bc[0] == -11:
        mat = stencil_value(nx, xi, xo, -1)
    elif bc[0] == -12:
        mat = stencil_value(nx, xi, xo, 1)
    elif bc[0] == -13:
        mat = stencil_value(nx, xi, xo, -1)
    elif bc[0] == -20:
        mat = stencil_value(nx, xi, xo, 0)
    elif bc[0] == -21:
        mat = stencil_diff(nx, xi, xo, 0)
    elif bc[0] == -23:
        mat = stencil_diff2(nx, xi, xo, 0)
    elif bc[0] == -40:
        mat = stencil_value_diff(nx, xi, xo, 0)
    elif bc[0] == -41:
        mat = stencil_value_diff2(nx, xi, xo, 0)

    return mat

def apply_galerkin(mat, xi, xo, bc):
    """Apply a Galerkin stencil on the matrix"""

    nx = mat.shape[0]
    mat = mat*stencil(nx, xi, xo, bc)
    return mat

def restrict_eye(nx, t, q):
    """Create the non-square identity to restrict matrix"""

    if t == 'rt':
        offsets = [q]
        diags = [[1]*(nx-q)]
        nrows = nx - q
        ncols = nx
    elif t == 'rb':
        offsets = [0]
        diags = [[1]*(nx-q)]
        nrows = nx - q
        ncols = nx
    elif t == 'cl':
        offsets = [-q]
        diags = [[1]*(nx-q)]
        nrows = nx
        ncols = nx - q
    elif t == 'cr':
        offsets = [0]
        diags = [[1]*(nx-q)]
        nrows = nx
        ncols = nx - q

    return spsp.diags(diags, offsets, (nrows, ncols))

def stencil_value(nx, xi, xo, pos):
    """Create stencil matrix for a zero boundary value"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        return galerkin_c(n+offsets[0])*c(n+offsets[0])*sgn

    # Generate diagonal
    def d0(n):
        return c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))

def stencil_diff(nx, xi, xo, pos):
    """Create stencil matrix for a zero 1st derivative"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos

    def c(n):
        return np.ones(n.shape)

    # Generate subdiagonal
    def d_1(n):
        return sgn*c(n+offsets[0])*(n+offsets[0])**2/n**2

    # Generate diagonal
    def d0(n):
        return c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))

def stencil_diff2(nx, xi, xo, pos):
    """Create stencil matrix for a zero 2nd derivative"""

    def c(n):
        return np.ones(n.shape)

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]

        # Generate subdiagonal
        def d_1(n):
            return -c(n-2.0)*(n - 3.0)*(n - 2.0)**2/(n**2*(n + 1.0))

    else:
        offsets = [-1, 0]

        # Generate subdiagonal
        def d_1(n):
            return -c(n-1.0)*pos*(n - 2.0)*(n - 1.0)/(n*(n + 1.0))

    # Generate diagonal
    def d0(n):
        return c(n)*np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))

def stencil_value_diff(nx, xi, xo, pos):
    """Create stencil matrix for a zero boundary value and a zero 1st derivative"""

    assert(pos == 0)

    ns = np.arange(0,nx,1)
    offsets = [-4, -2, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate 2nd subdiagonal
    def d_2(n):
        val = (n - 3.0)/(n - 1.0)
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/6.0
            if j > 4:
                break

        return c(n-4.0)*val

    # Generate 1st subdiagonal
    def d_1(n):
        val = -2.0*n/(n + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -2.0/3.0
            if j > 2:
                break

        return c(n-2.0)*val

    # Generate diagonal
    def d0(n):
        return c(n)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))

def stencil_value_diff2(nx, xi, xo, pos):
    """Create stencil matrix for a zero boundary value and a zero 2nd derivative"""

    assert(pos == 0)

    ns = np.arange(0,nx,1)
    offsets = [-4, -2, 0]

    def c(n):
        return np.ones(n.shape)

    # Generate 2nd subdiagonal
    def d_2(n):
        val_num = (n - 3.0)*(2.0*n**2 - 12.0*n + 19.0)
        val_den = (n - 1.0)*(2.0*n**2 - 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/38.0
            if j > 4:
                break

        return c(n-4.0)*val

    # Generate 1st subdiagonal
    def d_1(n):
        val_num = -2.0*n*(2.0*n**2 + 7.0)
        val_den = (n + 1.0)*(2.0*n**2 + 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -10.0/19.0
            if j > 2:
                break

        return c(n-2.0)*val

    # Generate diagonal
    def d0(n):
        return c(n)*np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))

def tau_c():
    """Compute the chebyshev normalisation c factor for tau boundary"""

    return 2.0

def galerkin_c(n):
    """Compute the chebyshev normalisation c factor for galerkin boundary"""

    val = np.ones(n.shape)

    for i, n in enumerate(n):
        if n == 0:
            val[i] = 0.5
            break
        if n > 0:
            break

    return val

def coeff_iterator(coeffs, pos):
    """Return an iterator over the constants"""

    if coeffs is None:
        it = itertools.cycle([1.0])
    else:
        try:
            if len(coeffs) == (1 + (pos == 0)):
                it = iter(coeffs)
            elif len(coeffs) == 1:
                it = itertools.cycle(coeffs)
            else:
                raise RuntimeError
        except:
            it = itertools.cycle([coeffs])

    return it

def alt_ones(nr, parity):
    """Get array of alternating 1 and -1. Parity is the parity of the -1"""

    if parity == 0:
        return np.cumprod(-np.ones(nr))
    else:
        return -np.cumprod(-np.ones(nr))

