"""Module provides functions to generate the radial boundary conditions in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import quicc.base.utils as utils


default_use_parity = False

def linear_map(ri, ro):
    """Calculat a and b for linear map r = a*x + b"""

    a = (ro - ri)/2.0;
    b = (ro + ri)/2.0;

    return (a, b)

def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, ri, ro, bc, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, ri, ro, bc, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, ri, ro, bc)
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

def apply_tau(mat, ri, ro, bc, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[1], ri, ro, 1, bc)
    elif bc[0] == 11:
        cond = tau_value(mat.shape[1], ri, ro, -1, bc)
    elif bc[0] == 12:
        cond = tau_diff(mat.shape[1], ri, ro, 1, bc)
    elif bc[0] == 13:
        cond = tau_diff(mat.shape[1], ri, ro, -1, bc)
    elif bc[0] == 20:
        cond = tau_value(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 21:
        cond = tau_diff(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 22:
        cond = tau_rdiffdivr(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 23:
        cond = tau_insulating(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 24:
        cond = tau_couette(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 25:
        cond1 = tau_diff(mat.shape[1], ri, ro, -1, bc)
        cond2 = tau_value(mat.shape[1], ri, ro, 1, bc)
        cond = np.vstack([cond1, cond2])
    elif bc[0] == 26:
        cond1 = tau_value(mat.shape[1], ri, ro, -1, bc)
        cond2 = tau_diff(mat.shape[1], ri, ro, 1, bc)
        cond = np.vstack([cond1, cond2])
    elif bc[0] == 40:
        cond = tau_value_diff(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 41:
        cond = tau_value_diff2(mat.shape[1], ri, ro, 0, bc)
    elif bc[0] == 42:
        cond1 = tau_value_diff2(mat.shape[1], ri, ro, -1, bc)
        cond2 = tau_value_diff(mat.shape[1], ri, ro, 1, bc)
        cond = np.vstack([cond1, cond2])
    elif bc[0] == 43:
        cond1 = tau_value_diff2(mat.shape[1], ri, ro, 1, bc)
        cond2 = tau_value_diff(mat.shape[1], ri, ro, -1, bc)
        cond = np.vstack([cond1, cond2])
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

def tau_value(nr, ri, ro, pos, bc):
    """Create the boundary value tau line(s)"""

    it = coeff_iterator(bc.get('c', None), pos)

    cond = []
    c = next(it)
    if pos >= 0:
        cnst = c*tau_c()
        cond.append(cnst*np.ones(nr))
        cond[-1][0] /= tau_c()
        c = next(it)

    if pos <= 0:
        cnst = c*tau_c()
        cond.append(cnst*alt_ones(nr, 1))
        cond[-1][0] /= tau_c()

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_diff(nr, ri, ro, pos, bc):
    """Create the first derivative tau line(s)"""

    it = coeff_iterator(bc.get('c', None), pos)
    a,b = linear_map(ri, ro)

    cond = []
    c = next(it)
    ns = np.arange(0,nr)
    if pos >= 0:
        cnst = c*(tau_c()/a)
        cond.append(cnst*ns**2)
        c = next(it)

    if pos <= 0:
        cnst = c*(tau_c()/a)
        cond.append(cnst*ns**2*alt_ones(nr, 0))

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_diff2(nr, ri, ro, pos, bc):
    """Create the second deriviative tau line(s)"""

    it = coeff_iterator(bc.get('c', None), pos)
    a,b = linear_map(ri, ro)

    cond = []
    c = next(it)
    ns = np.arange(0,nr)
    if pos >= 0:
        cnst = c*(tau_c()/a**2)
        cond.append((cnst/3.0)*(ns**4 - ns**2))
        c = next(it)

    if pos <= 0:
        cnst = c*(tau_c()/a**2)
        cond.append((cnst/3.0)*(ns**4 - ns**2)*alt_ones(nr, 1))

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_rdiffdivr(nr, ri, ro, pos, bc):
    """Create the r D 1/r tau line(s)"""

    it = coeff_iterator(bc.get('c', None), pos)
    a,b = linear_map(ri, ro)

    cond = []
    c = next(it)
    ns = np.arange(0,nr)
    if pos >= 0:
        cond.append(c*((1.0/a)*ns**2 - (1.0/(a+b)))*tau_c())
        cond[-1][0] /= tau_c()
        c = next(it)

    if pos <= 0:
        cond.append(c*((1.0/a)*ns**2 + (1.0/(-a+b)))*tau_c()*alt_ones(nr, 0))
        cond[-1][0] /= tau_c()

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_insulating(nr, ri, ro, pos, bc):
    """Create the insulating boundray tau line(s)"""

    assert(bc.get('l', None) is not None)

    it = coeff_iterator(bc.get('c', None), pos)
    a,b = linear_map(ri, ro)
    l = bc['l']

    cond = []
    c = next(it)
    ns = np.arange(0,nr)
    if pos >= 0:
        cond.append(c*((2.0/a)*ns**2 + ((l+1.0)/(a+b))*tau_c()))
        cond[-1][0] /= tau_c()
        c = next(it)

    if pos <= 0:
        cond.append(c*((2.0/a)*ns**2 + (l/(-a+b))*tau_c())*alt_ones(nr, 0))
        cond[-1][0] /= tau_c()

    if bc.get('use_parity', default_use_parity) and pos == 0:
        t = cond[0].copy()
        cond[0] = (cond[0] + cond[1])/2.0
        cond[1] = (t - cond[1])/2.0

    return np.array(cond)

def tau_couette(nr, ri, ro, pos, bc):
    """Create the toroidal Couette boundray tau line(s)"""

    assert(bc.get('c', None) is not None)
    #TODO: think of the ordering
    #assert(bc.get('l', None) is not None)
    assert(pos == 0)

    return tau_value(nr, ri, ro, 0, None)

def tau_value_diff(nr, ri, ro, pos, bc):
    """Create the no penetration and no-slip tau line(s)"""

    itv = coeff_iterator(bc.get('cv', None), pos)
    itd = coeff_iterator(bc.get('cd', None), pos)

    cond = []
    bcv = dict()
    bcd = dict()
    bcv['c'] = next(itv)
    bcd['c'] = next(itd)
    if pos >= 0:
        cond.append(tau_value(nr, ri, ro, 1, bcv)[0])
        cond.append(tau_diff(nr, ri, ro, 1, bcd)[0])
        bcv['c'] = next(itv)
        bcd['c'] = next(itd)

    if pos <= 0:
        cond.append(tau_value(nr, ri, ro, -1, bcv)[0])
        cond.append(tau_diff(nr, ri, ro, -1, bcd)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    return np.array(cond)

def tau_value_diff2(nr, ri, ro, pos, bc):
    """Create the no penetration and stress-free tau line(s)"""

    itv = coeff_iterator(bc.get('cv', None), pos)
    itd2 = coeff_iterator(bc.get('cd2', None), pos)

    cond = []
    bcv = dict()
    bcd2 = dict()
    bcv['c'] = next(itv)
    bcd2['c'] = next(itd2)
    if pos >= 0:
        cond.append(tau_value(nr, ri, ro, 1, bcv)[0])
        cond.append(tau_diff2(nr, ri, ro, 1, bcd2)[0])
        bcv['c'] = next(itv)
        bcd2['c'] = next(itd2)

    if pos <= 0:
        cond.append(tau_value(nr, ri, ro, -1, bcv)[0])
        cond.append(tau_diff2(nr, ri, ro, -1, bcd2)[0])

    if bc.get('use_parity', default_use_parity) and pos == 0:
        tv = cond[0].copy()
        td = cond[1].copy()
        cond[0] = (cond[0] + cond[2])/2.0
        cond[1] = (cond[1] + cond[3])/2.0
        cond[2] = (tv - cond[2])/2.0
        cond[3] = (td - cond[3])/2.0

    return np.array(cond)

def tau_last(nr, nrow):
    """Create the last modes to zero tau line(s)"""

    cond = np.zeros((nrow, nr))
    for j in range(0, nrow):
        cond[j,nr-nrow+j] = tau_c()

    return cond

def stencil(nr, ri, ro, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nr, ri, ro, 1, bc.get('c',None))
    elif bc[0] == -11:
        mat = stencil_value(nr, ri, ro, -1, bc.get('c',None))
    elif bc[0] == -12:
        mat = stencil_value(nr, ri, ro, 1, bc.get('c',None))
    elif bc[0] == -13:
        mat = stencil_value(nr, ri, ro, -1, bc.get('c',None))
    elif bc[0] == -20:
        mat = stencil_value(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] == -21:
        mat = stencil_diff(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] == -22:
        mat = stencil_rdiffdivr(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] == -23:
        mat = stencil_insulating(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] == -40:
        mat = stencil_value_diff(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] == -41:
        mat = stencil_value_diff2(nr, ri, ro, 0, bc.get('c',None))
    elif bc[0] < -1 and bc[0] > -5:
        mat = restrict_eye(nr, 'cr', -bc[0])

    return mat

def apply_galerkin(mat, ri, ro, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = mat*stencil(nr, ri, ro, bc)
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

def stencil_value(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero boundary value"""

    assert(coeffs.get('c', None) is None)

    ns = np.arange(0,nr,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos

    # Generate subdiagonal
    def d_1(n):
        return galerkin_c(n+offsets[0])*sgn

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero 1st derivative"""

    assert(coeffs.get('c', None) is None)

    ns = np.arange(0,nr,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos

    # Generate subdiagonal
    def d_1(n):
        return sgn*(n+offsets[0])**2/n**2

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_rdiffdivr(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero r D 1/r derivative"""

    assert(coeffs.get('c', None) is None)
    assert(pos == 0)

    a,b = linear_map(ri, ro)

    ns = np.arange(0,nr,1)
    offsets = [-2, -1, 0]

    # Generate 2nd subdiagonal
    def d_2(n):
        #val_num = a**2*((n - 3.0)**2*n**2 - 2.0) - b**2*(n**2 - 3.0*n + 2.0)**2
        #val_den = -a**2*((n - 2.0)*n - 1.0)*(n**2 - 2.0) + b**2*(n - 1.0)**2*n**2
        val_num = (n - 2.0)*(b**2*(n - 2.0)*(n - 1.0) - a**2*(n - 3.0)*n)
        val_den = n*(a**2*(n - 2.0)*(n + 1.0) - b**2*(n - 1.0)*n)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 2:
                #val[i] = a**2/(2.0*a**2 + 4.0*b**2)
                val[i] = 0
            if j > 2:
                break

        return val

    # Generate 1st subdiagonal
    def d_1(n):
        #val_num = -8.0*a*b*n
        #val_den = a**2*(n**2 - 2.0)*(n*(n + 2.0) - 1.0) - b**2*n**2*(n + 1.0)**2
        val_num = -4.0*a*b
        val_den = (n + 1.0)*(a**2*(n**2 + n - 2.0) - b**2*n*(n + 1.0))
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 1:
                #val[i] = 2.0*a*b/(a**2 + 2.0*b**2)
                val[i] = a/(2.0*b)
            if j > 1:
                break

        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_insulating(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for an insulating boundary"""

    assert(coeffs.get('c', None) is None)
    assert(coeffs.get('l', None) is not None)
    assert(pos == 0)

    a,b = linear_map(ri, ro)
    l = coeffs['l']

    ns = np.arange(0,nr,1)
    offsets = [-2, -1, 0]

    # Generate 2nd subdiagonal
    def d_2(n):
        val_num = a**2*(2.0*l*(l+1.0)-2.0*(n-3.0)*n*((n-3.0)*n+5.0)-13.0)+a*b*(2.0*l+1.0)*(2.0*(n-3.0)*n+5.0)+2.0*b**2*(n**2-3.0*n+2.0)**2
        val_den = a**2*(2.0*l*(l+1.0)-2.0*(n-1.0)*n*((n-1.0)*n+1.0)-1.0)+a*b*(2.0*l+1.0)*(2.0*(n-1.0)*n+1.0)+2.0*b**2*(n-1.0)**2.0*n**2
        val = -val_num/val_den
        for i,j in enumerate(n):
            if j == 2:
                corr_num = a*(a*(2.0*l**2+2.0*l-1.0)+2.0*b*l+b)
                corr_den = 2.0*(a**2*(2.0*l**2+2.0*l-13.0)+5.0*a*(2.0*b*l+b)+8.0*b**2)
                val[i] = -corr_num/corr_den
            if j > 2:
                break

        return val

    # Generate 1st subdiagonal
    def d_1(n):
        val_num = 4.0*a*n*((2.0*l+1.0)*a - b)
        val_den = a**2*(2.0*l*(l+1.0)-2.0*n*(n+1.0)*(n**2+n+1.0)-1.0)+a*b*(2.0*l+1.0)*(2.0*n*(n+1.0)+1.0)+2.0*b**2*n**2*(n+1.0)**2
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 1:
                corr_num = 2.0*a*(2.0*a*l+a-b)
                corr_den = a**2*(2.0*l**2+2.0*l-13.0)+5.0*a*(2.0*b*l+b)+8.0*b**2
                val[i] = corr_num/corr_den
            if j > 1:
                break

        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff2(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero 2nd derivative"""

    assert(coeffs.get('c', None) is None)

    ns = np.arange(0,nr,1)
    if pos == 0:
        offsets = [-2, 0]

        # Generate subdiagonal
        def d_1(n):
            return -(n - 3.0)*(n - 2.0)**2/(n**2*(n + 1.0))

    else:
        offsets = [-1, 0]

        # Generate subdiagonal
        def d_1(n):
            return -pos*(n - 2.0)*(n - 1.0)/(n*(n + 1.0))

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags, offsets = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero boundary value and a zero 1st derivative"""

    assert(coeffs.get('c', None) is None)
    assert(pos == 0)

    ns = np.arange(0,nr,1)
    offsets = [-4, -2, 0]

    # Generate 2nd subdiagonal
    def d_2(n):
        val = (n - 3.0)/(n - 1.0)
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/6.0
            if j > 4:
                break

        return val

    # Generate 1st subdiagonal
    def d_1(n):
        val = -2.0*n/(n + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -2.0/3.0
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

def stencil_value_diff2(nr, ri, ro, pos, coeffs = None):
    """Create stencil matrix for a zero boundary value and a zero 2nd derivative"""

    assert(coeffs.get('c', None) is None)
    assert(pos == 0)

    ns = np.arange(0,nr,1)
    offsets = [-4, -2, 0]

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

        return val

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
