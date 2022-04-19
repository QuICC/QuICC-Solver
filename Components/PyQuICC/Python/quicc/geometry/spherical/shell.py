"""Module provides functions to generate sparse operators in a spherical shell using a spherical harmonics expansion in the angular directions."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import quicc.geometry.spherical.shell_radius as rad
import quicc.geometry.spherical.shell_sh as sh
import quicc.geometry.spherical.shell_boundary as sphbc


def convert_bc(bc):
    """Convert boundary condition to be suitable for kronecker product boundaries"""

    if bc[0] < 0:
        bcr = bc
    else:
        bcr = rad.radbc.no_bc()
        for k, v in bc.items():
            if k != 0:
                bcr[k] = v

    return bcr

def sh_coeff(coeff):
    """Compute coefficients from spherical harmonic expansion"""

    if coeff == 'laplh':
        def fct(x):
            return x*(x + 1.0)
    elif coeff == 'invlaplh':
        def fct(x):
            if x == 0:
                return 0
            else:
                return 1.0/(x*(x + 1.0))
    else:
        def fct(x):
            return 1.0

    return fct

def fix_l_zero(nr, ri, ro, m, mat, bc, fix):
    """Fix problems with unused l = 0 modes"""

    if m > 0 or not fix:
        return mat
    elif fix == 'zero':
        return rad.zblk(nr, ri, ro, bc)
    elif fix == 'set':
        if bc[0] < 0:
            nn = nr - bc['rt']
        else:
            nn = nr
        return rad.qid(nn, ri, ro, 0, {0:0})
    else:
        raise RuntimeError("Unkown l=0 fix!")

def make_sh_zero(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Generic function to create a coupled spherical harmonics operator"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    base_mat = rad.zblk(nr, ri, ro, rad.radbc.no_bc())
    if restriction is None or m in restriction:
        bcr = sphbc.ldependent_bc(bcr, m)
        mat = coeff*shc(m)*rad.radbc.constrain(base_mat, ri, ro, bcr)
        mat = fix_l_zero(nr, ri, ro, m, mat, bcr, l_zero_fix)
    else:
        mat = rad.zblk(nr, ri, ro, bcr)
    blocks = [mat]

    for l in range(m+1, maxnl):
        if restriction is None or l in restriction:
            bcr = sphbc.ldependent_bc(bcr, l)
            mat = coeff*shc(l)*rad.radbc.constrain(base_mat, ri, ro, bcr)
        else:
            mat = rad.zblk(nr, ri, ro, bcr)
        blocks.append(mat)
    mat = spsp.block_diag(blocks, format = 'coo')

    return sphbc.constrain(mat, nr, ri, ro, maxnl, m, bc, l_zero_fix, restriction = restriction)

def make_sh_operator(op, nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Generic function to create a coupled spherical harmonics operator"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    base_mat = op(nr, ri, ro, rad.radbc.no_bc())
    if restriction is None or m in restriction:
        bcr = sphbc.ldependent_bc(bcr, m)
        mat = coeff*shc(m)*rad.radbc.constrain(base_mat, ri, ro, bcr)
        mat = fix_l_zero(nr, ri, ro, m, mat, bcr, l_zero_fix)
    else:
        mat = rad.zblk(nr, ri, ro, bcr)
    blocks = [mat]

    for l in range(m+1, maxnl):
        if restriction is None or l in restriction:
            bcr = sphbc.ldependent_bc(bcr, l)
            mat = coeff*shc(l)*rad.radbc.constrain(base_mat, ri, ro, bcr)
        else:
            mat = rad.zblk(nr, ri, ro, bcr)
        blocks.append(mat)
    mat = spsp.block_diag(blocks, format = 'coo')

    return sphbc.constrain(mat, nr, ri, ro, maxnl, m, bc, l_zero_fix, restriction = restriction)

def make_sh_loperator(op, nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Generic function to create a coupled l dependent spherical harmonics operator"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    if restriction is None or m in restriction:
        bcr = sphbc.ldependent_bc(bcr, m)
        mat = coeff*shc(m)*op(nr, ri, ro, m, bcr)
        mat = fix_l_zero(nr, ri, ro, m, mat, bcr, l_zero_fix)
    else:
        mat = rad.zblk(nr, ri, ro, bcr)
    blocks = [mat]

    for l in range(m+1, maxnl):
        if restriction is None or l in restriction:
            bcr = sphbc.ldependent_bc(bcr, l)
            mat = coeff*shc(l)*op(nr, ri, ro, l, bcr)
        else:
            mat = rad.zblk(nr, ri, ro, bcr)
        blocks.append(mat)
    mat = spsp.block_diag(blocks, format = 'coo')

    return sphbc.constrain(mat, nr, ri, ro, maxnl, m, bc, l_zero_fix, restriction = restriction)

def make_sh_qoperator(opl, opr, nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create the coupled operator for the coriolis Q term"""

    # Only compute if there are at least 2 harmonic degrees
    if maxnl - m > 1:
        cor_r = sh.coriolis_r(maxnl, m)
        cordr = sh.coriolisdr(maxnl, m)

        if restriction is not None:
            cor_r = cor_r.tolil()
            cordr = cordr.tolil()
            res_cols = list(set(range(0,cor_r.shape[1])) - set([i - m for i in restriction]))
            cor_r[:,res_cols] = 0
            cordr[:,res_cols] = 0

        cor_r = cor_r.tocsr()
        cordr = cordr.tocsr()

        bcr = convert_bc(bc)
        shc = sh_coeff(with_sh_coeff)

        base_opl = opl(nr, ri, ro, rad.radbc.no_bc())
        base_opr = opr(nr, ri, ro, rad.radbc.no_bc())

        assert(l_zero_fix == 'zero')
        bcr = sphbc.ldependent_bc(bcr, m)
        rmatl = rad.radbc.constrain(base_opl, ri, ro, bcr)
        rmatr = rad.radbc.constrain(base_opr, ri, ro, bcr)
        rmatl = fix_l_zero(nr, ri, ro, m, rmatl, bcr, l_zero_fix)
        rmatr = fix_l_zero(nr, ri, ro, m, rmatr, bcr, l_zero_fix)
        mat = coeff*shc(m)*spsp.kron(cor_r[0,:],rmatl, format = 'coo') + coeff*shc(m)*spsp.kron(cordr[0,:], rmatr, format = 'coo')
        rows = [mat]
        for ir,l in enumerate(range(m+1, maxnl)):
            bcr = sphbc.ldependent_bc(bcr, l)
            rmatl = rad.radbc.constrain(base_opl, ri, ro, bcr)
            rmatr = rad.radbc.constrain(base_opr, ri, ro, bcr)
            mat = coeff*shc(l)*spsp.kron(cor_r[ir+1,:],rmatl, format = 'coo') + coeff*shc(l)*spsp.kron(cordr[ir+1,:],rmatr, format = 'coo')
            rows.append(mat)
        mat = spsp.vstack(rows)
    else:
        mat = rad.zblk(nr, ri, ro, bc)

    return sphbc.constrain(mat, nr, ri, ro, maxnl, m, bc, l_zero_fix, restriction = restriction)

def zblk(nr, ri, ro, maxnl, m, bc, l_zero_fix = False, restriction = None):
    """Create a block of zeros"""

    return make_sh_zero(nr, ri, ro, maxnl, m, bc, 0.0, with_sh_coeff = None, l_zero_fix = l_zero_fix, restriction = restriction)

def i2(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2 radial operator kronecker with an identity"""

    return make_sh_operator(rad.i2, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i2r2(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2r2 radial operator kronecker with an identity"""

    return make_sh_operator(rad.i2r2, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i2r3(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2r3 radial operator kronecker with an identity"""

    return make_sh_operator(rad.i2r3, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i2r2coriolis(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2r2 radial operator kronecker with coriolis Q term"""

    return make_sh_qoperator(rad.i2r1, rad.i2r2d1, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i2r2lapl(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2r2lapl radial operator kronecker with an identity"""

    return make_sh_loperator(rad.i2r2lapl, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i2r3lapl(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i2r3lapl radial operator kronecker with an identity"""

    return make_sh_loperator(rad.i2r3lapl, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i4r4(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i4r4 radial operator kronecker with an identity"""

    return make_sh_operator(rad.i4r4, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i4r4coriolis(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i4r4 radial operator kronecker with coriolis Q term"""

    return make_sh_qoperator(rad.i4r3, rad.i4r4d1, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i4r4lapl(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i4r4lapl radial operator kronecker with an identity"""

    return make_sh_loperator(rad.i4r4lapl, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def i4r4lapl2(nr, ri, ro, maxnl, m, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False, restriction = None):
    """Create a i4r4lapl2 radial operator kronecker with an identity"""

    return make_sh_loperator(rad.i4r4lapl2, nr, ri, ro, maxnl, m, bc, coeff, with_sh_coeff = with_sh_coeff, l_zero_fix = l_zero_fix, restriction = restriction)

def qid(nr, ri, ro, maxnl, m, qr, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr = convert_bc(bc)

    nl = maxnl - m
    mat = coeff*spsp.kron(rad.qid(nl, ri, ro, 0, rad.radbc.no_bc()), rad.qid(nr, ri, ro, qr, bcr), format = 'coo')
    return sphbc.constrain(mat, nr, ri, ro, maxnl, m, bc)

def stencil(nr, ri, ro, maxnl, m, bc, make_square):
    """Create a galerkin stencil matrix"""

    bcr = convert_bc(bc)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = rad.stencil(nr, ri, ro, bcr, make_square)

    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,rad.stencil(nr, ri, ro, bcr, make_square)))

    return mat
