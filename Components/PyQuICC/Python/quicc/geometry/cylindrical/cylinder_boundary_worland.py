"""Module provides functions to generate the boundary conditions in a cylinder with Worland expansion in radius"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import quicc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import quicc.geometry.cylindrical.cylinder_radius_boundary_worland as radbc
import quicc.base.utils as utils


def no_bc():
    """Get a no boundary condition flag"""

    return {'r':radbc.no_bc(), 'z':c1dbc.no_bc()}

def brid(n, m, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(n-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = radbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return radbc.constrain(mat, m, tbc)

def bzid(nz, zi, zo, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(nz-q, nz-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(nz-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = c1dbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return c1dbc.constrain(mat, zi, zo, tbc)

def convert_priority(priority, qr, qz):
    """Convert priority flag into shifts"""

    sr = 0
    dr = 0
    sz = 0
    dz = 0
    if priority == 'r':
        sr = qr
    elif priority == 'z':
        sz = qz
    elif priority == 'n':
        sr = qr
        sz = qz
    elif priority == 'sr':
        sr = qr
        dr = -qr
    elif priority == 'sz':
        sz = qz
        dz = -qz
    else:
        raise RuntimeError("Unknown boundary condition priority!")

    return (sr, dr, sz, dz)

def constrain(mat, nr, m, nz, zi, zo, qr, qz, bc, location = 't', restriction = None):
    """Contrain the matrix with the Tau boundary condition"""

    sr, dr, sz, dz = convert_priority(bc.get('priority', 'r'), qr, qz)

    bc_mat = mat
    if bc['r'][0] > 0 or bc['r'].get('mixed',None) is not None:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, m, bc['r'], location = location)
        bc_mat = bc_mat + utils.restricted_kron_2d(bzid(nz, zi, zo, sz, dz, bc['z'], location = location), bcMat, restriction = restriction)
        # Implement mixed boundary conditions
        if bc['r'].get('mixed', None) is not None:
            mix = bc['r']['mixed']
            pad = mix.get('pad',0)
            s = mix.get('kron_shift',0)
            for mbc, mc, mkron in mixed_iterator(mix):
                bcMat = spsp.lil_matrix((nr,nr))
                bcMat = radbc.constrain(bcMat, m, {0:mbc, 'c':mc}, pad_zeros = pad, location = location)
                bc_mat += utils.restricted_kron_2d(bzid(nz, zi, zo, sz, dz, bc['z'], location = location)*mkron(nz+s, zi, zo, {0:0, 'rt':s, 'cr':s}), bcMat, restriction = restriction)

    if bc['z'][0] > 0 or bc['z'].get('mixed',None) is not None:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, zi, zo, bc['z'], location = location)
        bc_mat = bc_mat + utils.restricted_kron_2d(bcMat, brid(nr,m,sr,dr,bc['r'], location = location), restriction = restriction)
        # Implement mixed boundary conditions
        if bc['z'].get('mixed', None) is not None:
            mix = bc['z']['mixed']
            pad = mix.get('pad',0)
            s = mix.get('kron_shift',0)
            for mbc, mc, mkron in mixed_iterator(mix):
                bcMat = spsp.lil_matrix((nz,nz))
                bcMat = c1dbc.constrain(bcMat, zi, zo, mix, pad_zeros = pad, location = location)
                bc_mat += utils.restricted_kron_2d(bcMat, brid(nr,m,sr,dr,bc['r'], location = location)*mkron(nr+s, m, {0:0, 'rt':s, 'cr':s}), restriction = restriction)

    return bc_mat

def mixed_iterator(mixed):
    """Return an iterator over the constants"""

    try:
        if len(mixed[0]) == len(mixed['kron']) and len(mixed[0]) == len(mixed['c']):
            it = itertools.izip(mixed[0], mixed['c'], mixed['kron'])
        else:
            raise RuntimeError
    except:
        it = itertools.izip([mixed[0]], [mixed.get('c',1)], [mixed['kron']])

    return it
