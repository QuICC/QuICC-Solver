"""Module provides functions to transform Worland expansions for the radius of a sphere between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import functools
import scipy.special as sp
import quicc.geometry.worland.wnl as wnl

from quicc.transform.spherical import thgrid, totphys, totleg, eqgrid as eqgrid_full

eqgrid = functools.partial(eqgrid_full, phi = 2*np.pi)

min_r_points = 500

def nrgrid(nn):
    return max(min_r_points, 2*nn)

def rgrid(nr):
    """Create the radial Chebyshev grid"""

    return wnl.get_grid(nrgrid(nr))

grid_1d = rgrid
grid_fast = rgrid
grid_slow = thgrid

def grid_2d(nr, maxl, m):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr)
    th = thgrid(maxl, m)
    rmesh, thmesh = np.meshgrid(r, th)
    X = rmesh * np.sin(thmesh)
    Y = rmesh * np.cos(thmesh)

    return (X, Y)

def grid_fast_per(nr, m):
    """Compute the 2D grid for the equatorial contours"""

    r = rgrid(nr)
    phi = eqgrid(m)
    rmesh, phimesh = np.meshgrid(r, phi)
    X = rmesh * np.cos(phimesh)
    Y = rmesh * np.sin(phimesh)

    return (X, Y)

def grid_slow_per(maxl, tm, m):
    """Compute the 2D grid for the equatorial contours"""

    th = thgrid(maxl, tm)
    phi = eqgrid(m)
    X, Y = np.meshgrid(th, phi)

    return (X, Y)

def torphys(spec, l, nr):
    """Transform R spectral coefficients to R physical values"""

    mat = []
    n = spec.shape[0]
    for i in range(0, n):
        mat.append(wnl.eval_poly(i, l, nrgrid(nr)))
    mat = np.matrix(mat)

    return mat.T*np.matrix(spec)

toprofile = torphys

def torspec(phys, l, n):
    """Transform R physical values to R spectral coefficients"""

    mat = []
    nr = phys.shape[0]
    w = wnl.get_weights(nr)
    for i in range(0, n):
        mat.append(w*wnl.eval_poly(i, l, nr))
    mat = np.matrix(mat)

    return mat*np.matrix(phys).reshape((phys.shape[0], 1))

def torphys2D(spec, m, nr):
    """Transform 2D R spectral coefficients to 2D R physical values"""

    phys = np.zeros((nrgrid(nr), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j:j+1] = torphys(spec[:,j:j+1], m+j, nr)

    return phys

def torspec2D(phys, m):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = torspec(phys[:,j], m+j)

    return phys

def toslice(spec, nr, maxl, m):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, maxl-m+1), order = 'F')

    rphys = torphys2D(spec.real, m, nr) + 1j*torphys2D(spec.imag, m, nr)

    rspec = np.transpose(rphys)
    phys = totphys(rspec, maxl, m)
    return phys
