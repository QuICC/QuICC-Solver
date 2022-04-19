"""Module provides functions to transform chebyshev expansions for the radius of a sphere between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import functools

from quicc.transform.spherical import thgrid, totphys, totleg, eqgrid as eqgrid_full

eqgrid = functools.partial(eqgrid_full, phi = np.pi)

min_r_points = 500

def rgrid(nr):
    """Create the radial Chebyshev grid"""

    gN = max(min_r_points, 2*nr)

    return np.cos(np.pi*(np.arange(0,gN)+0.5)/gN)

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

def torphys(spec, parity):
    """Transform R spectral coefficients to R physical values"""

    n = 2*len(spec)
    full = np.array([0.0]*n)
    full[np.arange(parity,n,2)] = spec

    if len(full) < min_r_points:
        data = np.hstack((full, np.zeros(min_r_points - len(full))))
    else:
        data = full

    return fftpack.dct(data,3)

toprofile = torphys

def torcheb(phys, parity):
    """Transform R physical values to R spectral coefficients"""

    n = len(phys)
    spec = fftpack.dct(phys,2)/(2*n)

    return spec[np.arange(parity,n,2)]

def torphys2D(spec, m):
    """Transform 2D R spectral coefficients to 2D R physical values"""

    phys = np.zeros((max(2*spec.shape[0],min_r_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j] = torphys(spec[:,j], (m+j)%2)

    return phys

def torcheb2D(phys, m):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = torcheb(phys[:,j], (m+j)%2)

    return phys

def toslice(spec, nr, maxl, m):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, maxl-m+1), order = 'F')

    rphys = torphys2D(spec.real, m) + 1j*torphys2D(spec.imag, m)

    rspec = np.transpose(rphys)
    phys = totphys(rspec, maxl, m)
    return phys
