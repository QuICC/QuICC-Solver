"""Module provides functions to transform chebyshev expansions for the radius of a shell between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import numpy.polynomial.legendre as leg
import scipy.special as spe

from quicc.transform.spherical import thgrid, totphys, totleg, eqgrid

min_r_points = 2000

def rgrid(nr, a, b):
    """Create the radial Chebyshev grid"""

    gN = max(min_r_points, nr)

    return a*np.cos(np.pi*(np.arange(0,gN)+0.5)/gN) + b

grid_1d = rgrid
grid_fast = rgrid
grid_slow = thgrid

def grid_2d(nr, a, b, maxl, m):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr, a, b)
    th = thgrid(maxl, m)
    rmesh, thmesh = np.meshgrid(r, th)
    X = rmesh * np.sin(thmesh)
    Y = rmesh * np.cos(thmesh)

    return (X, Y)

def grid_fast_per(nr, a, b, m):
    """Compute the 2D grid for the equatorial contours"""

    r = rgrid(nr, a, b)
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

def torphys(spec):
    """Transform radial spectral coefficients to physical values"""

    if len(spec) < min_r_points:
        data = np.hstack((spec, np.zeros(min_r_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data,3)

toprofile = torphys

def torcheb(phys):
    """Transform radial physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def torphys2D(spec):
    """Transform 2D R spectral coefficients to 2D R physical values"""

    phys = np.zeros((max(spec.shape[0],min_r_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j] = torphys(spec[:,j])

    return phys

def torcheb2D(phys):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = torcheb(phys[:,j])

    return phys

def toslice(spec, nr, a, b, maxl, m):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, maxl-m+1), order = 'F')

    rphys = torphys2D(spec.real) + 1j*torphys2D(spec.imag)

    rspec = rphys.T
    phys = totphys(rspec, maxl, m)
    return phys
