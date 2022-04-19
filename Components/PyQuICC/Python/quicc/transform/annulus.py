"""Module provides functions to transform chebyshev expansions for the radius of an annulus between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np

min_r_points = 2000
min_th_points = 2000
min_z_points = 2000

def rgrid(nr, a, b):
    """Create the radial Chebyshev grid"""

    gN = max(min_r_points, nr)

    return a*np.cos(np.pi*(np.arange(0,gN)+0.5)/gN) + b

def eqgrid(m):
    """Create a equatorial (theta) grid for given order"""

    return np.linspace(0, 2*np.pi, max(min_th_points,3*m))

def zgrid(nz):
    """Create the z Chebyshev grid"""

    gN = max(min_z_points, nz)

    return np.cos(np.pi*(np.arange(0,gN)+0.5)/gN)

def grid_2d(nr, a, b, nz):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr, a, b)
    z = zgrid(nz)
    X, Y = np.meshgrid(r, z)

    return (X, Y)

def grid_eq(nr, a, b, m):
    """Compute the 2D grid for the equatorial contours"""

    r = rgrid(nr, a, b)
    th = eqgrid(m)
    rmesh, thmesh = np.meshgrid(r, th)
    X = rmesh * np.cos(thmesh)
    Y = rmesh * np.sin(thmesh)

    return (X, Y)

def torphys(spec):
    """Transform R spectral coefficients to physical values"""

    if len(spec) < min_r_points:
        data = np.hstack((spec, np.zeros(min_r_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data,3)

def torcheb(phys):
    """Transform R physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tozphys(spec):
    """Transform Z spectral coefficients to physical values"""

    if len(spec) < min_z_points:
        data = np.hstack((spec, np.zeros(min_z_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data,3)

def tozcheb(phys):
    """Transform Z physical values to spectral coefficients"""

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

def tozphys2D(spec):
    """Transform 2D Z spectral coefficients to 2D Z physical values"""

    phys = np.zeros((max(spec.shape[0],min_z_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j] = tozphys(spec[:,j])

    return phys

def tozcheb2D(phys):
    """Transform 2D Z physical values to 2D Z spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = tozcheb(phys[:,j])

    return phys

def toslice(spec, nr, nz):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, nz), order = 'F')

    rphys = torphys2D(spec.real) + 1j*torphys2D(spec.imag)

    rspec = np.transpose(rphys)
    phys = tozphys2D(rspec)
    return phys

def tophys2d(spec):
    """Transform 2D spectral coefficients to 2D physical values"""

    phys = spec.copy()

    if spec.dtype == 'complex_':
        for i in range(spec.shape[0]):
            phys.real[i,:] = tozphys(spec.real[i,:])
            phys.imag[i,:] = tozphys(spec.imag[i,:])
        for j in range(spec.shape[1]):
            phys.real[:,j] = torphys(phys.real[:,j])
            phys.imag[:,j] = torphys(phys.imag[:,j])
    else:
        for i in range(spec.shape[0]):
            phys[i,:] = tozphys(spec[i,:])
        for j in range(spec.shape[1]):
            phys[:,j] = torphys(phys[:,j])

    return phys

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    spec = phys.copy()

    if spec.dtype == 'complex_':
        for j in range(phys.shape[1]):
            spec.real[:,j] = torcheb(phys.real[:,j])
            spec.imag[:,j] = torcheb(phys.imag[:,j])
        for i in range(phys.shape[0]):
            spec.real[i,:] = tozcheb(spec.real[i,:])
            spec.imag[i,:] = tozcheb(spec.imag[i,:])
    else:
        for j in range(phys.shape[1]):
            spec[:,j] = torcheb(phys[:,j])
        for i in range(phys.shape[0]):
            spec[i,:] = tozcheb(spec[i,:])

    return spec
