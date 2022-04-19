"""Module provides functions transforms from spectral to physical space in a cylinder (Worland(r) + Chebyshev(z))."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import scipy.special as sp
import quicc.geometry.worland.wnl as wnl

min_r_points = 150
min_th_points = 200
min_z_points = 150

def nrgrid(nr):
    return max(min_r_points, nr)

def rgrid(nr, m):
    """Create the Chebyshev grid"""

    return wnl.get_grid(nrgrid(nr))

def zgrid(nz):
    """Create the z Chebyshev grid"""

    gN = max(min_z_points, nz)

    return np.cos(np.pi*(np.arange(0,gN)+0.5)/gN)

def eqgrid(m):
    """Create a equatorial (theta) grid for given order"""

    return np.linspace(0, 2*np.pi, max(min_th_points,3*m))

grid_1d = rgrid
grid_fast = rgrid
grid_slow = zgrid

def grid_2d(nr, m, nz):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr, m)
    z = zgrid(nz)
    X, Y = np.meshgrid(r, z)

    return (X, Y)

def grid_fast_per(nr, m):
    """Create 2D grid"""

    phi = eqgrid(m)
    r = rgrid(nr, m)
    rmesh, phimesh = np.meshgrid(r, phi)
    X = rmesh * np.cos(phimesh)
    Y = rmesh * np.sin(phimesh)

    return (X, Y)

def grid_slow_per(nz, m):
    """Compute the 2D grid for the equatorial contours"""

    z = zgrid(nz)
    phi = eqgrid(m)
    X, Y = np.meshgrid(z, phi)

    return (X, Y)

def torphys(spec, m, nr):
    """Transform R spectral coefficients to physical values"""

    mat = []
    n = spec.shape[0]
    for i in range(0, n):
        mat.append(wnl.eval_poly(i, m, nrgrid(nr)))
    mat = np.array(mat)

    return mat.T*np.matrix(spec)

def torspec(phys, m, n):
    """Transform R physical values to spectral coefficients"""

    mat = []
    nr = phys.shape[0]/2
    w = wnl.get_weights(nr)
    for i in range(0, n):
        mat.append(w*wnl.eval_poly(i, m, nr))
    mat = np.matrix(mat)

    return mat*np.matrix(phys).reshape((phys.shape[0], 1))

def tozphys(spec):
    """Transform Z spectral coefficients to physical values"""

    spec = np.array(spec).flatten()
    if len(spec) < min_z_points:
        data = np.hstack((spec, np.zeros(min_z_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data,3)

def tozspec(phys):
    """Transform Z physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tozphys2D(spec):
    """Transform 2D Z spectral coefficients to 2D Z physical values"""

    phys = np.zeros((max(spec.shape[0],min_z_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j] = tozphys(spec[:,j])

    return phys

def toslice(spec, nr, m, nz):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, nz), order = 'F')

    rphys = torphys(spec.real, m, nr) + 1j*torphys(spec.imag, m, nr)

    rspec = np.transpose(rphys)
    phys = tozphys2D(rspec)
    return phys

def tophys2d(spec, m):
    """Transform 2D spectral coefficients to 2D physical values"""

    tmp = spec.copy()

    if spec.dtype == 'complex_':
        phys = 1j*np.zeros((2*spec.shape[0], spec.shape[1]))
        for i in range(spec.shape[0]):
            tmp.real[i,:] = tozphys(spec.real[i,:])
            tmp.imag[i,:] = tozphys(spec.imag[i,:])
        for j in range(phys.shape[1]):
            phys.real[:,j] = torphys(tmp.real[:,j], m)
            phys.imag[:,j] = torphys(tmp.imag[:,j], m)
    else:
        phys = np.zeros((2*spec.shape[0], spec.shape[1]))
        for i in range(spec.shape[0]):
            tmp[i,:] = tozphys(spec[i,:])
        for j in range(phys.shape[1]):
            phys[:,j] = torphys(tmp[:,j], m)

    return phys

def tospec2d(phys, m):
    """Transform 2D physical array to 2D spectral coefficients"""

    spec = np.zeros((len(np.arange(0,phys.shape[0],2)), phys.shape[1]))

    for j in range(phys.shape[1]):
        spec[:,j] = torspec(phys[:,j], 0)
    for i in range(spec.shape[0]):
        spec[i,:] = tozspec(spec[i,:])

    return spec
