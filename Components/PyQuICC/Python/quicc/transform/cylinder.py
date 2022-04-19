"""Module provides functions to transform chebyshev expansions for the radius in a cylinder between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np


def rgrid(nr):
    """Create the Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,2*nr)+0.5)/(2*nr))

def zgrid(nz):
    """Create the z Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,nz)+0.5)/nz)

def torphys(spec, parity):
    """Transform R spectral coefficients to physical values"""

    n = 2*len(spec)
    full = np.array([0.0]*n)
    full[np.arange(parity,n,2)] = spec

    return fftpack.dct(full,3)

def torcheb(phys, parity):
    """Transform R physical values to spectral coefficients"""

    n = len(phys)
    spec = fftpack.dct(phys,2)/(2*n)

    return spec[np.arange(parity,n,2)]

def tozphys(spec):
    """Transform Z spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def tozcheb(phys):
    """Transform Z physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tophys2d(spec, parity):
    """Transform 2D spectral coefficients to 2D physical values"""

    tmp = spec.copy()

    if spec.dtype == 'complex_':
        phys = 1j*np.zeros((2*spec.shape[0], spec.shape[1]))
        for i in range(spec.shape[0]):
            tmp.real[i,:] = tozphys(spec.real[i,:])
            tmp.imag[i,:] = tozphys(spec.imag[i,:])
        for j in range(phys.shape[1]):
            phys.real[:,j] = torphys(tmp.real[:,j], parity)
            phys.imag[:,j] = torphys(tmp.imag[:,j], parity)
    else:
        phys = np.zeros((2*spec.shape[0], spec.shape[1]))
        for i in range(spec.shape[0]):
            tmp[i,:] = tozphys(spec[i,:])
        for j in range(phys.shape[1]):
            phys[:,j] = torphys(tmp[:,j], parity)

    return phys

def tocheb2d(phys, parity):
    """Transform 2D physical array to 2D spectral coefficients"""

    spec = np.zeros((len(np.arange(parity,phys.shape[0],2)), phys.shape[1]))

    for j in range(phys.shape[1]):
        spec[:,j] = torcheb(phys[:,j], parity)
    for i in range(spec.shape[0]):
        spec[i,:] = tozcheb(spec[i,:])

    return spec
