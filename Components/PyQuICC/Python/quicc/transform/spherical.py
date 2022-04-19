"""Module provides functions to transform spherical harmonics expansions between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.polynomial.legendre as leg
import scipy.special as spe

min_th_points = 200
min_phi_points = 200

def txgrid(maxl, m):
    """Create the theta Gauss-Legendre grid"""

    nt = max(min_th_points,3*(maxl - m + 1)//2)
    nt = nt + (nt+1)%2
    x, tmp = leg.leggauss(nt)

    return x

def thgrid(maxl, m):
    """Convert legendre grid to theta angle"""

    return np.arccos(txgrid(maxl, m))

def totphys(spec, maxl, m):
    """Tranform theta spectral coefficients to physical values"""

    mat = plm(maxl, m)
    phys = mat.dot(spec)

    return phys

def totleg(phys, maxl, m):
    """Tranform theta physical values to spectral coefficients"""

    nt = 3*(maxl - m + 1)//2
    x, w = leg.leggauss(nt)
    mat = plm(maxl, m).T
    spec = mat.dot(np.diag(w).dot(phys))

    return spec

def eqgrid(m, phi = 2*np.pi):
    """Create a equatorial (phi) grid for given harmonic order"""

    return np.linspace(0, phi, max(min_phi_points,3*m))

def plm(maxl, m):
    """Compute the normalized associated legendre polynomial projection matrix"""

    x = txgrid(maxl, m)
    mat = np.zeros((len(x), maxl - m + 1))
    mat[:,0] = pmm(maxl, m)[:,0]
    mat[:,1] = np.sqrt(2.0*m + 3.0)*x*mat[:,0]
    for i, l in enumerate(range(m+2, maxl + 1)):
        mat[:,i+2] = np.sqrt((2.0*l + 1)/(l - m))*np.sqrt((2.0*l - 1.0)/(l + m))*x*mat[:,i+1] - np.sqrt((2.0*l + 1)/(2.0*l - 3.0))*np.sqrt((l + m - 1.0)/(l + m))*np.sqrt((l - m - 1.0)/(l - m))*mat[:,i]

    return mat

def pmm(maxl, m):
    """Compute the normalized associated legendre polynomial of order and degree m"""

    x = txgrid(maxl, m)
    mat = np.zeros((len(x), 1))
    mat[:,0] = 1.0/np.sqrt(2.0)
    sx = np.sqrt(1.0 - x**2)
    for i in range(1, m+1):
        mat[:,0] = -np.sqrt((2.0*i + 1.0)/(2.0*i))*sx*mat[:,0]

    return mat
