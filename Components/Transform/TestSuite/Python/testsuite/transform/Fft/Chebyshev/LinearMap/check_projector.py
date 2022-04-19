#! /usr/bin/env python

import numpy as np
from quicc.geometry.chebyshev.chebyshev_linear import linear_map
import quicc.testsuite.utils as utils
import numpy.polynomial.chebyshev as cheby

base_dir = 'data/Transform/Fft/Chebyshev/LinearMap/Projector/'

def grid(size):
    """Get chebyshev grid"""

    x, w = cheby.chebgauss(size)
    return x

def check_d1(lower, upper):
    """Check D^1 projector"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    x = grid(data.shape[0])
    a, b = linear_map(lower, upper)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(c,1,1/a))
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1y1(lower, upper):
    """Check D^1 Y^1 projector"""

    print("\t" + "Checking D1Y1...")
    data = utils.read_complex(base_dir + "d1y1.dat")
    x = grid(data.shape[0])
    a, b = linear_map(lower, upper)
    ref = 1j*np.zeros(data.shape)
    cx = np.array([b, a])
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(cheby.chebmul(cx,c),1,1/a))
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d2(lower, upper):
    """Check D^2 projector"""

    print("\t" + "Checking D2...")
    data = utils.read_complex(base_dir + "d2.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(c,2,1/a))
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d3(lower, upper):
    """Check D^3 projector"""

    print("\t" + "Checking D3...")
    data = utils.read_complex(base_dir + "d3.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(c,3,1/a))
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d4(lower, upper):
    """Check D^4 projector"""

    print("\t" + "Checking D4...")
    data = utils.read_complex(base_dir + "d4.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(c,4,1/a))
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divy1d1y1(lower, upper):
    """Check 1/Y^1 D^1 Y^1 projector"""

    print("\t" + "Checking DivY1D1Y1...")
    data = utils.read_complex(base_dir + "divy1d1y1.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    y = a*x + b
    ref = 1j*np.zeros(data.shape)
    cx = np.array([b, a])
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, cheby.chebder(cheby.chebmul(cx,c),1,1/a))/y
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divy1(lower, upper):
    """Check 1/Y^1 D projector"""

    print("\t" + "Checking DivY1...")
    data = utils.read_complex(base_dir + "divy1.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    y = a*x + b
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, c)/y
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divy2(lower, upper):
    """Check 1/Y^2 projector"""

    print("\t" + "Checking DivY2...")
    data = utils.read_complex(base_dir + "divy2.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    y = a*x + b
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, c)/(y**2)
        ref[:,i] = (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_p(lower, upper):
    """Check P projector"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        t = (np.pi + i)*cheby.chebval(x, c)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_sphradlapl(lower, upper):
    """Check SphRadLapl | zero projector"""

    print("\t" + "Checking SphRadLapl...")
    data = utils.read_complex(base_dir + "sphradlapl.dat")
    a, b = linear_map(lower, upper)
    x = grid(data.shape[0])
    y = a*x + b
    ref = 1j*np.zeros(data.shape)
    cx = np.array([b, a])
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        cc = cheby.chebder(c,1,1/a)
        cc = cheby.chebmul(cx,cc)
        cc = cheby.chebmul(cx,cc)
        cc = cheby.chebder(cc,1,1/a)
        t = (np.pi + i)*cheby.chebval(x, cc)/y**2
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

print("Chebyshev linear map FFT projector")
yi = 7./13.
yo = 20./13.
code = 0
code += check_d1(yi, yo)
code += check_d1y1(yi, yo)
code += check_d2(yi, yo)
code += check_d3(yi, yo)
code += check_d4(yi, yo)
code += check_divy1d1y1(yi, yo)
code += check_divy1(yi, yo)
code += check_divy2(yi, yo)
code += check_p(yi, yo)
code += check_sphradlapl(yi, yo)

utils.test_summary(code)

import sys
sys.exit(code)
