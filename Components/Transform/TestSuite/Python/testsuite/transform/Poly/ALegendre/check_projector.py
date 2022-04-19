#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils
import numpy.polynomial.legendre as leg
import scipy.special as spec

base_dir = 'data/Transform/Poly/ALegendre/Projector/'

def idx2m(i):
    """Map column index to correct l,m"""
    
    d = {
            0:0, 
            1:1, 
            2:2,
            3:3,
            4:7,
            5:8,
            6:8,
            7:9,
        }
    return d[i]

def grid(size):
    """Get chebyshev grid"""

    x, w = leg.leggauss(size)
    return x

def aleg(m, l, x):
    """Normalized associated legendre function"""

    norm = (2.0*l+1)/(4.0*np.pi)
    if l > 0:
        norm *= np.exp(spec.gammaln(l-m+1) - spec.gammaln(l+m+1))
    return spec.lpmv(m, l, x)*np.sqrt(norm)

def daleg(m, l, x):
    """Normalized associated legendre function"""

    norm = (2.0*l+1)/(4.0*np.pi)
    if l > 0:
        norm *= np.exp(spec.gammaln(l-m+1) - spec.gammaln(l+m+1))
    if l > 0:
        d = -0.5*((l+m)*(l-m+1)*spec.lpmv(m-1, l, x) - spec.lpmv(m+1, l, x))
    else:
        d = 0.0
    return d*np.sqrt(norm)

def check_d1(n):
    """Check D^1 projector"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = (np.pi + i)*daleg(idx2m(i), l,x)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divs1(n):
    """Check 1/Sin projector"""

    print("\t" + "Checking DivS1...")
    data = utils.read_complex(base_dir + "divs1.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(1, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = (np.pi + i)*aleg(idx2m(i), l,x)/np.sin(np.arccos(x))
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divs1d1s1(n):
    """Check 1/Sin D^1 Sin projector"""

    print("\t" + "Checking DivS1D1S1...")
    data = utils.read_complex(base_dir + "divs1d1s1.dat")
    ref = -42.*np.ones(data.shape)
    return utils.transform_error(data, ref)

def check_divs1dp(n):
    """Check 1/Sin D_phi projector"""

    print("\t" + "Checking DivS1Dp...")
    data = utils.read_complex(base_dir + "divs1dp.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(1, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        m = idx2m(i)
        t = m*1j*(np.pi + i)*aleg(m, l,x)/np.sin(np.arccos(x))
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_ll(n):
    """Check l(l+1) projector"""

    print("\t" + "Checking Ll...")
    data = utils.read_complex(base_dir + "ll.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = l*(l+1)*(np.pi + i)*aleg(idx2m(i), l,x)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_lld1(n):
    """Check l(l+1) D^1 projector"""

    print("\t" + "Checking LlD1...")
    data = utils.read_complex(base_dir + "lld1.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = l*(l+1)*(np.pi + i)*daleg(idx2m(i), l,x)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_lldivs1(n):
    """Check l(l+1) 1/Sin projector"""

    print("\t" + "Checking LlDivS1...")
    data = utils.read_complex(base_dir + "lldivs1.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = l*(l+1)*(np.pi + i)*aleg(idx2m(i), l,x)/np.sin(np.arccos(x))
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_lldivs1dp(n):
    """Check l(l+1) 1/Sin D_phi projector"""

    print("\t" + "Checking LlDivS1Dp...")
    data = utils.read_complex(base_dir + "lldivs1dp.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(1, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        m = idx2m(i)
        t = l*(l+1)*m*1j*(np.pi + i)*aleg(m, l,x)/np.sin(np.arccos(x))
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_p(n):
    """Check P projector"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    x = grid(data.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        c = np.zeros((i+1,))
        c[-1] = 1
        l = idx2m(i)+i
        t = (np.pi + i)*aleg(idx2m(i), l,x)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

print("Associated Legendre projector")
code = 0
code += check_d1(20)
code += check_divs1(20)
#code += check_divs1d1s1(20)
code += check_divs1dp(20)
code += check_ll(20)
code += check_lld1(20)
code += check_lldivs1(20)
code += check_lldivs1dp(20)
code += check_p(20)

utils.test_summary(code)

import sys
sys.exit(code)
