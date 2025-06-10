#! /usr/bin/env python

import numpy as np
import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.spherical.shell_radius as rad
from quicc.geometry.chebyshev.chebyshev_linear import linear_map
import quicc.testsuite.utils as utils
import numpy.polynomial.chebyshev as cheby

base_dir = 'data/Transform/Fft/Chebyshev/LinearMap/Integrator/'

def check_i2d1(lower, upper):
    """Check I^2 Y^1 integrator"""

    print("\t" + "Checking I2D1...")
    data = utils.read_complex(base_dir + "i2d1.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        ref[:,i] = c1d.i2d1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2d1_i2(lower, upper):
    """Check I^2 Y^1 | I^2 integrator"""

    print("\t" + "Checking I2D1_I2...")
    data = utils.read_complex(base_dir + "i2d1_i2.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = c1d.i2(data.shape[0], lower, upper, {0:0})*ref[:,i]
        else:
            ref[:,i] = c1d.i2d1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2(lower, upper):
    """Check I^2 integrator"""

    print("\t" + "Checking I2...")
    data = utils.read_complex(base_dir + "i2.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        ref[:,i] = c1d.i2(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2y1d1y1_zero(lower, upper):
    """Check I^2 Y D Y integrator"""

    print("\t" + "Checking I2Y1D1Y1...")
    data = utils.read_complex(base_dir + "i2y1d1y1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i2r1d1r1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2y2d1y1_zero(lower, upper):
    """Check I^2 Y^2 D Y integrator"""

    print("\t" + "Checking I2Y2D1Y1...")
    data = utils.read_complex(base_dir + "i2y2d1y1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i2r2d1r1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2y1_zero(lower, upper):
    """Check I^2 Y^1  | Zero integrator"""

    print("\t" + "Checking I2Y1_Zero...")
    data = utils.read_complex(base_dir + "i2y1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i2r1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2y2_zero(lower, upper):
    """Check I^2 Y^2 | Zero integrator"""

    print("\t" + "Checking I2Y2_Zero...")
    data = utils.read_complex(base_dir + "i2y2_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i2r2(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4d1(lower, upper):
    """Check I^4 D integrator"""

    print("\t" + "Checking I4D1...")
    data = utils.read_complex(base_dir + "i4d1.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        ref[:,i] = c1d.i4d1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4d1_i2(lower, upper):
    """Check I^4 D^1 | I^2 integrator"""

    print("\t" + "Checking I4D1_I2...")
    data = utils.read_complex(base_dir + "i4d1_i2.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = c1d.i2(data.shape[0], lower, upper, {0:0})*ref[:,i]
        else:
            ref[:,i] = c1d.i4d1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4(lower, upper):
    """Check I^4 integrator"""

    print("\t" + "Checking I4...")
    data = utils.read_complex(base_dir + "i4.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        ref[:,i] = c1d.i4(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4y3d1y1_zero(lower, upper):
    """Check I^4 Y^3 D Y | zero integrator"""

    print("\t" + "Checking I4Y3D1Y1_Zero...")
    data = utils.read_complex(base_dir + "i4y3d1y1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i4r3d1r1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4y3_zero(lower, upper):
    """Check I^4 Y^3 | zero integrator"""

    print("\t" + "Checking I4Y3_Zero...")
    data = utils.read_complex(base_dir + "i4y3_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        if i == 0:
            ref[:,i] = 0
        else:
            ref[:,i] = rad.i4r3(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_p(lower, upper):
    """Check P integrator"""
   
    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_y1(lower, upper):
    """Check Y^1 integrator"""

    print("\t" + "Checking Y1...")
    data = utils.read_complex(base_dir + "y1.dat")
    ref = 1j*np.zeros(data.shape)
    a, b = linear_map(lower, upper)
    cx = np.array([b, a])
    for i in range(0,data.shape[1]):
        t = 0.5*(1 + (i == 0))
        ref[i,i] = (np.pi + i)*(1.0 - (i > 0)*3.0j)*t
        ref[:,i] = rad.r1(data.shape[0], lower, upper, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

print("Chebyshev linear map FFT integrator")
yi = 0.0
yo = 1.0
code = 0
code += check_i2d1(yi, yo)
code += check_i2d1_i2(yi, yo)
code += check_i2(yi, yo)
code += check_i2y1d1y1_zero(yi, yo)
code += check_i2y1_zero(yi, yo)
code += check_i2y2_zero(yi, yo)
code += check_i4d1(yi, yo)
code += check_i4d1_i2(yi, yo)
code += check_i4(yi, yo)
code += check_i4y3d1y1_zero(yi, yo)
code += check_i4y3_zero(yi, yo)
code += check_p(yi, yo)
code += check_y1(yi, yo)

utils.test_summary(code)

import sys
sys.exit(code)
