#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils
import quicc.geometry.spherical.sphere_radius_worland as sphere
import quicc.geometry.cylindrical.cylinder_radius_worland as cylinder
import quicc.geometry.worland.wnl as wnl

base_dir = 'data/Transform/Poly/Worland/Integrator/'

def idx2l(i):
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
    """Get Worland-Chebyshev grid"""
    
    return wnl.get_grid(size)[::-1]

def weights(size):
    """Get Worland-Chebyshev grid"""
    
    return wnl.get_weights(size)[::-1]

def worland(n, l, nr):
    """Normalized Worland polynomial"""

    return wnl.eval_poly(n, l, nr)[::-1]

def check_i2divr1d1r1_zero(n):
    """Check I^2 1/R^1 D^1 R^1 | Zero integrator"""

    print("\t" + "Checking I2DivR1D1R1_Zero...")
    data = utils.read_complex(base_dir + "i2divr1d1r1_zero.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    w = weights(rg.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            f = rg**l*(1.0 + rg**2 + rg**4 + rg**6)
            ref[0:4,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)*np.array([np.sum(f*w*worland(j, l, nr)) for j in range(0,4)])
            ref[:,i] = sphere.i2(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2divr1_zero(n):
    """Check I^2 1/R^1 | Zero integrator"""

    print("\t" + "Checking I2DivR1_Zero...")
    data = utils.read_complex(base_dir + "i2divr1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = sphere.i2(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i2_zero(n):
    """Check I^2 | Zero integrator"""

    print("\t" + "Checking I2_Zero...")
    data = utils.read_complex(base_dir + "i2_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = sphere.i2(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4divr1d1r1_i2(n):
    """Check I^4 1/R^1 D^1 R^1 | I^2 integrator"""

    print("\t" + "Checking I4DivR1D1R1_I2...")
    data = utils.read_complex(base_dir + "i4divr1d1r1_i2.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    w = weights(rg.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            l = 1
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - 3.0j)
            ref[:,i] = sphere.i2(data.shape[0], l, {0:0})*ref[:,i]
        else:
            f = rg**l*(1.0 + rg**2 + rg**4 + rg**6)
            ref[0:4,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)*np.array([np.sum(f*w*worland(j, l, nr)) for j in range(0,4)])
            ref[:,i] = sphere.i4(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4divr1d1r1_zero(n):
    """Check I^4 1/R^1 D^1 R^1 | Zero integrator"""

    print("\t" + "Checking I4DivR1D1R1_Zero...")
    data = utils.read_complex(base_dir + "i4divr1d1r1_zero.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    w = weights(rg.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            f = rg**l*(1.0 + rg**2 + rg**4 + rg**6)
            ref[0:4,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)*np.array([np.sum(f*w*worland(j, l, nr)) for j in range(0,4)])
            ref[:,i] = sphere.i4(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i4divr1_zero(n):
    """Check I^4 1/R^1 | Zero integrator"""

    print("\t" + "Checking I4DivR1_Zero...")
    data = utils.read_complex(base_dir + "i4divr1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = sphere.i4(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i6cyllaplh_i4d1r1(n):
    """Check I^6 CylLalph | I^4 D^1 R^1 integrator"""

    print("\t" + "Checking I6CylLaplh_I4D1R1...")
    data = utils.read_complex(base_dir + "i6cyllaplh_i4d1r1.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = cylinder.i4dr(data.shape[0], 1, {0:0})*ref[:,i]
        else:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = cylinder.i6laplh(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i6divr1d1r1_i4(n):
    """Check I^6 1/R^1 D^1 R1 | I^4 integrator"""

    print("\t" + "Checking I6DivR1D1R1_I4...")
    data = utils.read_complex(base_dir + "i6divr1d1r1_i4.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    w = weights(rg.shape[0])
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            l = 1
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - 3.0j)
            ref[:,i] = sphere.i4(data.shape[0], l, {0:0})*ref[:,i]
        else:
            f = rg**l*(1.0 + rg**2 + rg**4 + rg**6)
            ref[0:4,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)*np.array([np.sum(f*w*worland(j, l, nr)) for j in range(0,4)])
            ref[:,i] = cylinder.i6(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_i6divr1_zero(n):
    """Check I^6 1/R^1 | Zero  integrator"""

    print("\t" + "Checking I6DivR1_Zero...")
    data = utils.read_complex(base_dir + "i6divr1_zero.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        if l == 0:
            ref[:,i] = 0
        else:
            ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
            ref[:,i] = cylinder.i6(data.shape[0], l, {0:0})*ref[:,i]
    return utils.transform_error(data, ref)

def check_p(n):
    """Check P integrator"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
    return utils.transform_error(data, ref)

def check_r1(n):
    """Check R^1 integrator"""
   
    print("\t" + "Checking R1...")
    data = utils.read_complex(base_dir + "r1.dat")
    ref = 1j*np.zeros(data.shape)
    for i in range(0,data.shape[1]):
        l = idx2l(i)
        ref[i+3,i] = (np.pi + i)*(l+1.0)*(1.0 - (i > 0)*3.0j)
    return utils.transform_error(data, ref)

print("Worland integrator")
code = 0
code += check_i2divr1d1r1_zero(20)
code += check_i2divr1_zero(20)
code += check_i2_zero(20)
code += check_i4divr1d1r1_i2(20)
code += check_i4divr1d1r1_zero(20)
code += check_i4divr1_zero(20)
code += check_i6cyllaplh_i4d1r1(20)
code += check_i6divr1d1r1_i4(20)
code += check_i6divr1_zero(20)
code += check_p(20)
code += check_r1(20)

utils.test_summary(code)

import sys
sys.exit(code)
