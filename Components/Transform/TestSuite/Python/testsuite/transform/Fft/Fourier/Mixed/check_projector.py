#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Fft/Fourier/Mixed/Projector/'

def check_d1():
    """Check D^1 integrator"""

    print("\t" + "Checking D1...")
    data = utils.read_real(base_dir + "d1.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[:,i] = i*(np.pi + i)*(-np.sin(i*phi) + np.cos(i*phi))
    return utils.transform_error(data, ref)

def check_d2():
    """Check D^2 integrator"""

    print("\t" + "Checking D2...")
    data = utils.read_real(base_dir + "d2.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[:,i] = -i**2*(np.pi + i)*(np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_d3():
    """Check D^3 integrator"""

    print("\t" + "Checking D3...")
    data = utils.read_real(base_dir + "d3.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[:,i] = i**3*(np.pi + i)*(np.sin(i*phi) - np.cos(i*phi))
    return utils.transform_error(data, ref)

def check_d4():
    """Check D^4 integrator"""

    print("\t" + "Checking D4...")
    data = utils.read_real(base_dir + "d4.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[:,i] = i**4*(np.pi + i)*(np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_p():
    """Check P integrator"""

    print("\t" + "Checking P...")
    data = utils.read_real(base_dir + "p.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[:,i] = (np.pi + i)*(np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

print("Fourier mixed FFT projector")
code = 0
code += check_d1()
code += check_d2()
code += check_d3()
code += check_d4()
code += check_p()

utils.test_summary(code)

import sys
sys.exit(code)
