#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Fft/Fourier/Mixed/Integrator/'

def check_d1():
    """Check D^1 integrator"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[i,i] = (i*1j)*(np.pi + i)*(0.5-0.5j)
    return utils.transform_error(data, ref)

def check_d1_neg():
    """Check D^1 | -P integrator"""

    print("\t" + "Checking D1_Neg...")
    data = utils.read_complex(base_dir + "d1_neg.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[i,i] = (i*1j)*(np.pi + i)*(0.5-0.5j)
    ref[0,0] = -np.pi
    return utils.transform_error(data, ref)

def check_d1_p():
    """Check D^1 | P integrator"""

    print("\t" + "Checking D1_P...")
    data = utils.read_complex(base_dir + "d1_p.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[i,i] = (i*1j)*(np.pi + i)*(0.5-0.5j)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

def check_d2():
    """Check D^2 integrator"""

    print("\t" + "Checking D2...")
    data = utils.read_complex(base_dir + "d2.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(0, ref.shape[1]):
        ref[i,i] = (i*1j)**2*(np.pi + i)*(0.5-0.5j)
    return utils.transform_error(data, ref)

def check_p():
    """Check P integrator"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = (np.pi + i)*(0.5-0.5j)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

print("Fourier mixed FFT integrator")
code = 0
code += check_d1()
code += check_d1_neg()
code += check_d1_p()
code += check_d2()
code += check_p()

utils.test_summary(code)

import sys
sys.exit(code)
