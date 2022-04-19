#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Fft/Fourier/Complex/Integrator/'

def check_d1():
    """Check D^1 integrator"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = i*1j*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -i*1j*(np.pi + i)*(0.5+0.5j)
    return utils.transform_error(data, ref)

def check_d1_neg():
    """Check D^1 | -P integrator"""

    print("\t" + "Checking D1_Neg...")
    data = utils.read_complex(base_dir + "d1_neg.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = i*1j*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -i*1j*(np.pi + i)*(0.5+0.5j)
    ref[0,0] = -np.pi

    return utils.transform_error(data, ref)

def check_d1_p():
    """Check D^1 | P integrator"""

    print("\t" + "Checking D1_P...")
    data = utils.read_complex(base_dir + "d1_p.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = i*1j*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -i*1j*(np.pi + i)*(0.5+0.5j)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

def check_d2():
    """Check D^2 integrator"""

    print("\t" + "Checking D2...")
    data = utils.read_complex(base_dir + "d2.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = -i**2*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -i**2*(np.pi + i)*(0.5+0.5j)
    return utils.transform_error(data, ref)

def check_df1invlapl2d():
    """Check D(fast) Inv Lapl2D integrator"""

    print("\t" + "Checking DfInvLapl2D...")
    data = utils.read_complex(base_dir + "df1invlapl2d.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = -(i*1j/(i**2 + i**2))*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = (i*1j/(i**2 + i**2))*(np.pi + i)*(0.5+0.5j)
    return utils.transform_error(data, ref)

def check_invlapl2d():
    """Check Inv Lapl2D integrator"""

    print("\t" + "Checking InvLapl2D...")
    data = utils.read_complex(base_dir + "invlapl2d.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = -(1.0/(i**2 + i**2))*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -(1.0/(i**2 + i**2))*(np.pi + i)*(0.5+0.5j)
    return utils.transform_error(data, ref)

def check_lapl2d():
    """Check Lapl2D D integrator"""

    print("\t" + "Checking Lapl2D...")
    data = utils.read_complex(base_dir + "lapl2d.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = -(i**2 + i**2)*(np.pi + i)*(0.5-0.5j)
        ref[-i,i] = -(i**2 + i**2)*(np.pi + i)*(0.5+0.5j)
    return utils.transform_error(data, ref)

def check_mean():
    """Check mean integrator"""

    print("\t" + "Checking Mean...")
    data = utils.read_complex(base_dir + "mean.dat")
    ref = 0j*np.zeros(data.shape)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

def check_p_clean():
    """Check P clean integrator"""

    print("\t" + "Checking P_Clean...")
    data = utils.read_complex(base_dir + "p_clean.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = (np.pi + i)*(0.5-0.5j)
        ref[-i,i] = (np.pi + i)*(0.5+0.5j)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

def check_p():
    """Check P integrator"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[i,i] = (np.pi + i)*(0.5-0.5j)
        ref[-i,i] = (np.pi + i)*(0.5+0.5j)
    ref[0,0] = np.pi
    return utils.transform_error(data, ref)

print("Fourier complex FFT integrator")
code = 0
code += check_d1()
code += check_d1_neg()
code += check_d1_p()
code += check_d2()
code += check_df1invlapl2d()
code += check_invlapl2d()
code += check_lapl2d()
code += check_mean()
code += check_p_clean()
code += check_p()

utils.test_summary(code)

import sys
sys.exit(code)
