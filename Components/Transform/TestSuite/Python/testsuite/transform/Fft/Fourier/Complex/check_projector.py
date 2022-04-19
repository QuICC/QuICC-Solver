#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Fft/Fourier/Complex/Projector/'

def check_d1():
    """Check D^1 projector"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = i*(np.pi + i)*(1.0 - 1j)*(-3.0*np.sin(i*phi) + np.cos(i*phi))
    return utils.transform_error(data, ref)

def check_d2():
    """Check D^2 projector"""

    print("\t" + "Checking D2...")
    data = utils.read_complex(base_dir + "d2.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = -i**2*(np.pi + i)*(1.0 - 1j)*(3.0*np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_d3():
    """Check D^3 projector"""

    print("\t" + "Checking D3...")
    data = utils.read_complex(base_dir + "d3.dat")
    ref = -42.*np.ones(data.shape)
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = -i**3*(np.pi + i)*(1.0 - 1j)*(-3.0*np.sin(i*phi) + np.cos(i*phi))
    return utils.transform_error(data, ref)

def check_d4():
    """Check D^4 projector"""

    print("\t" + "Checking D4...")
    data = utils.read_complex(base_dir + "d4.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = i**4*(np.pi + i)*(1.0 - 1j)*(3.0*np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_df1lapl2d():
    """Check D(fast) Lapl2D projector"""

    print("\t" + "Checking DfLapl2D...")
    data = utils.read_complex(base_dir + "df1lapl2d.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = -i*(i**2 + i**2)*(np.pi + i)*(1.0 - 1j)*(-3.0*np.sin(i*phi) + np.cos(i*phi))
    return utils.transform_error(data, ref)

def check_ds1lapl2d():
    """Check D(slow) Lapl2D projector"""

    print("\t" + "Checking DsLapl2D...")
    data = utils.read_complex(base_dir + "ds1lapl2d.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = -i*1j*(i**2 + i**2)*(np.pi + i)*(1.0 - 1j)*(3.0*np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_lapl2d():
    """Check Lapl2D projector"""

    print("\t" + "Checking Lapl2D...")
    data = utils.read_complex(base_dir + "lapl2d.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = -(i**2 + i**2)*(np.pi + i)*(1.0 - 1j)*(3.0*np.cos(i*phi) + np.sin(i*phi))
    return utils.transform_error(data, ref)

def check_mean():
    """Check mean projector"""

    print("\t" + "Checking Mean...")
    data = utils.read_complex(base_dir + "mean.dat")
    ref = 0j*np.zeros(data.shape)
    ref[:,0] = np.pi
    return utils.transform_error(data, ref)

def check_p():
    """Check P projector"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    n = data.shape[0]
    phi = np.linspace(0, n-1, n)*2.0*np.pi/n
    ref = 0j*np.zeros(data.shape)
    for i in range(1, ref.shape[1]):
        ref[:,i] = (np.pi + i)*(1.0 - 1j)*(3.0*np.cos(i*phi) + np.sin(i*phi))
    ref[:,0] = np.pi
    return utils.transform_error(data, ref)

print("Chebyshev linear map FFT projector")
code = 0
code += check_d1()
code += check_d2()
code += check_d3()
code += check_d4()
code += check_df1lapl2d()
code += check_ds1lapl2d()
code += check_lapl2d()
code += check_mean()
code += check_p()

utils.test_summary(code)

import sys
sys.exit(code)
