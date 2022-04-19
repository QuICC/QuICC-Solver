#! /usr/bin/env python

import numpy as np
from quicc.geometry.chebyshev.chebyshev_linear import linear_map
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Fft/Chebyshev/LinearMap/Reductor/'

def check_energyd1y1(n):
    """Check Energy D^1 Y^1 reductor"""

    print("\t" + "Checking EnergyD1Y1...")
    data = utils.read_real(base_dir + "energyd1y1.dat")
    # Hard coded energy for T_18 + T_19
    nrg = 14062.57655968394
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg
    return utils.transform_error(data, ref)

def check_energyy2(n):
    """Check Energy Y^2 reductor"""

    print("\t" + "Checking EnergyY2...")
    data = utils.read_real(base_dir + "energyy2.dat")
    # Hard coded energy for T_18 + T_19
    nrg = 1.506154767924143
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg
    return utils.transform_error(data, ref)

def check_energy(n):
    """Check Energy reductor"""

    print("\t" + "Checking Energy...")
    data = utils.read_real(base_dir + "energy.dat")
    # Hard coded energy for T_18 + T_19
    nrg = 0.9992673992673993
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg
    return utils.transform_error(data, ref)

print("Chebyshev linear map FFT reductor")
code = 0
code += check_energyd1y1(20)
code += check_energyy2(20)
code += check_energy(20)

utils.test_summary(code)

import sys
sys.exit(code)
