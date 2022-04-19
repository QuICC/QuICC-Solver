#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Poly/Worland/Reductor/'

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

def check_energyd1r1(n):
    """Check Energy D^1 R^1 reductor"""

    print("\t" + "Checking EnergyD1R1...")
    data = utils.read_real(base_dir + "energyd1r1.dat")
    # Hard coded energy for W_18^l + W_19^l
    nrg = {
            0:14792.327713621122,
            1:15700.379351851454,
            2:16621.772926829326,
            3:17556.427346355278,
            4:18504.25119372133,
            5:19465.145512212654,
            6:20439.005949396676,
            7:21425.724417468493,
            8:22425.19038433554,
            9:23437.291880585213,
            10:24461.916286215142,
            }
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        l = idx2l(i)
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg[l]
    return utils.transform_error(data, ref)

def check_energyr2(n):
    """Check Energy R^2 reductor"""

    print("\t" + "Checking EnergyR2...")
    data = utils.read_real(base_dir + "energyr2.dat")
    # Hard coded energy for W_18^l + W_19^l
    nrg = {
            0:0.5088296847735425,
            1:0.5088539500103614,
            2:0.5089884758439588,
            3:0.5092145364363704,
            4:0.5095146220953009,
            5:0.5098727383977952,
            6:0.5102745228725756,
            7:0.5107072485026539,
            8:0.5111597594244834,
            9:0.5116223684464954,
            10:0.5120867355716504,
            }
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        l = idx2l(i)
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg[l]
    return utils.transform_error(data, ref)

def check_energy(n):
    """Check Energy reductor"""

    print("\t" + "Checking Energy...")
    data = utils.read_real(base_dir + "energy.dat")
    # Hard coded energy for W_18^l + W_19^l
    nrg = {
            0:0.8483607423649700,
            1:0.8483849556325944,
            2:0.8478462192995043,
            3:0.8468272436939758,
            4:0.8453980195166429,
            5:0.8436180475230560,
            6:0.8415381214873829,
            7:0.8392017656663448,
            8:0.8366464023707334,
            9:0.8339043067894033,
            10:0.8310033927202925,
            }
    ref = np.zeros(data.shape)
    for i in range(0, data.shape[0]):
        l = idx2l(i)
        ref[i] = (np.pi + i)**2*np.abs(1.0 - (i > 0)*3.0j)**2*nrg[l]
    return utils.transform_error(data, ref)

print("Worland reductor")
code = 0
code += check_energyd1r1(20)
code += check_energyr2(20)
code += check_energy(20)

utils.test_summary(code)

import sys
sys.exit(code)
