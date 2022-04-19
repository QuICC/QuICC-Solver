#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils

base_dir = 'data/Transform/Poly/ALegendre/Integrator/'

def check_d1(n):
    """Check D^1 integrator"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 2
        cm = (l - 1.0)/(2.0*np.pi)*np.sqrt((l + m)*(l - m)/((2.0*l - 1.0)*(2.0*l+1.0)))
        cp = -(l + 2.0)/(2.0*np.pi)*np.sqrt((l + m + 1.0)*(l - m + 1.0)/((2.0*l + 1.0)*(2.0*l+3.0)))
        ref[1,i] = cm*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        ref[3,i] = cp*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divlld1(n):
    """Check 1/l(l+1) D^1 integrator"""

    print("\t" + "Checking DivLlD1...")
    data = utils.read_complex(base_dir + "divlld1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 2
        cm = (1.0/((l-1.0)*l))*(l - 1.0)/(2.0*np.pi)*np.sqrt((l + m)*(l - m)/((2.0*l - 1.0)*(2.0*l+1.0)))
        cp = -(1.0/((l+1.0)*(l+2.0)))*(l + 2.0)/(2.0*np.pi)*np.sqrt((l + m + 1.0)*(l - m + 1.0)/((2.0*l + 1.0)*(2.0*l+3.0)))
        ref[1,i] = cm*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        ref[3,i] = cp*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divlldivs1dp(n):
    """Check 1/l(l+1) 1/Sin D_phi integrator"""

    print("\t" + "Checking DivLlDivS1Dp...")
    data = utils.read_complex(base_dir + "divlldivs1dp.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 1
        c = 1.0/(l*(l+1))
        ref[1,i] = -1j*m*c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divlldivs1(n):
    """Check 1/l(l+1) 1/Sin integrator"""

    print("\t" + "Checking DivLlDivS1...")
    data = utils.read_complex(base_dir + "divlldivs1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(1,data.shape[1]):
        m = ms[i]
        l = m + 1
        c = 1.0/(l*(l+1))
        ref[1,i] = c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divll(n):
    """Check 1/l(l+1) integrator"""

    print("\t" + "Checking DivLl...")
    data = utils.read_complex(base_dir + "divll.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 1
        if l == 0:
            c = 0.0
        else:
            c = 1.0/(l*(l+1))
        ref[1,i] = c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divs1dp(n):
    """Check 1/Sin D_phi integrator"""

    print("\t" + "Checking DivS1Dp...")
    data = utils.read_complex(base_dir + "divs1dp.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        ref[1,i] = -1j*m*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_divs1(n):
    """Check 1/Sin integrator"""

    print("\t" + "Checking DivS1...")
    data = utils.read_complex(base_dir + "divs1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(1,data.shape[1]):
        m = ms[i]
        ref[1,i] = (np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

def check_ll2(n):
    """Check (l(l+1))^2 integrator"""

    print("\t" + "Checking Ll2...")
    data = utils.read_complex(base_dir + "ll2.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 1
        c = (l*(l+1))**2
        ref[1,i] = c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref, True)

def check_lld1(n):
    """Check l(l+1) D^1 integrator"""

    print("\t" + "Checking LlD1...")
    data = utils.read_complex(base_dir + "lld1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 2
        cm = ((l-1.0)*l)*(l - 1.0)/(2.0*np.pi)*np.sqrt((l + m)*(l - m)/((2.0*l - 1.0)*(2.0*l+1.0)))
        cp = -((l+1.0)*(l+2.0))*(l + 2.0)/(2.0*np.pi)*np.sqrt((l + m + 1.0)*(l - m + 1.0)/((2.0*l + 1.0)*(2.0*l+3.0)))
        ref[1,i] = cm*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        ref[3,i] = cp*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref, True)

def check_lldivs1dp(n):
    """Check l(l+1) 1/Sin D_phi integrator"""

    print("\t" + "Checking LlDivS1Dp...")
    data = utils.read_complex(base_dir + "lldivs1dp.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m + 1
        c = (l*(l+1))
        ref[1,i] = -1j*m*c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref, True)

def check_lldivs1(n):
    """Check l(l+1) 1/Sin integrator"""
   
    print("\t" + "Checking l(l+1) 1/Sin...")
    data = utils.read_complex(base_dir + "lldivs1.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(1,data.shape[1]):
        m = ms[i]
        l = m + 1
        c = (l*(l+1))
        ref[1,i] = c*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref, True)

def check_ll(n):
    """Check l(l+1) integrator"""

    print("\t" + "Checking Ll...")
    data = utils.read_complex(base_dir + "ll.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        l = m+1
        ref[1,i] = l*(l+1)*(np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref, True)

def check_p(n):
    """Check P integrator"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    ref = 1j*np.zeros(data.shape)
    ms = [0, 1, 2, 3, 7, 8 ,8 ,9]
    for i in range(0,data.shape[1]):
        m = ms[i]
        ref[2,i] = (np.pi + i)*(m+1.0)*(1.0 - (i > 0)*3.0j)/(2.0*np.pi)
        if i > 0:
            ref[-m::,i] = -42424242.4242
    return utils.transform_error(data, ref)

print("Associated Legendre integrator")
code = 0
code += check_d1(20)
code += check_divlld1(20)
code += check_divlldivs1dp(20)
code += check_divlldivs1(20)
code += check_divll(20)
code += check_divs1dp(20)
code += check_divs1(20)
code += check_ll2(20)
code += check_lld1(20)
code += check_lldivs1dp(20)
code += check_lldivs1(20)
code += check_ll(20)
code += check_p(20)

utils.test_summary(code)

import sys
sys.exit(code)
