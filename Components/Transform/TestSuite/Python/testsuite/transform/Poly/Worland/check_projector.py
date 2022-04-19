#! /usr/bin/env python

import numpy as np
import quicc.testsuite.utils as utils
import quicc.geometry.worland.wnl as wnl

base_dir = 'data/Transform/Poly/Worland/Projector/'

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

def worland(n, l, nr):
    """Normalized Worland polynomial"""

    return wnl.eval_poly(n, l, nr)[::-1]

def dworland(n, l, k, nr):
    """Normalized Worland polynomial"""

    return wnl.worland_djacobi(n, l, k, nr)[::-1]

def check_cyllaplh(n):
    """Check CylLaplh projector"""

    print("\t" + "Checking CylLaplh...")
    data = utils.read_complex(base_dir + "cyllaplh.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
        if i == 0:
            t = np.zeros(rg.shape)
        elif i == 1:
            t = (np.pi + i)*4.0*(i + l)*((l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg)
        else:
            c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
            t = (np.pi + i)*4.0*(i + l)*((l + i + 1.0)*c2*dworland(i-2, l, 2, nr) + (l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_cyllaplh_divr1d1r1(n):
    """Check CylLaplh | 1/R^1 D^1 R^1 projector"""

    print("\t" + "Checking CylLaplh_DivR1D1R1...")
    data = utils.read_complex(base_dir + "cyllaplh_divr1d1r1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if l == 0:
            l = 1
            if i == 0:
                t = (np.pi + i)*(l+1.0)*worland(i, l, nr)/rg
            else:
                c = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
                t = (np.pi + i)*(2.0*(i + l)*c*rg*dworland(i-1, l, 1, nr) + (l + 1.0)*worland(i, l, nr))/rg
        else:
            c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
            if i == 0:
                t = np.zeros(rg.shape)
            elif i == 1:
                t = (np.pi + i)*4.0*(i + l)*((l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg)
            else:
                c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
                t = (np.pi + i)*4.0*(i + l)*((l + i + 1.0)*c2*dworland(i-2, l, 2, nr) + (l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1(n):
    """Check D^1 projector"""

    print("\t" + "Checking D1...")
    data = utils.read_complex(base_dir + "d1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        c = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
        if i == 0:
            t = (np.pi + i)*(l*worland(i, l, nr)/rg)
        else:
            t = (np.pi + i)*(2.0*(i + l)*c*dworland(i-1, l, 1, nr) + l*worland(i, l, nr)/rg)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1cyllaplh(n):
    """Check D^1 CylLaplh projector"""

    print("\t" + "Checking D1CylLaplh...")
    data = utils.read_complex(base_dir + "d1cyllaplh.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
        if i == 0:
            t = np.zeros(rg.shape)
        elif i == 1:
            t = (np.pi + i)*4.0*(i + l)*(l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
        elif i == 2:
            c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
            t = (np.pi + i)*4.0*(i + l)*((3.0*l + 4.0)*(l + i + 1.0)*c2*dworland(i-2, l, 2, nr)/rg + l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
        else:
            c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
            c3 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-3,l)
            t = (np.pi + i)*4.0*(i + l)*(2.0*(l + i + 1.0)*(l + i + 2.0)*c3*dworland(i-3, l, 3, nr) + (3.0*l + 4.0)*(l + i + 1.0)*c2*dworland(i-2, l, 2, nr)/rg + l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1cyllaplh_d1divr1d1r1(n):
    """Check D^1 CylLaplh | D^1 1/R^1 D^1 R^1 projector"""

    print("\t" + "Checking D1CylLaplh_D1DivR1D1R1...")
    data = utils.read_complex(base_dir + "d1cyllaplh_d1divr1d1r1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if l == 0:
            l = 1
            if i == 0:
                t = (np.pi + i)*((l + 1.0)*(l - 1.0)*worland(i, l, nr)/rg**2)
            elif i == 1:
                c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
                t = (np.pi + i)*(4.0*(l + 1.0)*(i + l)*c1*rg*dworland(i-1, l, 1, nr)/rg + (l + 1.0)*(l - 1.0)*worland(i, l, nr)/rg**2)
            else:
                c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
                c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
                t = (np.pi + i)*(4.0*(i + l)*(i + l + 1.0)*c2*rg*dworland(i-2, l, 2, nr) + 4.0*(l + 1.0)*(i + l)*c1*rg*dworland(i-1, l, 1, nr)/rg + (l + 1.0)*(l - 1.0)*worland(i, l, nr)/rg**2)
        else:
            c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
            if i == 0:
                t = np.zeros(rg.shape)
            elif i == 1:
                t = (np.pi + i)*4.0*(i + l)*(l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
            elif i == 2:
                c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
                t = (np.pi + i)*4.0*(i + l)*((3.0*l + 4.0)*(l + i + 1.0)*c2*dworland(i-2, l, 2, nr)/rg + l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
            else:
                c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
                c3 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-3,l)
                t = (np.pi + i)*4.0*(i + l)*(2.0*(l + i + 1.0)*(l + i + 2.0)*c3*dworland(i-3, l, 3, nr) + (3.0*l + 4.0)*(l + i + 1.0)*c2*dworland(i-2, l, 2, nr)/rg + l*(l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1r1(n):
    """Check D^1 R^1 projector"""

    print("\t" + "Checking D1R1...")
    data = utils.read_complex(base_dir + "d1r1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if i == 0:
            t = (np.pi + i)*((l+1.0)*worland(i, l, nr))
        else:
            c = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
            t = (np.pi + i)*(2.0*(i + l)*c*rg*dworland(i-1, l, 1, nr) + (l + 1.0)*worland(i, l, nr))
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_d1_p(n):
    """Check D^1 | P projector"""

    print("\t" + "Checking D1_P...")
    data = utils.read_complex(base_dir + "d1_p.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if l == 0:
            t = (np.pi + i)*worland(i, 1, nr)
        else:
            if i == 0:
                t = (np.pi + i)*(l*worland(i, l, nr)/rg)
            else:
                c = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
                t = (np.pi + i)*(2.0*(i + l)*c*dworland(i-1, l, 1, nr) + l*worland(i, l, nr)/rg)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divr1(n):
    """Check 1/R^1 projector"""

    print("\t" + "Checking DivR1...")
    data = utils.read_complex(base_dir + "divr1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        t = (np.pi + i)*worland(i, l, nr)/rg
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divr1cyllaplh_zero(n):
    """Check 1/R^1 CylLaplh | Zero projector"""

    print("\t" + "Checking DivR1CylLaplh_Zero...")
    data = utils.read_complex(base_dir + "divr1cyllaplh_zero.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if l == 0:
            t = np.zeros(rg.shape)
        else:
            c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
            if i == 0:
                t = np.zeros(rg.shape)
            elif i == 1:
                t = (np.pi + i)*4.0*(i + l)*((l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
            else:
                c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
                t = (np.pi + i)*4.0*(i + l)*((l + i + 1.0)*c2*dworland(i-2, l, 2, nr)/rg + (l + 1.0)*c1*dworland(i-1, l, 1, nr)/rg**2)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divr1d1r1(n):
    """Check 1/R^1 D^1 R^1 projector"""

    print("\t" + "Checking DivR1D1R1...")
    data = utils.read_complex(base_dir + "divr1d1r1.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        if i == 0:
            t = (np.pi + i)*(l+1.0)*worland(i, l, nr)/rg
        else:
            c = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
            t = (np.pi + i)*(2.0*(i + l)*c*rg*dworland(i-1, l, 1, nr) + (l + 1.0)*worland(i, l, nr))/rg
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_divr1_zero(n):
    """Check 1/R^1 | Zero projector"""

    print("\t" + "Checking DivR1_Zero...")
    data = utils.read_complex(base_dir + "divr1_zero.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(1, data.shape[1]):
        l = idx2l(i)
        t = (np.pi + i)*worland(i, l, nr)/rg
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_p(n):
    """Check P projector"""

    print("\t" + "Checking P...")
    data = utils.read_complex(base_dir + "p.dat")
    nr = data.shape[0]/2
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        t = (np.pi + i)*worland(i, l, nr)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

def check_sphlapl(n):
    """Check SphLapl projector"""

    print("\t" + "Checking SphLapl...")
    data = utils.read_complex(base_dir + "sphlapl.dat")
    nr = data.shape[0]/2
    rg = grid(nr)
    ref = 1j*np.zeros(data.shape)
    for i in range(0, data.shape[1]):
        l = idx2l(i)
        c1 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-1,l)
        if i == 0:
            t = np.zeros(rg.shape)
        elif i == 1:
            t = (np.pi + i)*2.0*(i + l)*((2.0*l + 3.0)*c1*dworland(i-1, l, 1, nr)/rg)
        else:
            c2 = wnl.get_invnorm(i,l)/wnl.get_invnorm(i-2,l)
            t = (np.pi + i)*2.0*(i + l)*(2.0*(l + i + 1.0)*c2*dworland(i-2, l, 2, nr) + (2.0*l + 3.0)*c1*dworland(i-1, l, 1, nr)/rg)
        ref[:,i] =  (1.0 - (i > 0)*3.0j)*t
    return utils.transform_error(data, ref)

print("Worland projector")
code = 0
code += check_cyllaplh(20)
code += check_cyllaplh_divr1d1r1(20)
code += check_d1(20)
code += check_d1cyllaplh(20)
code += check_d1cyllaplh_d1divr1d1r1(20)
code += check_d1r1(20)
code += check_d1_p(20)
code += check_divr1(20)
code += check_divr1cyllaplh_zero(20)
code += check_divr1d1r1(20)
code += check_divr1_zero(20)
code += check_p(20)
code += check_sphlapl(20)

utils.test_summary(code)

import sys
sys.exit(code)
