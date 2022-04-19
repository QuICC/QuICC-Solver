#!/usr/bin/env python

from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function

import sys, getopt

import h5py
import numpy as np
import re

def getOrdering(scheme):
    """Get data ordering for scheme"""

    if scheme in [b'TTT']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'TFT']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'TFF']:
        idSlow = 'z'
        idMid = 'x'
        idFast = 'y'
    elif scheme in [b'FFF']:
        idSlow = 'x'
        idMid = 'y'
        idFast = 'z'
    elif scheme in [b'AFT', b'WFT']:
        idSlow = 'r'
        idMid = 'theta'
        idFast = 'z'
    elif scheme in [b'SLFm', b'SLFl', b'WLFl', b'WLFm']:
        idSlow = 'r'
        idMid = 'theta'
        idFast = 'phi'
    elif scheme in [b'TT']:
        idSlow = None 
        idMid = 'x'
        idFast = 'z'
    elif scheme in [b'TF']:
        idSlow = None 
        idMid = 'z'
        idFast = 'x'

    return (idFast, idMid, idSlow)

def getGrid(scheme, h5_file):
    """Get the grid for scheme"""

    idFast, idMid, idSlow = getOrdering(scheme)

    gFast = h5_file['mesh']['grid_'+idFast]
    gMid = h5_file['mesh']['grid_'+idMid]
    if idSlow is not None:
        gSlow = h5_file['mesh']['grid_'+idSlow]
    else:
        gSlow = None

    return (gFast, gMid, gSlow)

def getGridSize(scheme, h5_file):
    """Get the grid for scheme"""

    idFast, idMid, idSlow = getOrdering(scheme)

    nFast = h5_file['mesh']['grid_'+idFast].size
    nMid = h5_file['mesh']['grid_'+idMid].size
    if idSlow is not None:
        nSlow = h5_file['mesh']['grid_'+idSlow].size
    else:
        nSlow = 0

    return (nFast, nMid, nSlow)

def checkCartesianExact(h5_file):
    """Check initial state for exact polynomial/trigonometric fields"""

    scheme = h5_file['/'].attrs['type']
    idFast, idMid, idSlow = getOrdering(scheme)
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    def check_values(data, tag, sol):
        scale = h5_file['/physical/scale1d'][()]
        err = -1
        if idSlow is not None:
            gF = gFast[()]
            itS = np.nditer(gSlow, flags=['c_index'])
            while not itS.finished:
                s = itS[0]
                itM = np.nditer(gMid, flags=['c_index'])
                while not itM.finished:
                    m = itM[0]
#                    if tag == 'velocity_curl_x':
#                        print(data[itS.index, itM.index,:])
#                        print('--------------------------------')
#                        print(sol(s, m, gF, tag,scale))
                    err = max(err, np.abs(data[itS.index, itM.index,:] - sol(s, m, gF, tag, scale)).max())

                    # Check velocity field
                    itM.iternext()
                itS.iternext()
        else:
            gF = gFast[()]
            itM = np.nditer(gMid, flags=['c_index'])
            while not itM.finished:
                m = itM[0]
                if scheme in [b'TT']:
                    scale_f = h5_file['/physical/scale1d'][()]
                    scale_m = h5_file['/physical/scale2d'][()]

                    # Check temperature field
                    err = max(err, np.abs(data[itM.index,:] - sol(m,gF)).max())
                elif scheme in [b'TF']:
                    scale_m = h5_file['/physical/scale1d'][()]

                    # Check temperature field
                    err = max(err, np.abs(dataT[itM.index,:] - sol(m,gF)).max())
                itM.iternext()
        print("\t\tError: {:g}".format(err))

    def exact_sol(s, m, gF, tag, scale):
        aT = [-3.0,2.0,3.0,-1.0,5.0]
        mT = [2.0,1.0,2.0,1.0,3.0]
        aP = [6.0,-7.0,5.0,2.0,5.0]
        mP = [-5.0,4.0,3.0,-3.0,1.0]
        # Chebyshevs order 0 to 4
        T = [np.cos(i*np.arccos(s)) for i in range(0,5)]

        if tag == 'temperature':
            return s**10*np.sin(10.*m)*np.sin(10.*gF)

        elif tag == 'temperature_grad_x':
            return s**2*m**3*gF

        elif tag == 'temperature_grad_y':
            return s**2*m**3*gF

        elif tag == 'temperature_grad_z':
            return s**2*m**3*gF

        elif tag == 'velocity_x':
            val = 0.0
            for j in range(0,5):
                valJT = (np.cos(j*m) + np.sin(j*m))
                valJP = j*(-np.sin(j*m) + np.cos(j*m))
                for i in range(0,5):
                    # Toroidal component
                    valIT = i*(-np.sin(i*gF) + np.cos(i*gF))
                    val += aT[0]*T[0]*valIT*valJT
                    val += aT[1]*T[1]*valIT*valJT
                    val += aT[2]*T[2]*valIT*valJT
                    val += aT[3]*T[3]*valIT*valJT
                    val += aT[4]*T[4]*valIT*valJT

                    # Poloidal component
                    valIP = scale*(np.cos(i*gF) + np.sin(i*gF))
                    val += aP[1]*valIP*valJP
                    val += (4.0*aP[2]*s)*valIP*valJP
                    val += aP[3]*(-3.0 + 12.0*s*s)*valIP*valJP
                    val += aP[4]*(-16.0*s + 32.0*s*s*s)*valIP*valJP

            # Mean X component
            val += mT[0]*T[0]
            val += mT[1]*T[1]
            val += mT[2]*T[2]
            val += mT[3]*T[3]
            val += mT[4]*T[4]

            return val

        elif tag == 'velocity_y':
            val = 0.0
            for j in range(0,5):
                valJT = j*(-np.sin(j*m) + np.cos(j*m))
                valJP = (np.cos(j*m) + np.sin(j*m))
                for i in range(0,5):
                    # Toroidal component
                    valIT = (np.cos(i*gF) + np.sin(i*gF))
                    val -= aT[0]*T[0]*valIT*valJT
                    val -= aT[1]*T[1]*valIT*valJT
                    val -= aT[2]*T[2]*valIT*valJT
                    val -= aT[3]*T[3]*valIT*valJT
                    val -= aT[4]*T[4]*valIT*valJT

                    # Poloidal component
                    valIP = scale*i*(-np.sin(i*gF) + np.cos(i*gF))
                    val += aP[1]*valIP*valJP
                    val += (4.0*aP[2]*s)*valIP*valJP
                    val += aP[3]*(-3.0 + 12.0*s*s)*valIP*valJP
                    val += aP[4]*(-16.0*s + 32.0*s*s*s)*valIP*valJP

            # Mean Y component
            val += mP[0]*T[0]
            val += mP[1]*T[1]
            val += mP[2]*T[2]
            val += mP[3]*T[3]
            val += mP[4]*T[4]

            return val

        elif tag == 'velocity_z':
            val = 0.0
            for j in range(0,5):
                valJP = (np.cos(j*m) + np.sin(j*m))
                for i in range(0,5):
                    # poloidal component
                    valIP = -(i*i + j*j)*(np.cos(i*gF) + np.sin(i*gF))
                    val -= aP[0]*T[0]*valIP*valJP
                    val -= aP[1]*T[1]*valIP*valJP
                    val -= aP[2]*T[2]*valIP*valJP
                    val -= aP[3]*T[3]*valIP*valJP
                    val -= aP[4]*T[4]*valIP*valJP

            return val

        elif tag == 'velocity_curl_x':
            val = 0.0
            for j in range(0,5):
                valJT = j*(-np.sin(j*m) + np.cos(j*m))
                valJP = (np.cos(j*m) + np.sin(j*m))
                for i in range(0,5):
                    # Toroidal component
                    valIT = scale*(np.cos(i*gF) + np.sin(i*gF))
                    val += aT[1]*valIT*valJT
                    val += (4.0*aT[2]*s)*valIT*valJT
                    val += aT[3]*(-3.0 + 12.0*s*s)*valIT*valJT
                    val += aT[4]*(-16.0*s + 32.0*s*s*s)*valIT*valJT

                    # Poloidal component, part 1
                    valIP1 = -i*(i*i + j*j)*(-np.sin(i*gF) + np.cos(i*gF))
                    val -= aP[0]*T[0]*valIP1*valJP
                    val -= aP[1]*T[1]*valIP1*valJP
                    val -= aP[2]*T[2]*valIP1*valJP
                    val -= aP[3]*T[3]*valIP1*valJP
                    val -= aP[4]*T[4]*valIP1*valJP

                    # Poloidal component, part 2
                    valIP2 = scale*scale*i*(-np.sin(i*gF) + np.cos(i*gF))
                    val -= (4.0*aP[2])*valIP2*valJP
                    val -= (24.0*aP[3]*s)*valIP2*valJP
                    val -= (-16.0 + 96*s*s)*aP[4]*valIP2*valJP

            # Mean Y component
            val -= scale*mP[1]
            val -= scale*(4.0*mP[2]*s)
            val -= scale*mP[3]*(-3.0 + 12.0*s*s)
            val -= scale*mP[4]*(-16.0*s + 32.0*s*s*s)

            return val

        elif tag == 'velocity_curl_y':
            val = 0.0
            for j in range(0,5):
                valJT = (np.cos(j*m) + np.sin(j*m))
                valJP = j*(-np.sin(j*m) + np.cos(j*m))
                for i in range(0,5):
                    # Toroidal component
                    valIT = scale*i*(-np.sin(i*gF) + np.cos(i*gF))
                    val += aT[1]*valIT*valJT
                    val += (4.0*aT[2]*s)*valIT*valJT
                    val += aT[3]*(-3.0 + 12.0*s*s)*valIT*valJT
                    val += aT[4]*(-16.0*s + 32.0*s*s*s)*valIT*valJT

                    # Poloidal component, part 1
                    valIP1 = -(i*i + j*j)*(np.cos(i*gF) + np.sin(i*gF))
                    val += aP[0]*T[0]*valIP1*valJP
                    val += aP[1]*T[1]*valIP1*valJP
                    val += aP[2]*T[2]*valIP1*valJP
                    val += aP[3]*T[3]*valIP1*valJP
                    val += aP[4]*T[4]*valIP1*valJP

                    # Poloidal component, part 2
                    valIP2 = scale*scale*(np.cos(i*gF) + np.sin(i*gF))
                    val += (4.0*aP[2])*valIP2*valJP
                    val += (24.0*aP[3]*s)*valIP2*valJP
                    val += (-16.0 + 96*s*s)*aP[4]*valIP2*valJP

            # Mean X component
            val += scale*mT[1]
            val += scale*(4.0*mT[2]*s)
            val += scale*mT[3]*(-3.0 + 12.0*s*s)
            val += scale*mT[4]*(-16.0*s + 32.0*s*s*s)

            return val

        elif tag == 'velocity_curl_z':
            val = 0.0
            for j in range(0,5):
                valJT = (np.cos(j*m) + np.sin(j*m))
                for i in range(0,5):
                    # Toroidal component
                    valIT = -(i*i + j*j)*(np.cos(i*gF) + np.sin(i*gF))
                    val -= aT[0]*T[0]*valIT*valJT
                    val -= aT[1]*T[1]*valIT*valJT
                    val -= aT[2]*T[2]*valIT*valJT
                    val -= aT[3]*T[3]*valIT*valJT
                    val -= aT[4]*T[4]*valIT*valJT

            return val

    for f in ['temperature', 'temperature_grad', 'velocity', 'velocity_grad', 'velocity_curl']:
        if f in h5_file:
            print("Found " + f)
            if f+'/'+f in h5_file:
                print("\tChecking " + f)
                data = h5_file[f+'/'+f]
                check_values(data,f,exact_sol)
            else: 
                for id in [idFast, idMid, idSlow]:
                    if f+'/'+f+'_'+id in h5_file:
                        print("\tChecking " + f + '_' + id)
                        data = h5_file[f+'/'+f+'_'+id]
                        check_values(data,f+'_'+id,exact_sol)

def checkShellBenchmarkC0(h5_file):
    """Check initial state for sphere benchmark C0"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']

    phi = gFast[()]
    ro = h5_file['/physical/ro'][()]
    ri = h5_file['/physical/rratio'][()]*ro
    itS = np.nditer(gSlow, flags=['c_index'])
    errT, errVr, errVth, errVph = -np.ones(4)
    while not itS.finished:
        r = itS[0]
        x = 2.0*r - ri  - ro
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = max(errT, np.abs(dataT[itS.index, itM.index,:] - (21.0/np.sqrt(17920*np.pi)*(1.0 - 3*x**2 + 3.0*x**4 - x**6)*np.sin(theta)**4*np.cos(4.0*phi))).max())

            # Check velocity field
            errVr = max(errVr, np.abs(dataVr[itS.index, itM.index,:]).max())
            errVth = max(errVth, np.abs(dataVth[itS.index, itM.index,:]).max())
            errVph = max(errVph, np.abs(dataVph[itS.index, itM.index,:]).max())
            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in temperature field:     {:g}".format(errT))

def checkShellBenchmarkC1(h5_file):
    """Check initial state for sphere benchmark C2"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']
    dataBr = h5_file['/magnetic/magnetic_r']
    dataBth = h5_file['/magnetic/magnetic_theta']
    dataBph = h5_file['/magnetic/magnetic_phi']

    phi = gFast[()]
    ro = h5_file['/physical/ro'][()]
    ri = h5_file['/physical/rratio'][()]*ro
    itS = np.nditer(gSlow, flags=['c_index'])
    errT, errVr, errVth, errVph, errBr, errBth, errBph = -np.ones(7)
    while not itS.finished:
        r = itS[0]
        x = 2.0*r - ri  - ro
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = max(errT, np.abs(dataT[itS.index, itM.index,:] - (21.0/np.sqrt(17920*np.pi)*(1.0 - 3*x**2 + 3.0*x**4 - x**6)*np.sin(theta)**4*np.cos(4.0*phi))).max())

            # Check velocity field
            errVr = max(errVr, np.abs(dataVr[itS.index, itM.index,:]).max())
            errVth = max(errVth, np.abs(dataVth[itS.index, itM.index,:]).max())
            errVph = max(errVph, np.abs(dataVph[itS.index, itM.index,:]).max())

            # Check magnetic field
            scale = 1.0/np.sqrt(2.0)
            errBr = max(errBr, np.abs(dataBr[itS.index, itM.index,:] - (scale*5.0/8.0*(8.0*ro - 6.0*r - 2.0*ri**4/r**3)*np.cos(theta))).max())
            errBth = max(errBth, np.abs(dataBth[itS.index, itM.index,:] - (scale*5.0/8.0*(9.0*r - 8.0*ro - ri**4/r**3)*np.sin(theta))).max())
            errBph = max(errBph, np.abs(dataBph[itS.index, itM.index,:] - (scale*5.0*np.sin(np.pi*(r - ri))*np.sin(2.0*theta))).max())

            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in radial magnetic field: {:g}".format(errBr))
    print("Error in theta magnetic field:  {:g}".format(errBth))
    print("Error in phi magnetic field:    {:g}".format(errBph))
    print("Error in temperature field:     {:g}".format(errT))

def checkSphereBenchmarkC1(h5_file):
    """Check initial state for sphere benchmark C1"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    dataT = h5_file['/temperature/temperature']
    dataDTr = h5_file['/temperature_grad/temperature_grad_r']
    dataDTth = h5_file['/temperature_grad/temperature_grad_theta']
    dataDTph = h5_file['/temperature_grad/temperature_grad_phi']
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']

    phi = gFast[()]
    itS = np.nditer(gSlow, flags=['c_index'])
    errT, errDTr, errDTth, errDTph, errVr, errVth, errVph = -np.ones(7)
    while not itS.finished:
        r = itS[0]
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = max(errT, np.abs(dataT[itS.index, itM.index,:] - (0.5*(1 - r**2) + 1e-5/8.*np.sqrt(35/np.pi)*r**3*(1 - r**2)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())

            # Check temperature gradient
            errDTr = max(errDTr, np.abs(dataDTr[itS.index, itM.index,:] - (-r + 1e-5/8.*np.sqrt(35/np.pi)*(3.0*r**2 - 5.0*r**4)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())
            errDTth = max(errDTth, np.abs(dataDTth[itS.index, itM.index,:] - (-3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.cos(theta)*np.sin(theta)**2*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())
            errDTph = max(errDTph, np.abs(dataDTph[itS.index, itM.index,:] - (3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.sin(theta)**2*(np.cos(3.0*phi) - np.sin(3.0*phi)))).max())

            # Check velocity field
            errVr = max(errVr, np.abs(dataVr[itS.index, itM.index,:]).max())
            errVth = max(errVth, np.abs(dataVth[itS.index, itM.index,:]).max())
            errVph = max(errVph, np.abs(dataVph[itS.index, itM.index,:]).max())
            itM.iternext()
        itS.iternext()
    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in temperature field:     {:g}".format(errT))
    print("Error in radial temperature gradient:     {:g}".format(errDTr))
    print("Error in theta temperature gradient:     {:g}".format(errDTth))
    print("Error in phi temperature gradient:     {:g}".format(errDTph))

def checkSphereBenchmarkC2(h5_file):
    """Check initial state for sphere benchmark C2"""

    scheme = h5_file['/'].attrs['type']
    gFast, gMid, gSlow = getGrid(scheme, h5_file)

    # Temperature field and gradient
    dataT = h5_file['/temperature/temperature']
    dataDTr = h5_file['/temperature_grad/temperature_grad_r']
    dataDTth = h5_file['/temperature_grad/temperature_grad_theta']
    dataDTph = h5_file['/temperature_grad/temperature_grad_phi']
    # Velocity field and curl
    dataVr = h5_file['/velocity/velocity_r']
    dataVth = h5_file['/velocity/velocity_theta']
    dataVph = h5_file['/velocity/velocity_phi']
    dataCVr = h5_file['/velocity_curl/velocity_curl_r']
    dataCVth = h5_file['/velocity_curl/velocity_curl_theta']
    dataCVph = h5_file['/velocity_curl/velocity_curl_phi']
    # Magnetic field and curl
    dataBr = h5_file['/magnetic/magnetic_r']
    dataBth = h5_file['/magnetic/magnetic_theta']
    dataBph = h5_file['/magnetic/magnetic_phi']
    dataCBr = h5_file['/magnetic_curl/magnetic_curl_r']
    dataCBth = h5_file['/magnetic_curl/magnetic_curl_theta']
    dataCBph = h5_file['/magnetic_curl/magnetic_curl_phi']

    phi = gFast[()]
    itS = np.nditer(gSlow, flags=['c_index'])
    errT, errDTr, errDTth, errDTph, errVr, errVth, errVph, errCVr, errCVth, errCVph, errBr, errBth, errBph, errCBr, errCBth, errCBph = -np.ones(16)
    while not itS.finished:
        r = itS[0]
        itM = np.nditer(gMid, flags=['c_index'])
        while not itM.finished:
            theta = itM[0]
            # Check temperature field
            errT = max(errT, np.abs(dataT[itS.index, itM.index,:] - (0.5*(1 - r**2) + 1e-5/8.*np.sqrt(35/np.pi)*r**3*(1 - r**2)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())

            # Check temperature gradient
            errDTr = max(errDTr, np.abs(dataDTr[itS.index, itM.index,:] - (-r + 1e-5/8.*np.sqrt(35/np.pi)*(3.0*r**2 - 5.0*r**4)*np.sin(theta)**3*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())
            errDTth = max(errDTth, np.abs(dataDTth[itS.index, itM.index,:] - (-3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.cos(theta)*np.sin(theta)**2*(np.cos(3.0*phi) + np.sin(3.0*phi)))).max())
            errDTph = max(errDTph, np.abs(dataDTph[itS.index, itM.index,:] - (3e-5/8.*np.sqrt(35/np.pi)*r**2*(r**2 - 1.0)*np.sin(theta)**2*(np.cos(3.0*phi) - np.sin(3.0*phi)))).max())

            # Check velocity field
            errVr = max(errVr, np.abs(dataVr[itS.index, itM.index,:]).max())
            errVth = max(errVth, np.abs(dataVth[itS.index, itM.index,:] - (-10*r**2/(7.0*np.sqrt(3.0))*np.cos(theta)*(3.0*(-147.0 + 343.0*r**2 - 217.0*r**4 + 29*r**6)*np.cos(phi) + 14.0*(-9.0 - 125.0*r**2 + 39*r**4 + 27*r**6)*np.sin(phi)))).max())
            errVph = max(errVph, np.abs(dataVph[itS.index, itM.index,:] - (-5.0*r/5544.*(7.0*((43700.0 - 58113*r**2 - 15345*r**4 + 1881*r**6 + 20790*r**8)*np.sin(theta) + 1485*r**2*(-9.0 + 115.0*r**2 - 167.0*r**4 + 70*r**6)*np.sin(3.0*theta)) + 528.0*np.sqrt(3)*r*np.cos(2.0*theta)*(14.0*(-9.0 - 125.0*r**2 + 39.0*r**4 + 27.0*r**6)*np.cos(phi) + 3.0*(147.0 - 343.0*r**2 + 217*r**4 - 29.0*r**6)*np.sin(phi))))).max())
            
            # Check curl of velocity field
            errCVr = max(errCVr, np.abs(dataCVr[itS.index, itM.index,:] - ((1.0/1386.0)*(-35.0*(21850.0 + 99.0*r**2*(-361.0 + 785.0*r**2 - 1243.0*r**4 + 630.0*r**6))*np.cos(theta) + 495.0*r*(-105.0*r*(-9.0 + 115.0*r**2 - 167.0*r**4 + 70.0*r**6)*np.cos(3.0*theta) + 4.0*np.sqrt(3)*np.sin(2.0*theta)*(14.0*(-9.0 - 125.0*r**2 + 39.0*r**4 + 27.0*r**6)*np.cos(phi) + 3.0*(147.0 - 343.0*r**2 + 217.0*r**4 - 29.0*r**6)*np.sin(phi)))))).max())
            errCVth = max(errCVth, np.abs(dataCVth[itS.index, itM.index,:] - ((54625.0*np.sin(theta))/99.0 + 5.0/84.0*r*(21.0*r*((-1174.0 - 465.0*r**2 + 76.0*r**4 + 1050.0*r**6)*np.sin(theta) + 15.0*(-18.0 + 345.0*r**2 - 668.0*r**4 + 350.0*r**6)*np.sin(3.0*theta)) + 8.0*np.sqrt(3)*np.cos(2.0*theta)*(14.0*(-27.0 - 625.0*r**2 + 273.0*r**4 + 243.0*r**6)*np.cos(phi) + 3.0*(441.0 - 1715.0*r**2 + 1519.0*r**4 - 261.0*r**6)*np.sin(phi))))).max())
            errCVph = max(errCVph, np.abs(dataCVph[itS.index, itM.index,:] - (10.0*r/(7.0*np.sqrt(3))*np.cos(theta)*((1323.0 - 5145.0*r**2 + 4557.0*r**4 - 783.0*r**6)*np.cos(phi) + 14.0*(27.0 + 625.0*r**2 - 273.0*r**4 - 243.0*r**6)*np.sin(phi)))).max())
            

            # Check magnetic field
            errBr = max(errBr, np.abs(dataBr[itS.index, itM.index,:]).max())
            errBth = max(errBth, np.abs(dataBth[itS.index, itM.index,:] - (-3.0/2.0*r*(-1.0 + 4.0*r**2 - 6.0*r**4 + 3.0*r**6)*(np.cos(phi) + np.sin(phi)))).max())
            errBph = max(errBph, np.abs(dataBph[itS.index, itM.index,:] - (-3.0/4.0*r*(-1.0 + r**2)*np.cos(theta)*(3.0*r*(2.0 - 5.0*r**2 + 4.0*r**4)*np.sin(theta) + 2.0*(1.0 - 3.0*r**2 + 3.0*r**4)*(np.cos(phi) - np.sin(phi))))).max())

            # Check curl of magnetic field
            errCBr = max(errCBr, np.abs(dataCBr[itS.index, itM.index,:] - (-9.0/8.0*r*(-2.0 + 7.0*r**2 - 9.0*r**4 + 4.0*r**6)*(1.0 + 3.0*np.cos(2.0*theta)) + 3.0*(-1.0 + 4*r**2 - 6.0*r**4 + 3.0*r**6)*np.sin(theta)*np.cos(phi) + 3.0*(1.0 - 4.0*r**2 + 6.0*r**4 - 3.0*r**6)*np.sin(theta)*np.sin(phi))).max())
            errCBth = max(errCBth, np.abs(dataCBth[itS.index, itM.index,:] - (3.0/4.0*np.cos(theta)*(3.0*r*(-6.0 + 35.0*r**2 - 63.0*r**4 + 36.0*r**6)*np.sin(theta) + 4.0*(-1.0 + 8.0*r**2 - 18.0*r**4 + 12.0*r**6)*(np.cos(phi) - np.sin(phi))))).max())
            errCBph = max(errCBph, np.abs(dataCBph[itS.index, itM.index,:] - (-3.0*(-1.0 + 8.0*r**2 - 18.0*r**4 + 12.0*r**6)*(np.cos(phi) + np.sin(phi)))).max())

            itM.iternext()
        itS.iternext()

    print("Error in radial velocity field: {:g}".format(errVr))
    print("Error in theta velocity field:  {:g}".format(errVth))
    print("Error in phi velocity field:    {:g}".format(errVph))
    print("Error in radial curl of velocity field: {:g}".format(errCVr))
    print("Error in theta curl of velocity field:  {:g}".format(errCVth))
    print("Error in phi curl of velocity field:    {:g}".format(errCVph))

    print("Error in radial magnetic field: {:g}".format(errBr))
    print("Error in theta magnetic field:  {:g}".format(errBth))
    print("Error in phi magnetic field:    {:g}".format(errBph))
    print("Error in radial curl of magnetic field: {:g}".format(errCBr))
    print("Error in theta curl of magnetic field:  {:g}".format(errCBth))
    print("Error in phi curl of magnetic field:    {:g}".format(errCBph))

    print("Error in temperature field:     {:g}".format(errT))
    print("Error in radial temperature gradient:     {:g}".format(errDTr))
    print("Error in theta temperature gradient:     {:g}".format(errDTth))
    print("Error in phi temperature gradient:     {:g}".format(errDTph))

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:")
    except getopt.GetoptError:
        print('check_physical.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_physical.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg

    # Read HDF5 file
    h5_file = h5py.File(inputfile, 'r')
    # Check state
    #checkShellBenchmarkC0(h5_file)
    #checkShellBenchmarkC1(h5_file)
    #checkSphereBenchmarkC1(h5_file)
    #checkSphereBenchmarkC2(h5_file)
    checkCartesianExact(h5_file)

    # Close file
    h5_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
