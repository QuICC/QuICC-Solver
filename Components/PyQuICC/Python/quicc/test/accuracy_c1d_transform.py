"""Check accuracy for cartesian 1D operators"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import mpmath
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl
import scipy.io as io

import quicc.transform.cartesian as transf
import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_generic_1d as cg1d
import timeit

repTime = 100

def gaussian(x0, s, xg):
    """Compute Chebyshev expansion for Gaussian"""

    func = (1 - xg**2)*np.exp(-s*(xg-x0)**2)
    spectrum = transf.tocheb(func)
    return spectrum

def diff_gaussian(x0, s, xg):
    """Compute Chebyshev expansion for 1st derivative of the Gaussian"""

    func = 2.0*np.exp(-s*(xg - x0)**2)*(-(1.0 + s)*xg + s*xg**3 + s*x0 - s*xg**2*x0)
    spectrum = transf.tocheb(func)
    return spectrum

def diff2_gaussian(x0, s, xg):
    """Compute Chebyshev expansion for 2nd derivative of the Gaussian"""

    func = -2.0*np.exp(-s*(xg - x0)**2)*(1.0 + 2.0*s**2*(-1 + xg**2)*(xg - x0)**2 + s*(1.0 - 5.0*xg**2 + 4.0*xg*x0))
    spectrum = transf.tocheb(func)
    return spectrum

def diff4_gaussian(x0, s, xg):
    """Compute Chebyshev expansion for 4th derivative of the Gaussian"""

    func = -4.0*np.exp(-s*(xg - x0)**2)*s*(-6.0 + 4.0*s**3*(-1.0 + xg**2)*(xg - x0)**4 - 4.0*s**2*(xg - x0)**2*(-3.0 + 7.0*xg**2 - 4.0*xg*x0) + 3.0*s*(-1.0 + 13.0*xg**2 - 16.0*xg*x0 + 4.0*x0**2))
    spectrum = transf.tocheb(func)
    return spectrum


def d1(spec):
    """Compute 1st derivative of spectrum by matrix product"""

    A = c1d.d1(spec.size, c1d.c1dbc.no_bc()).tocsr()

    val = A*spec
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A*spec
    end = timeit.default_timer()
    spTime = end-start
    print("Timing d1 (sparse): " + str(spTime))

    A = A.todense()
    val = A.dot(spec)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A.dot(spec)
    end = timeit.default_timer()
    deTime = end-start
    print("Timing d1 (dense): " + str(deTime))
    return (val,spTime,deTime)

def d2(spec):
    """Compute 2nd derivative of spectrum by matrix product"""

    A = c1d.d2(spec.size, c1d.c1dbc.no_bc()).tocsr()

    val = A*spec
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A*spec
    end = timeit.default_timer()
    spTime = end-start
    print("Timing d2 (sparse): " + str(spTime))

    A = A.todense()
    val = A.dot(spec)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A.dot(spec)
    end = timeit.default_timer()
    deTime = end-start
    print("Timing d2 (dense): " + str(deTime))
    return (val,spTime,deTime)

def d4(spec):
    """Compute 4th derivative of spectrum by matrix product"""

    A = c1d.d4(spec.size, c1d.c1dbc.no_bc()).tocsc()

    val = A*spec
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A*spec
    end = timeit.default_timer()
    spTime = end-start
    print("Timing d4 (sparse): " + str(spTime))

    A = A.todense()
    val = A.dot(spec)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = A.dot(spec)
    end = timeit.default_timer()
    deTime = end-start
    print("Timing d4 (dense): " + str(deTime))
    return (val,spTime,deTime)

def i1(spec):
    """Compute 1st derivative of spectrum by quasi inverse solve"""

    A = c1d.i1(spec.size, {0:991}).tocsc()
    LU = spsplin.factorized(A)
    rhs = spec.copy()
    rhs[0] = 0.0

    val = LU(rhs)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = LU(rhs)
    end = timeit.default_timer()
    spTime = end - start
    print("Timing i1: " + str(spTime))
    return (val,spTime)

def i2(spec):
    """Compute 2nd derivative of spectrum by quasi inverse solve"""

    A = c1d.i2(spec.size, {0:992}).tocsc()
    LU = spsplin.factorized(A)
    rhs = spec.copy()
    rhs[0:2] = 0.0

    val = LU(rhs)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = LU(rhs)
    end = timeit.default_timer()
    spTime = end - start
    print("Timing i2: " + str(spTime))
    return (val,spTime)

def i4(spec):
    """Compute 4th derivative of spectrum by quasi inverse solve"""

    A = c1d.i4(spec.size, {0:994}).tocsc()
    LU = spsplin.factorized(A)
    rhs = spec.copy()
    rhs[0:4] = 0.0

    val = LU(rhs)
    start = timeit.default_timer()
    for i in range(0,repTime):
        val = LU(rhs)
    end = timeit.default_timer()
    spTime = end - start
    print("Timing i4: " + str(spTime))
    return (val,spTime)

if __name__ == "__main__":

    runs = [64,94,128,192,256,384,512,768,896,1024,1152,1280,1408,1536,1664,1792,1920,2048]
    timings_d1 = np.zeros((len(runs), 3))
    timings_d2 = np.zeros((len(runs), 3))
    timings_d4 = np.zeros((len(runs), 3))
    # Set test parameters
    for i, nx in enumerate(runs):
        print("------------------------------------------")
        print("Running for nx = " + str(nx))
        xg = transf.grid(nx)
        ref_gauss = gaussian(0.1, 500, xg)
        ref_d1 = diff_gaussian(0.1, 500, xg)
        ref_d2 = diff2_gaussian(0.1, 500, xg)
        ref_d4 = diff4_gaussian(0.1, 500, xg)
        mult_d1 = d1(ref_gauss)
        solve_d1 = i1(ref_gauss)
        timings_d1[i,0] = mult_d1[1]
        timings_d1[i,1] = mult_d1[2]
        timings_d1[i,2] = solve_d1[1]

        mult_d2 = d2(ref_gauss)
        solve_d2 = i2(ref_gauss)
        timings_d2[i,0] = mult_d2[1]
        timings_d2[i,1] = mult_d2[2]
        timings_d2[i,2] = solve_d2[1]

        mult_d4 = d4(ref_gauss)
        solve_d4 = i4(ref_gauss)
        timings_d4[i,0] = mult_d4[1]
        timings_d4[i,1] = mult_d4[2]
        timings_d4[i,2] = solve_d4[1]

#        pl.semilogy(abs(ref_gauss))
#        pl.semilogy(abs(ref_d1))
#        pl.semilogy(abs(ref_d2))
#        pl.semilogy(abs(ref_d4))
#        pl.show()
#
#        pl.semilogy(abs(ref_d1 - mult_d1[0]))
#        pl.semilogy(abs(ref_d1 - solve_d1[0]))
#        pl.show()
#
#        pl.semilogy(abs(ref_d2 - mult_d2[0]))
#        pl.semilogy(abs(ref_d2 - solve_d2[0]))
#        pl.show()
#
#        pl.semilogy(abs(ref_d4 - mult_d4[0]))
#        pl.semilogy(abs(ref_d4 - solve_d4[0]))
#        pl.show()

    aruns = np.array(runs)
    pl.plot(aruns/runs[0], timings_d1[:,0]/timings_d1[0,0], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d1[:,1]/timings_d1[0,1], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d1[:,2]/timings_d1[0,2], linewidth=3.0)
    pl.legend(['Sparse D','Dense D','Sparse D^{-1}'])
    pl.show()

    pl.plot(aruns/runs[0], timings_d2[:,0]//timings_d2[0,0], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d2[:,1]//timings_d2[0,1], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d2[:,2]//timings_d2[0,2], linewidth=3.0)
    pl.legend(['Sparse D^2','Dense D^2','Sparse D^{-2}'])
    pl.show()

    pl.plot(aruns/runs[0], timings_d4[:,0]/timings_d4[0,0], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d4[:,1]/timings_d4[0,1], linewidth=3.0)
    pl.plot(aruns/runs[0], timings_d4[:,2]/timings_d4[0,2], linewidth=3.0)
    pl.legend(['Sparse D^4','Dense D^4','Sparse D^{-4}'])
    pl.show()
