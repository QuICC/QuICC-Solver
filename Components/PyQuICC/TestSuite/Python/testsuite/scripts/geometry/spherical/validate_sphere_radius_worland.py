"""Validate Worland of Chebyshev type operators in sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import quicc.geometry.worland.setup as wsetup
#wsetup.type = "Chebyshev"
#wsetup.type = "Legendre"
#wsetup.type = "CylEnergy"
wsetup.type = "SphEnergy"

import quicc.transform.sphere_worland as transf
transf.min_r_points = 1
import quicc.geometry.spherical.sphere_radius_worland as sph

acceptable_error = 1e3*np.spacing(1)
acceptable_fwd_error = acceptable_error
acceptable_bwd_error = acceptable_error
input_spectrum_ratio = 3.0/4.0


def in_nn(nn):
    import math
    return math.ceil(nn*input_spectrum_ratio)

def rand_ls(nl, maxL):
    ls = list([])
    for i in range(0,nl):
        l = np.random.randint(1, maxL-1)
        l = l + (l+i)%2
        ls.append(l)
    return ls

def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    if expr == 0:
        arr = np.zeros(grid.shape[0])
    else:
        x = sy.Symbol('x')
        func = sy.utilities.lambdify(x, expr)
        arr = func(grid)
    return arr

def test_bc(op, l, res_expr, sol_expr, point, grid):
    """Perform a boundary condition test"""

    try:
        l_res, l_sol = l
    except:
        l_res = l
        l_sol = l_res

    n_sol, n_res = op.shape
    x = sy.Symbol('x')
    lhs = transf.torspec(x_to_phys(res_expr,grid), l_res, n_res)
    rhs = op*lhs
    sol = x_to_phys(sol_expr, point)
    err = np.abs(rhs[0,0] - sol)
    relerr = err/np.abs(1 + sol)
    print("\t"*2 + "Value: {:.12e}, Reference: {:.12e}".format(rhs[0,0],sol))
    print("\t"*2 + "Errors: {:.2e}, rel. Error: {:.2e}".format(err,relerr))

def test_forward(op, l, res_expr, sol_expr, grid, q, p):
    """Perform a forward operation test"""

    try:
        l_res, l_sol = l
    except:
        l_res = l
        l_sol = l_res

    n_sol, n_res = op.shape
    x = sy.Symbol('x')
    lhs = transf.torspec(x_to_phys(res_expr,grid), l_res, n_res)
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torspec(t, l_sol, n_sol+(p-q))[(p-q):]
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > acceptable_fwd_error:
        print("-"*30)
        print(err.T)
        print("-"*30)
    print("\t\tMax forward error: {:.2e}".format(np.max(err[q:])))
    if np.max(relerr[q:]) > acceptable_fwd_error:
        print("="*30)
        print(relerr.T)
        print("="*30)
    print("\t\tMax forward relative error: {:.2e}".format(np.max(relerr[q:])))

def test_backward_tau(opA, opB, parity, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    rhs = transf.torspec(x_to_phys(res_expr,grid), pres)
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torspec(x_to_phys(sol_expr,grid), psol)
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > acceptable_bwd_error:
        print("-"*30)
        print(err.T)
        print("-"*30)
    print("\t\tMax tau backward error: {:.2e}".format(np.max(err)))
    if np.max(relerr) > acceptable_bwd_error:
        print("="*30)
        print(relerr.T)
        print("="*30)
    print("\t\tMax tau backward relative error: {:.2e}".format(np.max(relerr)))

def zblk(nn, rg, ls = None):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, 0, sph.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,nn)])
        ssol = 0
        test_forward(A, l, sphys, ssol, rg, 0, 0)

def bc(nn, rg, ls = None):
    """Accuracy test for boundary conditions"""

    print("BC 10:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:10})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        test_bc(A, l, sphys, sphys, 1.0, rg)

    print("BC 11:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:11})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys))
        test_bc(A, l, sphys, ssol, 1.0, rg)

    print("BC 12:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:12})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys) - 1/x*sphys)
        test_bc(A, l, sphys, ssol, 1.0, rg)

    print("BC 13:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:13})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys) + (l+1.0)/x*sphys)
        test_bc(A, l, sphys, ssol, 1.0, rg)

    print("BC 14:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:14})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys,x,x))
        test_bc(A, l, sphys, ssol, 1.0, rg)

    print("BC 15:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:15})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys,x,x,x))
        test_bc(A, l, sphys, ssol, 1.0, rg)

    print("BC 16:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.zblk(nn, l, {0:16})
        sphys = np.sum([np.random.ranf()*x**(l+2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sy.expand(sy.diff(sphys,x,x,x,x))
        test_bc(A, l, sphys, ssol, 1.0, rg)

def r2(nn, rg, ls = None):
    """Accuracy test for r2 operator"""

    print("r2:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.r2(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = x*x*spoly*x**l
        test_forward(A, l, sphys, ssol, rg, 0, 0)

def i1(nn, rg, ls = None):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i1(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.integrate(4*x*spoly,x)*x**l
        test_forward(A, l, sphys, ssol, rg, 0, 1)
#
#    print("solve D:")
#    for i in range(0,2):
#        l = np.random.randint(1, nn-1)
#        l = l + (l+i)%2
#        print("\tTest for l = " + str(l))
#        A = sph.i1(nn, l, {0:991}).tocsr()
#        B = sph.qid(nn, l, 0, sph.radbc.no_bc())
#        sphys = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
#        ssol = sy.expand(sy.diff(sphys,x))
#        test_backward_tau(A, B, l, sphys, ssol, rg)

def i2(nn, rg, ls = None):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i2(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.integrate(4*x*sy.integrate(4*x*spoly,x),x)*x**l
        test_forward(A, l, sphys, ssol, rg, 1, 2)
#
#    print("solve D^2:")
#    for i in range(0,2):
#        l = np.random.randint(1, nn-1)
#        l = l + (l+i)%2
#        print("\tTest for l = " + str(l))
#        A = sph.i2(nn, l, {0:992}).tocsr()
#        B = sph.qid(nn, l, 1, sph.radbc.no_bc())
#        spoly = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,in_nn(nn),2)])
#        sphys = spoly*x**l
#        ssol = sy.expand(sy.diff(sphys,x,x))
#        test_backward_tau(A, B, l, sphys, ssol, rg)

def i2lapl(nn, rg, ls = None):
    """Accuracy test for zblk operator"""

    print("i2lapl:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i2lapl(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.expand(sy.diff(sphys,x,x) + 2*sy.diff(sphys,x)/x - l*(l+1)*sphys/x**2)*x**(-l)
        ssol = sy.integrate(4*x*sy.integrate(4*x*ssol,x),x)*x**l
        test_forward(A, l, sphys, ssol, rg, 1, 2)

def i4lapl(nn, rg, ls = None):
    """Accuracy test for i4lapl operator"""

    print("i4lapl:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i4lapl(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.expand(sy.diff(sphys,x,x) + 2*sy.diff(sphys,x)/x - l*(l+1)*sphys/x**2)*x**(-l)
        ssol = sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*ssol,x),x),x),x)*x**l
        test_forward(A, l, sphys, ssol, rg, 2, 4)

def i4lapl2(nn, rg, ls = None):
    """Accuracy test for i4lapl2 operator"""

    print("i4lapl2:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i4lapl2(nn, l, sph.radbc.no_bc())
        spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.expand(sy.diff(sphys,x,x,x,x) + 4*sy.diff(sphys,x,x,x)/x - 2*l*(l+1)*sy.diff(sphys,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*sphys/x**4)*x**(-l)
        ssol = sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*ssol,x),x),x),x)*x**l
        test_forward(A, l, sphys, ssol, rg, 2, 4)

def i4(nn, rg, ls = None, fixed_coeffs = False):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sph.i4(nn, l, sph.radbc.no_bc())
        if fixed_coeffs:
            spoly = np.sum([x**(2*i) for i in np.arange(0,in_nn(nn))])
        else:
            spoly = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        sphys = spoly*x**l
        ssol = sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*sy.integrate(4*x*spoly,x),x),x),x)*x**l
        test_forward(A, l, sphys, ssol, rg, 2, 4)

def qid(nn, rg, ls = None):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    if ls is None:
        ls = rand_ls(2, nn)
    for l in ls:
        A = sph.qid(nn, l, 3, sph.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(2*i) for i in np.arange(0,in_nn(nn))])
        ssol = sphys
        test_forward(A, l, sphys, ssol, rg, 3)


if __name__ == "__main__":
    # Set test parameters
    nn = 20
    nr = transf.nrgrid(nn)
    rg = transf.rgrid(nr)
    ls = None
    ls = [0,1,2]

    # run tests
#    zblk(nn, rg, ls)
    bc(nn, rg, ls)
    r2(nn, rg, ls)
    i1(nn, rg, ls)
#    i1qm(nn, rg, ls)
#    i1qp(nn, rg, ls)
    i2(nn, rg, ls)
    i2lapl(nn, rg, ls)
#    i2qm(nn, rg, ls)
#    i2qp(nn, rg, ls)
    i4(nn, rg, ls)
    i4lapl(nn, rg, ls)
    i4lapl2(nn, rg, ls)
#    i4qm(nn, rg, ls)
#    i4qp(nn, rg, ls)
#    qid(nn, rg, ls)
