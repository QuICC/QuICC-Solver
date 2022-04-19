"""Check accuracy for radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import quicc.transform.sphere as transf
transf.min_r_points = 1
import quicc.geometry.spherical.sphere_radius as sphere


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    try:
        n = len(expr.free_symbols)
        if n == 0:
            return expr*np.ones(grid.shape)
        else:
            x = sy.Symbol('x')
            func = sy.utilities.lambdify(x, expr)
            return func(grid)
    except:
        return expr*np.ones(grid.shape)

def test_bc(op, parity, res_expr, sol_expr, point, grid):
    """Perform a boundary condition test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    rhs = op*lhs
    sol = x_to_phys(sol_expr, point)
    err = np.abs(rhs[0] - sol)
    relerr = err/np.abs(1 + sol)
    print((rhs[0], sol, err, relerr))

def test_integral(op, parity, res_expr, grid):
    """Perform a forward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    rhs = op*lhs
    lhs = lhs[0:op.shape[1]]
    val = (op*lhs)[0]
    ref = sy.integrate(res_expr, (x,0,1))
    print("\t\tIntegral error: " + str(np.abs(val-ref)))

def test_forward(op, parity, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torcheb(t, psol)
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, parity, res_expr, sol_expr, grid, fix = None):
    """Perform a tau backward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    if fix == 'DIVR':
        if pres == 0:
            n = rhs.shape[0]
            rhs[0] = 0
            for i in range(1,n):
                rhs[0] += (-1)**(i+1)*2.0*rhs[i]
        lhs = spsplin.spsolve(opA,opB*rhs)
#    elif fix == 'DIVR2':
#        err = 0.0
#        n = lhs.shape[0]
#        if pres == 1:
#            for i,v in enumerate(lhs[::-1]):
#                err += (-1)**(n-i+1)*v*2.0
#            err -= lhs[0]
#            #rhs[0] -= err/2.0
#        lhs = spsplin.spsolve(opA[0],opB[0]*rhs)
#
#        if pres == 0:
#            n = lhs.shape[0]
#            for i in range(0,lhs.shape[0]):
#                lhs[i] += (-1)**(n+i)*lhs[-1]
#        lhs = spsplin.spsolve(opA[1],opB[1]*lhs)
#        if (pres+1)%2 == 0:
#            n = lhs.shape[0]
#            err = lhs[-2]
#            for i in range(0,lhs.shape[0]):
#                lhs[i] += (-1)**(n+i+1)*err
    else:
        lhs = spsplin.spsolve(opA,opB*rhs)
#    if fix == 'DIVR':
#        if pres == 0:
#            n = lhs.shape[0]
#            for i in range(0,lhs.shape[0]):
#                lhs[i] += (-1)**(n+i)*lhs[-1]
#    elif fix == 'DIVR2':
#        n = lhs.shape[0]
#        for i in range(0,n):
#            lhs[i] += (-1)**(i+n)*(n-i)*lhs[-1]
    print(lhs)
    print('-------------------------------------------------------------')
    sol = transf.torcheb(x_to_phys(sol_expr,grid), psol)
    print(sol)
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

def integral(nr, rg):
    """Accuracy test for integral operator"""

    print("integral:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.integral(nr, l)
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        test_integral(A, l%2, sphys, rg)

def zblk(nr, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, 0, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = 0.0
        test_forward(A, l%2, sphys, ssol, rg, 0)

def bc(nr, rg):
    """Accuracy test for boundary conditions"""

    print("BC 10:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, l, {0:10})
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        test_bc(A, l%2, sphys, sphys, 1.0, rg)

    print("BC 11:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, l, {0:11})
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys))
        test_bc(A, l%2, sphys, ssol, 1.0, rg)

    print("BC 12:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, l, {0:12})
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys) - 1/x*sphys)
        test_bc(A, l%2, sphys, ssol, 1, rg)

    print("BC 13:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, l, {0:13, 'c':{'l':l}})
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys) + (l+1.0)/x*sphys)
        test_bc(A, l%2, sphys, ssol, 1, rg)

    print("BC 14:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, l, {0:14})
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x,x))
        test_bc(A, l%2, sphys, ssol, 1.0, rg)

def d1(nr, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.diff(sphys,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)

def r1(nr, rg):
    """Accuracy test for r1 operator"""

    x = sy.Symbol('x')
    print("r1:")
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.r1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x*sphys)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 0)

    print("solve 1/r:")
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        l = i
        print("\tTest for l = " + str(l))
        A = sphere.r1(nr, (l+1)%2, sphere.radbc.no_bc()).tocsr()
        B = sphere.qid(nr, l, 0, sphere.radbc.no_bc())
        sphys = np.sum([1e6*np.random.ranf()*x**(i) for i in np.arange(2*((l+1)%2) + l%2,2*nr-2,2)])
        ssol = sy.expand(sphys/x)
        error = 1e-4
        sphys += error
        test_backward_tau(A, B, (l%2,(l+1)%2), sphys, ssol, rg, fix = 'DIVR')

def r2(nr, rg):
    """Accuracy test for r2 operator"""

    x = sy.Symbol('x')
#    print("r2:")
#    for i in range(0,2):
#        l = np.random.randint(1, nr-1)
#        l = l + (l+i)%2
#        print("\tTest for l = " + str(l))
#        A = sphere.r2(nr, l, sphere.radbc.no_bc())
#        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
#        ssol = sy.expand(x*x*sphys)
#        test_forward(A, (l%2,l%2), sphys, ssol, rg, 0)

    print("solve 1/r^2:")
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        l = 2*((i+1)%2) + i
        print("\tTest for l = " + str(l))
        #A = (sphere.r1(nr, (l+1)%2, sphere.radbc.no_bc()).tocsr(), sphere.r1(nr, l%2, sphere.radbc.no_bc()).tocsr())
        #B = (sphere.qid(nr, l, 0, sphere.radbc.no_bc()),sphere.qid(nr, l, 0, sphere.radbc.no_bc()))
        A = sphere.r2(nr, l%2, sphere.radbc.no_bc()).tocsr()
        B = sphere.qid(nr, l, 0, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(2*((l+1)%2) + 3*(l%2),2*nr-2,2)])
        sphys = x*(x**2 + x**4)
        ssol = sy.expand(sphys/(x*x))
        error = 0e5*((i+1)%2 + i*x)
        sphys += error
        test_backward_tau(A, B, (l%2,l%2), sphys, ssol, rg, fix = 'DIVR2')

def i1(nr, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.integrate(sphys,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, l%1)

    print("solve D:")
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i1(nr, l, {0:991}).tocsr()
        B = sphere.qid(nr, l, l%2, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange((l+1)%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x))
        test_backward_tau(A, B, ((l+1)%2,l%2), sphys, ssol, rg)

def i2(nr, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.integrate(sphys,x,x)
        test_forward(A, (l%2,l%2), sphys, ssol, rg, 1)

    print("solve D^2:")
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2(nr, l, {0:992}).tocsr()
        B = sphere.qid(nr, l, 1, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x,x))
        test_backward_tau(A, B, (l%2,l%2), sphys, ssol, rg)

def i2r1(nr, rg):
    """Accuracy test for i2r1 operator"""

    print("i2r1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2r1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)

def i2r1d1r1(nr, rg):
    """Accuracy test for i2r1d1r1 operator"""

    print("i2r1d1r1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2r1d1r1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x*sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)

def i2r2d1(nr, rg):
    """Accuracy test for i2r2d1 operator"""

    print("i2r2d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2r2d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)


def i2r2(nr, rg):
    """Accuracy test for i2r2 operator"""

    print("i2r2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2r2(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)


def i2r2lapl(nr, rg):
    """Accuracy test for zblk operator"""

    print("i2r2lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2r2lapl(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)

def i4r3(nr, rg):
    """Accuracy test for i4r3 operator"""

    print("i4r3:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r3(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**3*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)

def i4r3d1r1(nr, rg):
    """Accuracy test for i4r3d1r1 operator"""

    print("i4r3d1r1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r3d1r1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**3*sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)

def i4r4d1(nr, rg):
    """Accuracy test for i4r4d1 operator"""

    print("i4r4d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r4d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)

def i4r4(nr, rg):
    """Accuracy test for i4r4 operator"""

    print("i4r4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r4(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def i4r4lapl(nr, rg):
    """Accuracy test for i4r4lapl operator"""

    print("i4r4lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r4lapl(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def i4r4lapl2(nr, rg):
    """Accuracy test for i4r4lapl2 operator"""

    print("i4r4lapl2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4r4lapl2(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def qid(nr, rg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        A = sphere.qid(nr, l, 3, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sphys
        test_forward(A, l%2, sphys, ssol, rg, 3)


if __name__ == "__main__":
    # Set test parameters
    nr = 64
    rg = transf.rgrid(nr)

    # run tests
#    integral(nr, rg)
#    zblk(nr, rg)
#    bc(nr, rg)
#    d1(nr, rg)
    r1(nr, rg)
#    r2(nr, rg)
#    i1(nr, rg)
#    i2(nr, rg)
#    i2r1(nr, rg)
#    i2r1d1r1(nr, rg)
#    i2r2d1(nr, rg)
#    i2r2(nr, rg)
#    i2r2lapl(nr, rg)
#    i4r3(nr, rg)
#    i4r3d1r1(nr, rg)
#    i4r4d1(nr, rg)
#    i4r4(nr, rg)
#    i4r4lapl(nr, rg)
#    i4r4lapl2(nr, rg)
#    qid(nr, rg)
