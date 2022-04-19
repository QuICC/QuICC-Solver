"""Check accuracy for radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import quicc.transform.cylinder as transf
import quicc.geometry.cylindrical.cylinder_radius as cylinder


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

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
    print("solution:")
    print(sol)
    print("computation:")
    print(rhs)
    print("-----------------------------------------")
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, parity, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torcheb(x_to_phys(sol_expr,grid), psol)
    print("solution:")
    print(sol)
    print("computation:")
    print(lhs)
    print("-----------------------------------------")
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

def test_backward_galerkin(opA, opB, opS, parity, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torcheb(x_to_phys(sol_expr,grid), psol)
    print("solution:")
    print(sol)
    print("computation:")
    print(lhs)
    print("-----------------------------------------")
    err = np.abs(opS*lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax galerkin backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax galerkin backward relative error: " + str(np.max(relerr)))

def zblk(nr, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.zblk(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)])
        ssol = 0
        test_forward(A, parity, sphys, ssol, rg, 0)

def x1(nr, rg):
    """Accuracy test for x1 operator"""

    print("x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1(nr, parity, cylinder.radbc.no_bc(), zr = 0)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        sphys = 1e10*x**(2*nr-2+parity)
        ssol = sy.expand(x*sphys)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def x2(nr, rg):
    """Accuracy test for x2 operator"""

    print("x2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x2(nr, parity, cylinder.radbc.no_bc(), zr = 0)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        sphys = 1e10*x**(2*nr-2+parity)
        ssol = sy.expand(x**2*sphys)
        test_forward(A, parity, sphys, ssol, rg, 1)

def d1(nr, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x))
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def x1d1(nr, rg):
    """Accuracy test for x1d1 operator"""

    print("x1d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1d1(nr, parity, cylinder.radbc.no_bc(), zr = 0)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        print(sphys)
        ssol = sy.expand(x*sy.diff(sphys,x))
        test_forward(A, parity, sphys, ssol, rg, 1)

def x1div(nr, rg):
    """Accuracy test for x1div operator"""

    print("x1div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1div(nr, parity, cylinder.radbc.no_bc(), zr = 0)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        print(sphys)
        ssol = sy.expand(sy.diff(x*sphys,x))
        test_forward(A, parity, sphys, ssol, rg, 1)

def i1(nr, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.integrate(sphys,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1d1(nr, rg):
    """Accuracy test for i1x1d1 operator"""

    print("i1x1d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1d1(nr, parity, cylinder.radbc.no_bc())
        print(A.todense())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1div(nr, rg):
    """Accuracy test for i1x1div operator"""

    print("i1x1div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1div(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1(nr, rg):
    """Accuracy test for i1x1 operator"""

    print("i1x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i1x2(nr, rg):
    """Accuracy test for i1x2 operator"""

    print("i1x2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i2(nr, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.integrate(sphys,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i2x1(nr, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i2x2d1(nr, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 2)

def i2x2d2(nr, rg):
    """Accuracy test for i2x2d2 operator"""

    print("i2x2d2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2d2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i2x3d1(nr, rg):
    """Accuracy test for i2x3d1 operator"""

    print("i2x3d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x3d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x**3)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i2x3d1x_2(nr, rg):
    """Accuracy test for i2x3d1x_2 operator"""

    print("i2x3d1x_2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x3d1x_2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x - 2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i2x2(nr, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i2x2div(nr, rg):
    """Accuracy test for i2x2div operator"""

    print("i2x2div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2div(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.expand(x*sy.diff(ssol,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i2x2laplh(nr, rg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

    print("\tbc = 10:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\t\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, {0:10}).tocsr()
        B = cylinder.i2(nr, parity, cylinder.radbc.no_bc()).tocsr()
        ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*(nr-1),2)]))
        sphys = sy.expand(x**2*sy.diff(ssol,x,x) + x*sy.diff(ssol,x) - m**2*ssol)
        test_backward_tau(A, B, parity, sphys, ssol, rg)

    print("\tbc = 11:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\t\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, {0:11}).tocsr()
        B = cylinder.i2(nr, parity, cylinder.radbc.no_bc()).tocsr()
        ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*(nr-2),2)]))
        sphys = sy.expand(x**2*sy.diff(ssol,x,x) + x*sy.diff(ssol,x) - m**2*ssol)
        test_backward_tau(A, B, parity, sphys, ssol, rg)

    print("\tbc = -10:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\t\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, {0:-10, 'rt':1}).tocsr()
        B = cylinder.i2(nr, parity, cylinder.radbc.no_bc()).tocsr()
        B = B[1:,:]
        S = cylinder.radbc.stencil(nr, parity, {0:-10})
        ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*(nr-1),2)]))
        sphys = sy.expand(x**2*sy.diff(ssol,x,x) + x*sy.diff(ssol,x) - m**2*ssol)
        test_backward_galerkin(A, B, S, parity, sphys, ssol, rg)

    print("\tbc = -11:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\t\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, {0:-11, 'rt':1}).tocsr()
        B = cylinder.i2(nr, parity, cylinder.radbc.no_bc()).tocsr()
        B = B[1:,:]
        S = cylinder.radbc.stencil(nr, parity, {0:-11})
        ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*(nr-2),2)]))
        sphys = sy.expand(x**2*sy.diff(ssol,x,x) + x*sy.diff(ssol,x) - m**2*ssol)
        test_backward_galerkin(A, B, S, parity, sphys, ssol, rg)

def i2x3laplhx_1(nr, rg):
    """Accuracy test for i2x3laplhx_1 operator"""

    print("i2x3laplhx_1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x3laplhx_1(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) - x*sy.diff(sphys,x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i4(nr, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4(nr, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4laplh(nr, rg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4laplh(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4lapl2h(nr, rg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4lapl2h(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def qid(nr, rg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        A = cylinder.qid(nr, parity, 3, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(m%2,2*nr,2)])
        ssol = sphys
        test_forward(A, parity, sphys, ssol, rg, 3)

def mvss_1_x(nr, rg):
    """Accuracy test for 1/r operator vs r solve operator"""

    print("mvss_1_x:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        Am = cylinder.x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys/x)
        test_forward(Am, (parity,(parity+1)%2), sphys, ssol, rg, 0)
        As = cylinder.x1(nr, parity, cylinder.radbc.no_bc())
        As = As.tocsr()
        Bs = cylinder.qid(nr, parity, 0, cylinder.radbc.no_bc())
        test_backward_tau(As, Bs, ((parity+1)%2,parity), sphys, ssol, rg)

def mvss_d1(nr, rg):
    """Accuracy test for d1 operator vs quasi-inverse solve"""

    print("mvss_d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        Am = cylinder.d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x))
        test_forward(Am, (parity,(parity+1)%2), sphys, ssol, rg, 0)
        if (nr%2 == 0 and parity == 1) or (nr%2 == 1 and parity == 0):
            As = cylinder.i1(nr, parity, {0:99})
            Bs = cylinder.qid(nr, parity, 1, cylinder.radbc.no_bc())
        else:
            As = cylinder.i1(nr, parity, cylinder.radbc.no_bc())
            Bs = cylinder.qid(nr, parity, 0, cylinder.radbc.no_bc())
        As = As.tocsr()
        test_backward_tau(As, Bs, ((parity+1)%2,parity), sphys, ssol, rg)


if __name__ == "__main__":
    # Set test parameters
    nr = 16
    rg = transf.rgrid(nr)

    # run tests
    #zblk(nr, rg)
#    x1(nr, rg)
#    x2(nr, rg)
#    d1(nr, rg)
#    x1d1(nr, rg)
#    x1div(nr, rg)
#    i1(nr, rg)
#    i1x1d1(nr, rg)
#    i1x1div(nr, rg)
#    i1x1(nr, rg)
#    i1x2(nr, rg)
#    i2(nr, rg)
#    i2x1(nr, rg)
#    i2(nr, rg)
#    i2x1(nr, rg)
#    i2x2(nr, rg)
#    i2x2d1(nr, rg)
#    i2x2d2(nr, rg)
#    i2x3d1(nr, rg)
#    i2x3d1x_2(nr, rg)
#    i2x2div(nr, rg)
    i2x2laplh(nr, rg)
#    i2x3laplhx_1(nr, rg)
#    i4(nr, rg)
#    i4x4(nr, rg)
#    i4x4laplh(nr, rg)
#    i4x4lapl2h(nr, rg)
#    qid(nr, rg)

    # Matmult vs Solve operators
#    mvss_1_x(nr, rg)
#    mvss_d1(nr, rg)
