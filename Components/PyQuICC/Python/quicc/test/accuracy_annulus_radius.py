"""Check accuracy for radial direction in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import mpmath
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import quicc.transform.annulus as transf
import quicc.geometry.cylindrical.annulus_radius as annulus


rand_amplitude = 1e6

def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    mpmath.mp.dps = 200
    func = sy.lambdify(x, expr)
    return func(grid)

def test_forward(op, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid))
    lhs = lhs[0:op.shape[0]]
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torcheb(t)
    sol = sol[0:op.shape[0]]
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid))
    rhs = rhs[0:opA.shape[0]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torcheb(x_to_phys(sol_expr,grid))
    sol = sol[0:opA.shape[0]]
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

def zblk(nr, a, b, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.zblk(nr, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = 0
    test_forward(A, sphys, ssol, rg, 0)

def x1(nr, a, b, rg):
    """Accuracy test for x operator"""

    print("x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full chebyshev):")
    A = annulus.x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def x2(nr, a, b, rg):
    """Accuracy test for x^2 operator"""

    print("x2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full chebyshev):")
    A = annulus.x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def d1(nr, a, b, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.diff(sphys,x)
    test_forward(A, sphys, ssol, rg, 0)

def x1d1(nr, a, b, rg):
    """Accuracy test for x1d1 operator"""

    print("x1d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.x1d1(nr, a, b, annulus.radbc.no_bc(), zr = 0)
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    print(sphys)
    ssol = sy.expand(x*sy.diff(sphys,x))
    test_forward(A, sphys, ssol, rg, 0)

def x1div(nr, a, b, rg):
    """Accuracy test for x1div operator"""

    print("x1div:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.x1div(nr, a, b, annulus.radbc.no_bc(), zr = 0)
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    print(sphys)
    ssol = sy.expand(sy.diff(sphys*x,x))
    test_forward(A, sphys, ssol, rg, 0)

def i1(nr, a, b, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full chebyshev):")
    A = annulus.i1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

def i1x1d1(nr, a, b, rg):
    """Accuracy test for i1x1d1 operator"""

    print("i1x1d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i1x1d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x)
    ssol = sy.integrate(ssol,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full chebyshev):")
    A = annulus.i1x1d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x)
    ssol = sy.integrate(ssol,x)
    test_forward(A, sphys, ssol, rg, 1)

def i1x1div(nr, a, b, rg):
    """Accuracy test for i1x1div operator"""

    print("i1x1div:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i1x1div(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full chebyshev):")
    A = annulus.i1x1div(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x)
    test_forward(A, sphys, ssol, rg, 1)

def i1x1(nr, a, b, rg):
    """Accuracy test for i1x1 operator"""

    print("i1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i1x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.expand(sy.integrate(ssol,x))
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full chebyshev):")
    A = annulus.i1x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.expand(sy.integrate(ssol,x))
    test_forward(A, sphys, ssol, rg, 1)

def i1x2(nr, a, b, rg):
    """Accuracy test for i1x2 operator"""

    print("i1x2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i1x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.expand(sy.integrate(ssol,x))
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full chebyshev):")
    A = annulus.i1x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.expand(sy.integrate(ssol,x))
    test_forward(A, sphys, ssol, rg, 1)

def i2(nr, a, b, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x1(nr, a, b, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2d1(nr, a, b, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x2d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x2d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2d2(nr, a, b, rg):
    """Accuracy test for i2x2d2 operator"""

    print("i2x2d2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x2d2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x,x)*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x2d2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x,x)*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x3d1(nr, a, b, rg):
    """Accuracy test for i2x3d1 operator"""

    print("i2x3d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x3d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x**3)
    ssol = sy.expand(sy.integrate(ssol,x,x))
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x3d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x**3)
    ssol = sy.expand(sy.integrate(ssol,x,x))
    test_forward(A, sphys, ssol, rg, 2)

def i2x3d1x_2(nr, a, b, rg):
    """Accuracy test for i2x3d1x_2 operator"""

    print("i2x3d1x_2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x3d1x_2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x - 2*sphys)
    ssol = sy.expand(sy.integrate(ssol,x,x))
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x3d1x_2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x)*x - 2*sphys)
    ssol = sy.expand(sy.integrate(ssol,x,x))
    test_forward(A, sphys, ssol, rg, 2)

def i2x2(nr, a, b, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x2(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2div(nr, a, b, rg):
    """Accuracy test for i2x2div operator"""

    print("i2x2div:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i2x2div(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.expand(x*sy.diff(ssol,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x2div(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.expand(x*sy.diff(ssol,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2laplh(nr, a, b, rg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    m = np.random.randint(1, nr)
    A = annulus.i2x2laplh(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x2laplh(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x3laplhx_1(nr, a, b, rg):
    """Accuracy test for i2x3laplhx_1 operator"""

    print("i2x3laplhx_1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    m = np.random.randint(1, nr)
    A = annulus.i2x3laplhx_1(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) - x*sy.diff(sphys,x) - m**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full chebyshev):")
    A = annulus.i2x3laplhx_1(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) - x*sy.diff(sphys,x) - m**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i4x4(nr, a, b, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.i4x4(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full chebyshev):")
    A = annulus.i4x4(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4laplh(nr, a, b, rg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4laplh:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    m = np.random.randint(1, nr)
    A = annulus.i4x4laplh(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full chebyshev):")
    A = annulus.i4x4laplh(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl2h(nr, a, b, rg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    m = np.random.randint(1, nr)
    A = annulus.i4x4lapl2h(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full chebyshev):")
    A = annulus.i4x4lapl2h(nr, m, a, b, annulus.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def qid(nr, a, b, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = annulus.qid(nr, 3, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)

def mvss_1_x(nr, a, b, rg):
    """Accuracy test for 1/r operator vs r solve"""

    print("mvss_1_x:")
    x = sy.Symbol('x')
    Am = annulus.x1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys/x)
    test_forward(Am, sphys, ssol, rg, 0)
    As = annulus.x1(nr, a, b, annulus.radbc.no_bc())
    As = As.tocsr()
    Bs = annulus.qid(nr, 0, annulus.radbc.no_bc())
    test_backward_tau(As, Bs, sphys, ssol, rg)

def mvss_d1(nr, a, b, rg):
    """Accuracy test for derivative operator vs quasi inverse solve"""

    print("mvss_d1:")
    x = sy.Symbol('x')
    Am = annulus.d1(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    test_forward(Am, sphys, ssol, rg, 0)
    As = annulus.i1(nr, a, b, {0:99})
    As = As.tocsr()
    Bs = annulus.qid(nr, 1, annulus.radbc.no_bc())
    test_backward_tau(As, Bs, sphys, ssol, rg)

def mvss_x1div(nr, a, b, rg):
    """Accuracy test for radial divergence operator vs quasi inverse solve"""

    print("mvss_x1div:")
    x = sy.Symbol('x')
    Am = annulus.x1div(nr, a, b, annulus.radbc.no_bc())
    sphys = np.sum([rand_amplitude*np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sy.diff(sphys*x,x))
    test_forward(Am, sphys, ssol, rg, 0)
    As = annulus.i1(nr, a, b, {0:99})
    As = As.tocsr()
    Bs = annulus.qid(nr, 1, annulus.radbc.no_bc())*annulus.x1(nr, a, b, annulus.radbc.no_bc())
    test_backward_tau(As, Bs, sphys, ssol, rg)


if __name__ == "__main__":
    # Set test parameters
    nr = 14
    a, b = annulus.linear_r2x(1.0, 0.35)
    rg = transf.rgrid(2*nr, a, b)

    # run tests
    #zblk(nr, a, b, rg)
#    x1(nr, a, b, rg)
#    x2(nr, a, b, rg)
#    d1(nr, a, b, rg)
#    x1d1(nr, a, b, rg)
#    x1div(nr, a, b, rg)
    i1(nr, a, b, rg)
    i1x1d1(nr, a, b, rg)
    i1x1div(nr, a, b, rg)
    i1x1(nr, a, b, rg)
    i1x2(nr, a, b, rg)
    i2(nr, a, b, rg)
    i2x1(nr, a, b, rg)
    i2x2(nr, a, b, rg)
    i2x2d1(nr, a, b, rg)
    i2x2d2(nr, a, b, rg)
    i2x3d1(nr, a, b, rg)
    i2x3d1x_2(nr, a, b, rg)
    i2x2div(nr, a, b, rg)
    i2x2laplh(nr, a, b, rg)
    i2x3laplhx_1(nr, a, b, rg)
    i4x4(nr, a, b, rg)
    i4x4laplh(nr, a, b, rg)
    i4x4lapl2h(nr, a, b, rg)
#    qid(nr, a, b, rg)

    # Matmult vs Solve operators
#    mvss_1_x(nr, a, b, rg)
#    mvss_d1(nr, a, b, rg)
#    mvss_x1div(nr, a, b, rg)
