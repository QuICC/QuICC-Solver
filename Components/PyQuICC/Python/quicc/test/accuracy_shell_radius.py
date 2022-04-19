"""Check accuracy for radial direction in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
#import mpmath as mp
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl

import quicc.transform.shell as transf
transf.min_r_points = 1
import quicc.geometry.spherical.shell_radius as shell

#mp.mp.dps = 200


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.lambdify(x, expr, "numpy")
    return func(grid)

def test_bc(op, res_expr, sol, grid):
    """Perform a boundary condition test"""

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid))
    lhs = lhs[0:op.shape[1]]
    rhs = op*lhs
    err = np.abs(rhs - sol)
    if np.max(err) > 10*np.spacing(1):
        print(rhs)
        print(err)
    print("\t\tMax BC error: " + str(np.max(err)))

def test_integral(op, res_expr, a, b, grid):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid))
    lhs = lhs[0:op.shape[1]]
    val = (op*lhs)[0]
    ref = sy.integrate(res_expr, (x,b-a,a+b))
    print("\t\tIntegral error: " + str(np.abs(val-ref)))

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

def test_backward_projector(opA, opB, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid))
    rhs = rhs[0:opA.shape[0]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    pad = np.zeros(grid.shape)
    pad[0:opA.shape[0]] = lhs
    lhs = transf.torphys(pad)
    sol = x_to_phys(sol_expr,grid)
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    #if np.max(err) > 10*np.spacing(1):
    #    print(err)
    print("\t\tMax backward projector error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax backward projector relative error: " + str(np.max(relerr)))

def test_backward_galerkin(opA, opB, opS, res_expr, sol_expr, grid):
    """Perform a galerkin backward operation test"""

    print("\tBackward galkerin test")
    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid))
    rhs = rhs[0:opB.shape[1]]
    lhs = spsplin.spsolve(opA,opB*rhs, use_umfpack=True)
    sol = transf.torcheb(x_to_phys(sol_expr,grid))
    sol = sol[0:opB.shape[1]]
    tmpOp = opS[0:opS.shape[1],:]
    tmp = spsplin.spsolve(tmpOp, sol[0:opS.shape[1]])
    print(np.abs(opS*tmp - sol))
    err = np.abs(opS*lhs - sol)
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))

def all_bc(nr, a, b, rg):
    """Accuracy test for BC operators"""

    print("\t Test BC 20:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 2, {0:20}).tolil()[0:2,:]
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    #sphys = x**2
    fsol = sy.lambdify(x, sphys)
    sol = np.array([fsol(a+b), fsol(-a+b)])
    test_bc(A, sphys, sol, rg)

    print("\t Test BC 21:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 2, {0:21, 'c':{'a':a, 'b':b}}).tolil()[0:2,:]
    #sphys = x**2
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    fsol = sy.lambdify(x, sy.diff(sphys))
    sol = np.array([fsol(a+b), fsol(-a+b)])
    test_bc(A, sphys, sol, rg)

    print("\t Test BC 22:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 2, {0:22, 'c':{'a':a, 'b':b}}).tolil()[0:2,:]
    #sphys = x**5
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    fsol = sy.lambdify(x, x*sy.diff(sphys/x))
    sol = np.array([fsol(a+b), fsol(-a+b)])
    test_bc(A, sphys, sol, rg)

    print("\t Test BC 23:")
    x = sy.Symbol('x')
    l = 5
    A = shell.qid(nr, 2, {0:23, 'c':{'a':a, 'b':b, 'l':l}}).tolil()[0:2,:]
    #sphys = x**2
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    fsolO = sy.lambdify(x, sy.diff(sphys) + (l+1.0)*sphys/x)
    fsolI = sy.lambdify(x, sy.diff(sphys) - l*sphys/x)
    sol = np.array([fsolO(a+b),fsolI(-a+b)])
    test_bc(A, sphys, sol, rg)

    print("\t Test BC 40:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 4, {0:40, 'c':{'a':a, 'b':b}}).tolil()[0:4,:]
    #sphys = x**7
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    fsolV = sy.lambdify(x, sphys)
    fsolD = sy.lambdify(x, sy.diff(sphys))
    sol = np.array([fsolV(a+b), fsolD(a+b), fsolV(-a+b), fsolD(-a+b)])
    test_bc(A, sphys, sol, rg)

    print("\t Test BC 41:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 4, {0:41, 'c':{'a':a, 'b':b}}).tolil()[0:4,:]
    #sphys = x**7
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    fsolV = sy.lambdify(x, sphys)
    fsolD = sy.lambdify(x, sy.diff(sphys,x,x))
    sol = np.array([fsolV(a+b), fsolD(a+b), fsolV(-a+b), fsolD(-a+b)])
    test_bc(A, sphys, sol, rg)

def integral(nr, a, b, rg):
    """Test the integration operator"""

    print("integral:")
    x = sy.Symbol('x')
    A = shell.integral(nr, a, b)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    test_integral(A, sphys, a, b, rg)


def zblk(nr, a, b, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    A = shell.zblk(nr, 0, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = 0
    test_forward(A, sphys, ssol, rg, 0)

def d1(nr, a, b, rg):
    """Accuracy test for d1 operator"""

    print("solve i1 (perm):")
    x = sy.Symbol('x')
    A = shell.i1(nr+1, a, b, {0:0, 'rt':1, 'cr':1})
    B = shell.qid(nr+1, 0, {0:0, 'rt':1, 'cr':1})
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*np.random.ranf()**i*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,nr,1)])
    ssol = sy.diff(sphys,x)
    test_backward_projector(A, B, sphys, ssol, rg)

def d2(nr, a, b, rg):
    """Accuracy test for d^2 operator"""

    print("solve i2 (perm):")
    x = sy.Symbol('x')
    A = shell.i2(nr+2, a, b, {0:0, 'rt':2, 'cr':2})
    B = shell.qid(nr+2, 0, {0:0, 'rt':2, 'cr':2})
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*np.random.ranf()**i*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,nr,1)])
    ssol = sy.diff(sphys,x,x)
    test_backward_projector(A, B, sphys, ssol, rg)

def divr(nr, a, b, rg):
    """Accuracy test for 1/x operator"""

    print("solve r (perm):")
    x = sy.Symbol('x')
    A = shell.r1(nr, a, b, {0:0})
    B = shell.qid(nr, 0, {0:0})
    #sphys = np.sum([(-1.0)**np.round(np.random.ranf())*np.random.ranf()**i*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,nr,1)])
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,np.ceil(2.0*nr/3.0),1)])
    ssol = sphys/x
    test_backward_projector(A, B, sphys, ssol, rg)

def divr2(nr, a, b, rg):
    """Accuracy test for 1/x^2 operator"""

    print("solve r^2 (perm):")
    x = sy.Symbol('x')
    A = shell.r2(nr, a, b, {0:0})
    B = shell.qid(nr, 0, {0:0})
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*np.random.ranf()**i*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,nr,1)])
    ssol = sphys/(x*x)
    test_backward_projector(A, B, sphys, ssol, rg)

def divr4(nr, a, b, rg):
    """Accuracy test for 1/x^4 operator"""

    print("solve r^4 (perm):")
    x = sy.Symbol('x')
    A = shell.r4(nr, a, b, {0:0})
    B = shell.qid(nr, 0, {0:0})
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*np.random.ranf()**i*sy.cos(int(i)*sy.acos((x-b)/a)) for i in np.arange(0,nr,1)])
    ssol = sphys/(x*x*x*x)
    test_backward_projector(A, B, sphys, ssol, rg)

def r1(nr, a, b, rg):
    """Accuracy test for r1 operator"""

    print("r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full chebyshev):")
    A = shell.r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def r2(nr, a, b, rg):
    """Accuracy test for x^2 operator"""

    print("r2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.r2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full Chebyshev):")
    A = shell.r2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def r4(nr, a, b, rg):
    """Accuracy test for x^4 operator"""

    print("r4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.r4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*x*x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full Chebyshev):")
    A = shell.r4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*x*x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def i1(nr, a, b, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full Chebyshev):")
    A = shell.i1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.simplify(sy.integrate(sphys,x))
    test_forward(A, sphys, ssol, rg, 1)

def i1r1(nr, a, b, rg):
    """Accuracy test for i1r1 operator"""

    print("i1r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(x*sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full Chebyshev):")
    A = shell.i1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.simplify(sy.integrate(x*sphys,x))
    test_forward(A, sphys, ssol, rg, 1)

def i2(nr, a, b, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2(nr, a, b, shell.radbc.no_bc()).tocsr()
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2r1(nr, a, b, rg):
    """Accuracy test for i2r1 operator"""

    print("i2r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(x*sphys,x,x)
    test_forward(A, sphys, ssol, rg, 0)

def i2r2d1(nr, a, b, rg):
    """Accuracy test for i2r2d1 operator"""

    print("i2r2d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2r2d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r2d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2r1d1r1(nr, a, b, rg):
    """Accuracy test for i2r1d1r1 operator"""

    print("i2r1d1r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2r1d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r1d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2r2(nr, a, b, rg):
    """Accuracy test for i2r2 operator"""

    print("i2r2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2r2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2r3(nr, a, b, rg):
    """Accuracy test for i2r3 operator"""

    print("i2r3:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2r3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**3)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2r2lapl(nr, a, b, rg):
    """Accuracy test for zblk operator"""

    print("i2r2lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r2lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\tbc = 20:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:20}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand((x-(a+b))*(x-(-a+b))*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -20:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:-20, 'rt':2}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-20, 'rt':2}).tocsr()
    ssol = (x-(a+b))*(x-(-a+b))*np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr-2,1)])
#    ssol = sy.expand((x-(a+b))*(x-(-a+b))*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 21:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:21}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -21:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:-21, 'rt':2}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-21, 'rt':2}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 22:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:22, 'c':{'a':a, 'b':b}}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -22:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r2lapl(nr, l, a, b, {0:-22, 'rt':2, 'c':{'a':a, 'b':b}}).tocsr()
    B = shell.i2r2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-22, 'rt':2, 'c':{'a':a, 'b':b}}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

def i2r3lapl(nr, a, b, rg):
    """Accuracy test for i2r3lapl operator"""

    print("i2r3lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2r3lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x) + 2*x**2*sy.diff(sphys,x) - l*(l+1)*sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2r3lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x) + 2*x**2*sy.diff(sphys,x) - l*(l+1)*sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i4(nr, a, b, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4(nr, a, b, shell.radbc.no_bc()).tocsr()
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r1(nr, a, b, rg):
    """Accuracy test for i4r1 operator"""

    print("i4r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r1d1r1(nr, a, b, rg):
    """Accuracy test for i4r1d1r1 operator"""

    print("i4r1d1r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r1d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r1d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r3(nr, a, b, rg):
    """Accuracy test for i4r3 operator"""

    print("i4r3:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r4d1(nr, a, b, rg):
    """Accuracy test for i4r4d1 operator"""

    print("i4r4d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r4d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r4d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r3d1r1(nr, a, b, rg):
    """Accuracy test for i4r3d1r1 operator"""

    print("i4r3d1r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r3d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r3d1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r3d2(nr, a, b, rg):
    """Accuracy test for i4r3d2 operator"""

    print("i4r3d2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r3d2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r3d2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r4(nr, a, b, rg):
    """Accuracy test for i4r4 operator"""

    print("i4r4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r4laplrd1r1(nr, a, b, rg):
    """Accuracy test for i4r4laplrd1r1 operator"""

    print("i4r4laplrd1r1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4r4laplrd1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x) + 3.0*x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r4laplrd1r1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x) + 3.0*x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r4lapl(nr, a, b, rg):
    """Accuracy test for i4r4lapl operator"""

    print("i4r4lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r4lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4r4lapl2(nr, a, b, rg):
    """Accuracy test for i4r4lapl2 operator"""

    print("i4r4lapl2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl2(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4r4lapl2(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\tbc = 40:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl2(nr, l, a, b, {0:40}).tocsr()
    B = shell.i4r4(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -40:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl2(nr, l, a, b, {0:-40, 'rt':4}).tocsr()
    B = shell.i4r4(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[4:,:]
    S = shell.radbc.stencil(nr, {0:-40, 'rt':4}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 41:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl2(nr, l, a, b, {0:41}).tocsr()
    B = shell.i4r4(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**3*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-6,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -41:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4r4lapl2(nr, l, a, b, {0:-41, 'rt':4}).tocsr()
    B = shell.i4r4(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[4:,:]
    S = shell.radbc.stencil(nr, {0:-41, 'rt':4}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**3*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-3,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

def qid(nr, a, b, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 3, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)

if __name__ == "__main__":
    # Set test parameters
    nr = 512
    a, b = shell.linear_r2x(1.0, 0.35)
    #a, b = shell.linear_r2x(20.0/13.0, 0.35)
    print((a, b))
    print((a+b, b-a))
    rg = transf.rgrid(2*nr, a, b)

    # run tests
#    all_bc(nr, a, b, rg)
#    integral(nr, a, b, rg)
#    zblk(nr, a, b, rg)
#    d1(nr, a, b, rg)
#    d2(nr, a, b, rg)
    divr(nr, a, b, rg)
#    divr2(nr, a, b, rg)
#    divr4(nr, a, b, rg)
#    r1(nr, a, b, rg)
#    r2(nr, a, b, rg)
#    r4(nr, a, b, rg)
#    i1(nr, a, b, rg)
#    i1r1(nr, a, b, rg)
#    i2(nr, a, b, rg)
#    i2r1(nr, a, b, rg)
#    i2r1d1r1(nr, a, b, rg)
#    i2r2d1(nr, a, b, rg)
#    i2r2(nr, a, b, rg)
#    i2r3(nr, a, b, rg)
#    i2r2lapl(nr, a, b, rg)
#    i2r3lapl(nr, a, b, rg)
#    i4(nr, a, b, rg)
#    i4r1(nr, a, b, rg)
#    i4r1d1r1(nr, a, b, rg)
#    i4r3d2(nr, a, b, rg)
#    i4r3(nr, a, b, rg)
#    i4r3d1r1(nr, a, b, rg)
#    i4r4d1(nr, a, b, rg)
#    i4r4(nr, a, b, rg)
#    i4r4laplrd1r1(nr, a, b, rg)
#    i4r4lapl(nr, a, b, rg)
#    i4r4lapl2(nr, a, b, rg)
#    qid(nr, a, b, rg)

#    divr(nr, a, b, rg)
#    divr2(nr, a, b, rg)
#    divr4(nr, a, b, rg)
