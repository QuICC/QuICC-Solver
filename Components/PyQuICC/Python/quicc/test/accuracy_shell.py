"""Check accuracy for spherical shell operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
if True:
    import matplotlib.pylab as pl
    has_error_plot = True
else:
    has_error_plot = False

import quicc.transform.shell as transf
import quicc.geometry.spherical.shell as shell

import scipy.io as io

def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        if has_error_plot:
            pl.imshow(np.log10(np.abs(err)))
            pl.title(title)
            pl.colorbar()
            pl.show()

def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

def test_forward(op, res_expr, sol_expr, rg, qx):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    nr = len(rg)
    nl = len(res_expr)
    phys = np.zeros((nr,nl))
    for i,expr in enumerate(res_expr):
        phys[:,i] = x_to_phys(expr,rg)
    lhs = transf.torcheb2D(phys)
    lhs = lhs.reshape(nr*nl, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape(nr, nl, order = 'F')
    rhs_r = transf.torphys2D(rhs)
    sol_r = np.zeros((nr, nl))
    for i,expr in enumerate(sol_expr):
        sol_r[:,i] = x_to_phys(expr, rg)
    sol = transf.torcheb2D(sol_r)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward point error')
    print("\t\tMax forward error: " + str(np.max(err)))

def i2x2(nr, maxl, a, b, rg):
    """Accuracy test for i2x2 operator (2D)"""

    print("i2x2:")
    x = sy.Symbol('x')
    m = np.random.randint(0, maxl//2)
    A = shell.i2x2(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**2*tphys)
        tsol = sy.integrate(tsol,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2coriolis(nr, maxl, a, b, rg):
    """Accuracy test for i2j2x2d1 operator (2D)"""

    print("i2x2coriolis:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i2x2coriolis(nr, maxl, m, a, b, shell.sphbc.no_bc())
    io.mmwrite('matrix_Q.mtx', A)
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**2*sy.diff(tphys,x))
        tsol = sy.integrate(tsol,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2lapl(nr, maxl, a, b, rg):
    """Accuracy test for i2j2x2lapl operator (2D)"""

    print("i2x2lapl:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i2x2lapl(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**2*sy.diff(tphys,x,x) + 2*x*sy.diff(tphys,x) - l*(l+1)*tphys)
        tsol = sy.integrate(tsol,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 2)

def i4x4(nr, maxl, a, b, rg):
    """Accuracy test for i4j4x4 operator (2D)"""

    print("i4x4:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i4x4(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**4*tphys)
        tsol = sy.integrate(tsol,x,x,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4coriolis(nr, maxl, a, b, rg):
    """Accuracy test for i4x4coriolis operator (2D)"""

    print("i4x4coriolis:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i4x4coriolis(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**3*sy.diff(x*sy.diff(tphys,x),x) - m**2*x**2*tphys)
        tsol = sy.integrate(tsol,x,x,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl(nr, maxl, a, b, rg):
    """Accuracy test for i4x4lapl operator (2D)"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i4x4lapl(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**4*sy.diff(tphys,x,x) + 2*x**3*sy.diff(tphys,x) - l*(l+1)*x**2*tphys)
        tsol = sy.integrate(tsol,x,x,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl2(nr, maxl, a, b, rg):
    """Accuracy test for i4x4lapl2 operator (2D)"""

    print("i4x4lapl2:")
    x = sy.Symbol('x')
    m = np.random.randint(1, maxl//2)
    A = shell.i4x4lapl2(nr, maxl, m, a, b, shell.sphbc.no_bc())
    sphys = []
    ssol = []
    for l in range(m, maxl+1):
        tphys = np.random.ranf()*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        tsol = sy.expand(x**4*sy.diff(tphys,x,x,x,x) + 4*x**3*sy.diff(tphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(tphys,x,x) + (l-1)*l*(l+1)*(l+2)*tphys)
        tsol = sy.integrate(tsol,x,x,x,x)
        sphys.append(tphys)
        ssol.append(tsol)
    test_forward(A, sphys, ssol, rg, 4)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    maxl = 11
    a, b = shell.rad.linear_r2x(1.0, 0.35)
    rg = transf.rgrid(nr, a, b)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    i2x2(nr, maxl, a, b, rg)
    i2x2coriolis(nr, maxl, a, b, rg)
    i2x2lapl(nr, maxl, a, b, rg)
    i4x4(nr, maxl, a, b, rg)
    i4x4coriolis(nr, maxl, a, b, rg)
    i4x4lapl(nr, maxl, a, b, rg)
    i4x4lapl2(nr, maxl, a, b, rg)
