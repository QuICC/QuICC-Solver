"""Check accuracy for cylindrical annulus operators"""

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

import quicc.transform.annulus as transf
import quicc.geometry.cylindrical.annulus as annulus


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        if has_error_plot:
            pl.imshow(np.log10(np.abs(err)))
            pl.title(title)
            pl.colorbar()
            pl.show()

def xz_to_phys(expr, rg, zg):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z), expr)
    vx, vz = np.meshgrid(rg, zg, indexing = 'ij')
    return func(vx, vz)

def test_forward(op, res_expr, sol_expr, rg, zg, qx, qz):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nr = len(rg)
    nz = len(zg)
    mesh = xz_to_phys(res_expr, rg, zg)
    lhs = transf.tocheb2d(mesh)
    lhs = lhs.reshape(nr*nz, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape(nr,nz, order='F')
    mesh = xz_to_phys(sol_expr,rg, zg)
    sol = transf.tocheb2d(mesh)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def test_backward_tau(opA, opB, res_expr, sol_expr, rg, zg):
    """Perform a tau backward operation test"""

    print("\tBackward tau test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nr = len(rg)
    nz = len(rg)
    rhs = xz_to_phys(res_expr, rg, zg)
    rhs = transf.tocheb2d(rhs)
    rhs = rhs.reshape(nr*nz, order = 'F')
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = lhs.reshape(nr,nz, order='F')
    sol = xz_to_phys(sol_expr, rg, zg)
    sol = transf.tocheb2d(sol)
    err = np.abs(lhs - sol)
    vis_error(err, 'Tau backward error')
    print("\t\tMax tau backward error: " + str(np.max(err)))

def x1div(nr, nz, a, b, rg, zg):
    """Accuracy test for x1div operator"""

    print("x1div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.x1div(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys*x,x))
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def x1e1(nr, nz, a, b, rg, zg):
    """Accuracy test for x1e1 operator"""

    print("x1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.x1e1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys*x,z))
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def i1j1(nr, nz, a, b, rg, zg):
    """Accuracy test for i1j1 operator"""

    print("i1j1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i1j1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def i1j1x1d1(nr, nz, a, b, rg, zg):
    """Accuracy test for i1j1x1d1 operator"""

    print("i1j1x1d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i1j1x1d1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def i1j1x1div(nr, nz, a, b, rg, zg):
    """Accuracy test for i1j1x1div operator"""

    print("i1j1x1div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i1j1x1div(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def i1j1x1e1(nr, nz, a, b, rg, zg):
    """Accuracy test for i1j1x1e1 operator"""

    print("i1j1x1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i1j1x1e1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, rg, zg, 1, 1)

def i2j2(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2 operator"""

    print("i2j2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2j2x1(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x1 operator"""

    print("i2j2x1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2x1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2j2x2(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x2 operator"""

    print("i2j2x2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2x2(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2j2x2d1(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x2d1 operator"""

    print("i2j2x2d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2x2d1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2j2x2e1(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x2e1 operator"""

    print("i2j2x2e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2x2e1(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2j2x2div(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x2div operator"""

    print("i2j2x2div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i2j2x2div(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sy.diff(x*sphys,x))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i2x2laplh(nr, nz, a, b, rg, zg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    m = np.random.randint(1, nr)
    A = annulus.i2x2laplh(nr, nz, m, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sy.diff(x*sy.diff(sphys,x),x) - m**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 2, 0)

def i2j2x2lapl(nr, nz, a, b, rg, zg):
    """Accuracy test for i2j2x2lapl operator"""

    print("i2j2x2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    m = np.random.randint(1, nr)
    A = annulus.i2j2x2lapl(nr, nz, m, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x*sy.diff(x*sy.diff(sphys,x),x) - m**2*sphys + x**2*sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 2, 2)

def i4j4(nr, nz, a, b, rg, zg):
    """Accuracy test for i4j4 operator"""

    print("i4j4:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i4j4(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, rg, zg, 4, 4)

def i4j4x4(nr, nz, a, b, rg, zg):
    """Accuracy test for i4j4x4 operator"""

    print("i4j4x4:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = annulus.i4j4x4(nr, nz, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**4*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, rg, zg, 4, 4)

def i4x4laplh(nr, nz, a, b, rg, zg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    m = np.random.randint(1, nr)
    A = annulus.i4x4laplh(nr, nz, m, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**3*sy.diff(x*sy.diff(sphys,x),x) - m**2*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 4, 0)

def i4j4x4lapl(nr, nz, a, b, rg, zg):
    """Accuracy test for i4j4x4lapl operator"""

    print("i4j4x4lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    m = np.random.randint(1, nr)
    A = annulus.i4j4x4lapl(nr, nz, m, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**3*sy.diff(x*sy.diff(sphys,x),x) - m**2*x**2*sphys + x**4*sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, rg, zg, 4, 4)

def i4x4lapl2h(nr, nz, a, b, rg, zg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    m = np.random.randint(1, nr)
    A = annulus.i4x4lapl2h(nr, nz, m, a, b, annulus.cylbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2.0*x**3*sy.diff(sphys,x,x,x) - (1.0 + 2.0*m**2)*x**2*sy.diff(sphys,x,x) + (1.0 + 2.0*m**2)*x*sy.diff(sphys,x) + (m**2 - 4.0)*m**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, rg, zg, 4, 0)


if __name__ == "__main__":
    # Set test parameters
    nr = 24
    nz = 24
    a, b = annulus.rad.linear_r2x(1.0, 0.35)
    rg = transf.rgrid(nr, a, b)
    zg = transf.zgrid(nz)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    x1div(nr, nz, a, b, rg, zg)
    x1e1(nr, nz, a, b, rg, zg)
    i1j1(nr, nz, a, b, rg, zg)
    i1j1(nr, nz, a, b, rg, zg)
    i1j1x1d1(nr, nz, a, b, rg, zg)
    i1j1x1div(nr, nz, a, b, rg, zg)
    i1j1x1e1(nr, nz, a, b, rg, zg)
    i2j2(nr, nz, a, b, rg, zg)
    i2j2x1(nr, nz, a, b, rg, zg)
    i2j2x2(nr, nz, a, b, rg, zg)
    i2j2x2d1(nr, nz, a, b, rg, zg)
    i2j2x2e1(nr, nz, a, b, rg, zg)
    i2j2x2div(nr, nz, a, b, rg, zg)
    i2x2laplh(nr, nz, a, b, rg, zg)
    i2j2x2lapl(nr, nz, a, b, rg, zg)
    i4j4(nr, nz, a, b, rg, zg)
    i4j4x4(nr, nz, a, b, rg, zg)
    i4x4laplh(nr, nz, a, b, rg, zg)
    i4j4x4lapl(nr, nz, a, b, rg, zg)
    i4x4lapl2h(nr, nz, a, b, rg, zg)
