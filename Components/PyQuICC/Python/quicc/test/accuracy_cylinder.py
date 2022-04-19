"""Check accuracy for operators in a cylinder"""

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

import quicc.transform.cylinder as transf
import quicc.geometry.cylindrical.cylinder as cylinder


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

def test_forward(op, parity, res_expr, sol_expr, rg, zg, qr, qz):
    """Perform a forward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    print("\tForward test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    mesh = xz_to_phys(res_expr, rg, zg)
    lhs = transf.tocheb2d(mesh, pres)
    nr, nz = lhs.shape
    lhs = lhs.reshape(nr*nz, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape(nr,nz, order='F')
    mesh = xz_to_phys(sol_expr,rg, zg)
    sol = transf.tocheb2d(mesh, psol)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def test_backward_tau(opA, opB, parity, res_expr, sol_expr, rg, zg):
    """Perform a tau backward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    print("\tBackward tau test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    rhs = xz_to_phys(res_expr, rg, zg)
    rhs = transf.tocheb2d(rhs, pres)
    nr, nz = rhs.shape
    rhs = rhs.reshape(nr*nz, order = 'F')
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = lhs.reshape(nr,nz, order='F')
    sol = xz_to_phys(sol_expr, rg, zg)
    sol = transf.tocheb2d(sol, psol)
    err = np.abs(lhs - sol)
    vis_error(err, 'Tau backward error')
    print("\t\tMax tau backward error: " + str(np.max(err)))

def x1div(nr, nz, rg, zg):
    """Accuracy test for x1div operator"""

    print("x1div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1div(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(sy.diff(x*sphys,x))
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 1)

def x1e1(nr, nz, rg, zg):
    """Accuracy test for x1e1 operator"""

    print("x1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1e1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(sphys,z))
        test_forward(A, (parity, (parity+1)%2), sphys, ssol, rg, zg, 1, 1)

def i1j1(nr, nz, rg, zg):
    """Accuracy test for i1j1 operator"""

    print("i1j1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1j1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(sphys)
        ssol = sy.integrate(ssol,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z)
        test_forward(A, (parity, (parity+1)%2), sphys, ssol, rg, zg, 1, 1)

def i1j1x1d1(nr, nz, rg, zg):
    """Accuracy test for i1j1x1d1 operator"""

    print("i1j1x1d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1j1x1d1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z)
        test_forward(A, (parity, (parity+1)%2), sphys, ssol, rg, zg, 1, 1)

def i1j1x1div(nr, nz, rg, zg):
    """Accuracy test for i1j1x1div operator"""

    print("i1j1x1div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1j1x1div(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z)
        ssol = sy.expand(ssol)
        test_forward(A, (parity, (parity+1)%2), sphys, ssol, rg, zg, 1, 1)

def i1j1x1e1(nr, nz, rg, zg):
    """Accuracy test for i1j1x1e1 operator"""

    print("i1j1x1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1j1x1e1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(sphys,z))
        ssol = sy.integrate(ssol,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 1)

def i2j2(nr, nz, rg, zg):
    """Accuracy test for i2j2 operator"""

    print("i2j2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(sphys)
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 2)

def i2j2x1(nr, nz, rg, zg):
    """Accuracy test for i2j2x1 operator"""

    print("i2j2x1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sphys)
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, (parity, (parity+1)%2), sphys, ssol, rg, zg, 1, 2)

def i2j2x2(nr, nz, rg, zg):
    """Accuracy test for i2j2x2 operator"""

    print("i2j2x2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x2(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 2)

def i2j2x2d1(nr, nz, rg, zg):
    """Accuracy test for i2j2x2d1 operator"""

    print("i2j2x2d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x2d1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**2*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, zg, 1, 2)

def i2j2x2e1(nr, nz, rg, zg):
    """Accuracy test for i2j2x2e1 operator"""

    print("i2j2x2e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x2e1(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**2*sy.diff(sphys,z))
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 2)

def i2j2x2div(nr, nz, rg, zg):
    """Accuracy test for i2j2x2div operator"""

    print("i2j2x2div:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x2div(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, zg, 1, 2)

def i2x2laplh(nr, nz, rg, zg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, nz, m, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(x*sy.diff(sphys,x),x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 0)

def i2j2x2lapl(nr, nz, rg, zg):
    """Accuracy test for i2j2x2lapl operator"""

    print("i2j2x2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2j2x2lapl(nr, nz, m, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x*sy.diff(x*sy.diff(sphys,x),x) - m**2*sphys + x**2*sy.diff(sphys,z,z))
        ssol = sy.integrate(ssol,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z)
        ssol = sy.expand(ssol)
        test_forward(A, parity, sphys, ssol, rg, zg, 1, 2)

def i4j4(nr, nz, rg, zg):
    """Accuracy test for i4j4 operator"""

    print("i4j4:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4j4(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 2, 4)

def i4j4x4(nr, nz, rg, zg):
    """Accuracy test for i4j4x4 operator"""

    print("i4j4x4:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4j4x4(nr, nz, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**4*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 2, 4)

def i4x4laplh(nr, nz, rg, zg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4laplh(nr, nz, m, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**3*sy.diff(x*sy.diff(sphys,x),x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        ssol = sy.expand(ssol)
        test_forward(A, parity, sphys, ssol, rg, zg, 2, 0)

def i4j4x4lapl(nr, nz, rg, zg):
    """Accuracy test for i4j1e1 operator"""

    print("i4j4x4lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4j4x4lapl(nr, nz, m, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**3*sy.diff(x*sy.diff(sphys,x),x) - m**2*x**2*sphys + x**4*sy.diff(sphys,z,z))
        ssol = sy.integrate(ssol,x,x,x,x)
        ssol = sy.expand(ssol)
        ssol = sy.integrate(ssol,z,z,z,z)
        test_forward(A, parity, sphys, ssol, rg, zg, 2, 4)

def i4x4lapl2h(nr, nz, rg, zg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4lapl2h(nr, nz, m, parity, cylinder.cylbc.no_bc())
        sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)]) for j in np.arange(0,nz,1)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2.0*x**3*sy.diff(sphys,x,x,x) - (1.0 + 2.0*m**2)*x**2*sy.diff(sphys,x,x) + (1.0 + 2.0*m**2)*x*sy.diff(sphys,x) + (m**2 - 4.0)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        ssol = sy.expand(ssol)
        test_forward(A, parity, sphys, ssol, rg, zg, 2, 0)


if __name__ == "__main__":
    # Set test parameters
    nr = 24
    nz = 24
    rg = transf.rgrid(nr)
    zg = transf.zgrid(nz)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    x1div(nr, nz, rg, zg)
    x1e1(nr, nz, rg, zg)
    i1j1(nr, nz, rg, zg)
    i1j1x1d1(nr, nz, rg, zg)
    i1j1x1div(nr, nz, rg, zg)
    i1j1x1e1(nr, nz, rg, zg)
    i2j2(nr, nz, rg, zg)
    i2j2x1(nr, nz, rg, zg)
    i2j2x2(nr, nz, rg, zg)
    i2j2x2d1(nr, nz, rg, zg)
    i2j2x2e1(nr, nz, rg, zg)
    i2j2x2div(nr, nz, rg, zg)
    i2x2laplh(nr, nz, rg, zg)
    i2j2x2lapl(nr, nz, rg, zg)
    i4j4(nr, nz, rg, zg)
    i4j4x4(nr, nz, rg, zg)
    i4x4laplh(nr, nz, rg, zg)
    i4j4x4lapl(nr, nz, rg, zg)
    i4x4lapl2h(nr, nz, rg, zg)
