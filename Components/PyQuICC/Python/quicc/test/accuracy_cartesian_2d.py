"""Check accuracy for cartesian 2D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import scipy.io as io
if True:
    import matplotlib.pylab as pl
    has_error_plot = True
else:
    has_error_plot = False

import quicc.transform.cartesian as transf
transf.min_x_points = 1
import quicc.geometry.cartesian.cartesian_2d as c2d


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        if has_error_plot:
            pl.imshow(np.log10(np.abs(err)))
            pl.title(title)
            pl.colorbar()
            pl.show()

def xz_to_phys(expr, grid_x, grid_z):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z), expr)
    vx, vz = np.meshgrid(grid_x, grid_z, indexing = 'ij')
    return func(vx, vz)

def test_value(op, res_expr, sol, grid_x, grid_z):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    mesh = xz_to_phys(res_expr, grid_x, grid_z)
    lhs = transf.tocheb2d(mesh)
    lhs = lhs.reshape(nx*nz, order = 'F')
    rhs = op*lhs
    err = np.abs(rhs - sol)
    print("\t\tValue error: " + str(err))

def test_forward(op, res_expr, sol_expr, grid_x, grid_z, qx, qz):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    mesh = xz_to_phys(res_expr, grid_x, grid_z)
    lhs = transf.tocheb2d(mesh)
    lhs = lhs.reshape(nx*nz, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape(nx,nz, order='F')
    mesh = xz_to_phys(sol_expr,grid_x, grid_z)
    sol = transf.tocheb2d(mesh)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid_x, grid_z):
    """Perform a tau backward operation test"""

    print("\tBackward tau test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    rhs = xz_to_phys(res_expr, grid_x, grid_z)
    rhs = transf.tocheb2d(rhs)
    rhs = rhs.reshape(nx*nz, order = 'F')
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = lhs.reshape(nx,nz, order='F')
    sol = xz_to_phys(sol_expr, grid_x, grid_z)
    sol = transf.tocheb2d(sol)
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    vis_error(err, 'Tau backward error')
    vis_error(relerr, 'Tau backward relative error')
    print("\t\tMax tau backward error: " + str(np.max(err)))

def test_backward_galerkin(opA, opB, opS, res_expr, sol_expr, grid):
    """Perform a galerkin backward operation test"""

    print("\tBackward galkerin test")
    x = sy.Symbol('x')
    rhs = transf.tocheb(x_to_phys(res_expr,grid))
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.tocheb(x_to_phys(sol_expr,grid))
    err = np.abs(opS*lhs - sol)
    vis_error(err, 'Galerkin backward error')
    print("\t\tMax galerkin backward error: " + str(np.max(err)))

def d1(nx,nz, xg, zg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.d1(nx, nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    test_forward(A, sphys, ssol, xg, zg, 1, 0)

def e1(nx,nz, xg, zg):
    """Accuracy test for e1 operator"""

    print("e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.e1(nx, nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    test_forward(A, sphys, ssol, xg, zg, 0, 1)

def laplh(nx,nz, xg, zg):
    """Accuracy test for laplh operator"""

    print("laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.laplh(nx, nz, k, 0, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    test_forward(A, sphys, ssol, xg, zg, 2, 0)

def lapl2h(nx,nz, xg, zg):
    """Accuracy test for lapl2h operator"""

    print("lapl2h:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.lapl2h(nx, nz, k, 0, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) - 2*k**2*sy.diff(sphys,x,x) + k**4*sphys)
    test_forward(A, sphys, ssol, xg, zg, 2, 0)

def i2j1e1(nx,nz, xg, zg):
    """Accuracy test for i2j1e1 operator"""

    print("i2j1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j1e1(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 1)

def i2j2e2(nx,nz, xg, zg):
    """Accuracy test for i2j2e2 operator"""

    print("i2j2e2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2e2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

def i2(nx,nz, xg, zg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 0)

def i2j1(nx,nz, xg, zg):
    """Accuracy test for i2j1 operator"""

    print("i2j1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j1(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 1)

def i2j2(nx,nz, xg, zg):
    """Accuracy test for i2j2 operator"""

    print("i2j2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

def i2j2d2e2(nx,nz, xg, zg):
    """Accuracy test for i2j2d2e2 operator"""

    print("i2j2d2e2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2d2e2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sy.diff(sphys,x,x),z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 20, 20")
    k = np.random.ranf()*nx
    A = c2d.i2j2d2e2(nx,nz, {'x':{0:20}, 'z':{0:20}}).tocsr()
    B = c2d.i2j2(nx,nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(sy.diff(ssol,x,x),z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i2laplh(nx,nz, xg, zg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 0)

    print("\tbc = 20, 0")
    k = np.random.ranf()*nx
    A = c2d.i2laplh(nx, nz, k, {'x':{0:20}, 'z':{0:0}}).tocsr()
    B = c2d.i2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i2j1laplh(nx,nz, xg, zg):
    """Accuracy test for i2j1laplh operator"""

    print("i2j1laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j1laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 1)

def i2j2laplh(nx,nz, xg, zg):
    """Accuracy test for i2j2laplh operator"""

    print("i2j2laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j2laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

def i2j2lapl(nx,nz, xg, zg):
    """Accuracy test for i2j2lapl operator"""

    print("i2j2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j2lapl(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 20, 20")
    k = np.random.ranf()*nx
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}, 'priority':'x'}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = 1e3*(1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

    print("\tbc = 20, 20")
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}, 'priority':'z'}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    test_backward_tau(A, B, sphys, ssol, xg, zg)

    print("\tbc = 20, 20")
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}, 'priority':'sx'}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    test_backward_tau(A, B, sphys, ssol, xg, zg)

    print("\tbc = 20, 20")
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}, 'priority':'sz'}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i2j2lapl2D(nx,nz, xg, zg):
    """Accuracy test for i2j2lapl operator for 2D (k=0)"""

    print("i2j2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2lapl(nx, nz, 0, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) + sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tforward backward solve")
    bc = c2d.c2dbc.no_bc()
    s = 2
    bc['x']['rt'] = s
    bc['x']['cr'] = s
    bc['z']['rt'] = s
    bc['z']['cr'] = s
    A = c2d.i2j2(nx+s, nz+s, bc).tocsr()
    B = c2d.i2j2lapl(nx+s, nz+s, 0, bc).tocsr()
    ssol = sy.expand(sy.diff(sphys,x,x) + sy.diff(sphys,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i4j1(nx,nz, xg, zg):
    """Accuracy test for i4j1 operator"""

    print("i4j1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j1(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 1)

def i4j2(nx,nz, xg, zg):
    """Accuracy test for i4j2 operator"""

    print("i4j2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 2)

def i4j4(nx,nz, xg, zg):
    """Accuracy test for i4j4 operator"""

    print("i4j4:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j4(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 4)

def i4j1e1(nx,nz, xg, zg):
    """Accuracy test for i4j1e1 operator"""

    print("i4j1e1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j1e1(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 2)

def i4j2e2(nx,nz, xg, zg):
    """Accuracy test for i4j2e2 operator"""

    print("i4j2e2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j2e2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 2)

def i4j1laplh(nx,nz, xg, zg):
    """Accuracy test for i4j1laplh operator"""

    print("i4j1laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j1laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 1)

def i4j2laplh(nx,nz, xg, zg):
    """Accuracy test for i4j2laplh operator"""

    print("i4j2laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j2laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 2)

def i4j4d1(nx,nz, xg, zg):
    """Accuracy test for i4j2d1 operator"""

    print("i4j4d1:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i4j4d1(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, xg, zg, 4, 4)

def i4j4lapl(nx,nz, xg, zg):
    """Accuracy test for i4j4lapl operator"""

    print("i4j4lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 4)

def i4j1lapl2h(nx,nz, xg, zg):
    """Accuracy test for i4j1lapl2h operator"""

    print("i4j1lapl2h:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j1lapl2h(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 1)

def i4j2lapl2h(nx,nz, xg, zg):
    """Accuracy test for i4j2lapl2h operator"""

    print("i4j2lapl2h:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j2lapl2h(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 2)

def i4j4lapl2(nx,nz, xg, zg):
    """Accuracy test for i4j4lapl2 operator"""

    print("i4j4lapl2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    ssol = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 4, 4)

    print("\tbc = 40, 40 (1)")
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}).tocsr()
    B = c2d.i4j4(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    sphys = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

    print("\tbc = 40, 40 (2)")
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}).tocsr()
    B = c2d.i4j4lapl(nx, nz, k, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i4j1_lapl2he1_e1laplh(nx,nz, xg, zg):
    """Accuracy test for coupled system lapl2he1_e1laplh operator"""

    print("i4j1_lapl2he1_e1laplh:")
    print("\tbc = (40, 00), (20,20) (uncoupled)")
    import scipy.io as io
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A11 = c2d.i4j1lapl2h(nx, nz, k, {'x':{0:40}, 'z':{0:00}})
    A12 = c2d.i4j1e1(nx, nz, {'x':{0:0}, 'z':{0:10, 'c':1.0}})
    A21 = c2d.i2j1e1(nx, nz, {'x':{0:0}, 'z':{0:0}}, -1.0)
    A22 = c2d.i2j1laplh(nx, nz, k, {'x':{0:20}, 'z':{0:11}})
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B11 = c2d.i4j1(nx, nz, c2d.c2dbc.no_bc())
    B12 = c2d.zblk(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B21 = c2d.zblk(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B22 = c2d.i2j1(nx, nz, c2d.c2dbc.no_bc())
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    io.mmwrite('matrix_A.mtx', A)
    io.mmwrite('matrix_B.mtx', B)

    ssol_psi = (1.0 - x**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz,1)])
    ssol_w = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys_psi = sy.expand(sy.diff(ssol_psi,x,x) - k**2*ssol_psi)
    sphys_psi = sy.expand(sy.diff(sphys_psi,x,x) - k**2*sphys_psi) + sy.expand(sy.diff(ssol_w,z))
    sphys_w = sy.expand(sy.diff(ssol_w,x,x) - k**2*ssol_w) - sy.expand(sy.diff(ssol_psi,z))

    rhs_psi = xz_to_phys(sphys_psi, xg, zg)
    rhs_w = xz_to_phys(sphys_w, xg, zg)
    rhs_psi = transf.tocheb2d(rhs_psi)
    rhs_w = transf.tocheb2d(rhs_w)
    rhs_psi = rhs_psi.reshape(nx*nz, order = 'F')
    rhs_w = rhs_w.reshape(nx*nz, order = 'F')
    rhs = np.append(rhs_psi, rhs_w)
    lhs = spsplin.spsolve(A,B*rhs)
    lhs_psi = lhs[0:nx*nz]
    lhs_w = lhs[nx*nz:]
    lhs_psi = lhs_psi.reshape(nx,nz, order='F')
    lhs_w = lhs_w.reshape(nx,nz, order='F')
    sol_psi = xz_to_phys(ssol_psi, xg, zg)
    sol_w = xz_to_phys(ssol_w, xg, zg)
    sol_psi = transf.tocheb2d(sol_psi)
    sol_w = transf.tocheb2d(sol_w)
    err_psi = np.abs(lhs_psi - sol_psi)
    err_w = np.abs(lhs_w - sol_w)
    pl.subplot(1,2,1)
    pl.imshow(np.log10(np.abs(err_psi)))
    pl.title('Tau backward error Psi')
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(np.log10(np.abs(err_w)))
    pl.title('Tau backward error W')
    pl.colorbar()
    pl.show()

    print("i4j1_lapl2he1_e1laplh:")
    print("\tbc = (40, 20), (20,20) (coupled)")
    import scipy.io as io
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A11 = c2d.i4j1lapl2h(nx, nz, k, {'x':{0:40}, 'z':{0:10}})
    A12 = c2d.i4j1e1(nx, nz, {'x':{0:0}, 'z':{0:10, 'c':-1.0}})
    A21 = c2d.i2j1e1(nx, nz, {'x':{0:0}, 'z':{0:11}}, -1.0)
    A22 = c2d.i2j1laplh(nx, nz, k, {'x':{0:20}, 'z':{0:11}})
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B11 = c2d.i4j1(nx, nz, c2d.c2dbc.no_bc())
    B12 = c2d.zblk(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B21 = c2d.zblk(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B22 = c2d.i2j1(nx, nz, c2d.c2dbc.no_bc())
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    io.mmwrite('matrix_A.mtx', A)
    io.mmwrite('matrix_B.mtx', B)

    ssol_psi = (1.0 - x**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-1,1)])
    ssol_w = ssol_psi*z
    sphys_psi = sy.expand(sy.diff(ssol_psi,x,x) - k**2*ssol_psi)
    sphys_psi = sy.expand(sy.diff(sphys_psi,x,x) - k**2*sphys_psi) + sy.expand(sy.diff(ssol_w,z))
    sphys_w = sy.expand(sy.diff(ssol_w,x,x) - k**2*ssol_w) - sy.expand(sy.diff(ssol_psi,z))

    rhs_psi = xz_to_phys(sphys_psi, xg, zg)
    rhs_w = xz_to_phys(sphys_w, xg, zg)
    rhs_psi = transf.tocheb2d(rhs_psi)
    rhs_w = transf.tocheb2d(rhs_w)
    rhs_psi = rhs_psi.reshape(nx*nz, order = 'F')
    rhs_w = rhs_w.reshape(nx*nz, order = 'F')
    rhs = np.append(rhs_psi, rhs_w)
    lhs = spsplin.spsolve(A,B*rhs)
    lhs_psi = lhs[0:nx*nz]
    lhs_w = lhs[nx*nz:]
    lhs_psi = lhs_psi.reshape(nx,nz, order='F')
    lhs_w = lhs_w.reshape(nx,nz, order='F')
    sol_psi = xz_to_phys(ssol_psi, xg, zg)
    sol_w = xz_to_phys(ssol_w, xg, zg)
    sol_psi = transf.tocheb2d(sol_psi)
    sol_w = transf.tocheb2d(sol_w)
    err_psi = np.abs(lhs_psi - sol_psi)
    err_w = np.abs(lhs_w - sol_w)
    pl.subplot(1,2,1)
    pl.imshow(np.log10(np.abs(err_psi)))
    pl.title('Tau backward error Psi')
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(np.log10(np.abs(err_w)))
    pl.title('Tau backward error W')
    pl.colorbar()
    pl.show()

def lapl2he1_e1laplh(nx,nz, xg, zg):
    """Accuracy test for coupled system lapl2he1_e1laplh operator"""

    print("lapl2he1_e1laplh:")
    print("\tbc = (40, 00), (20,20) (uncoupled)")
    import scipy.io as io
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A11 = c2d.lapl2h(nx, nz, k, 1, {'x':{0:40}, 'z':{0:0}})
    A12 = c2d.e1(nx, nz, 4, {'x':{0:0}, 'z':{0:10, 'c':-1.0}})
    A21 = c2d.e1(nx, nz, 2, {'x':{0:0}, 'z':{0:0}}, -1.0)
    A22 = c2d.laplh(nx, nz, k, 1, {'x':{0:20}, 'z':{0:11}})
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B11 = c2d.sid(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B12 = c2d.zblk(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B21 = c2d.zblk(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B22 = c2d.sid(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    io.mmwrite('matrix_A.mtx', A)
    io.mmwrite('matrix_B.mtx', B)

    ssol_psi = (1.0 - x**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz,1)])
    ssol_w = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys_psi = sy.expand(sy.diff(ssol_psi,x,x) - k**2*ssol_psi)
    sphys_psi = sy.expand(sy.diff(sphys_psi,x,x) - k**2*sphys_psi) + sy.expand(sy.diff(ssol_w,z))
    sphys_w = sy.expand(sy.diff(ssol_w,x,x) - k**2*ssol_w) - sy.expand(sy.diff(ssol_psi,z))

    rhs_psi = xz_to_phys(sphys_psi, xg, zg)
    rhs_w = xz_to_phys(sphys_w, xg, zg)
    rhs_psi = transf.tocheb2d(rhs_psi)
    rhs_w = transf.tocheb2d(rhs_w)
    rhs_psi = rhs_psi.reshape(nx*nz, order = 'F')
    rhs_w = rhs_w.reshape(nx*nz, order = 'F')
    rhs = np.append(rhs_psi, rhs_w)
    lhs = spsplin.spsolve(A,B*rhs)
    lhs_psi = lhs[0:nx*nz]
    lhs_w = lhs[nx*nz:]
    lhs_psi = lhs_psi.reshape(nx,nz, order='F')
    lhs_w = lhs_w.reshape(nx,nz, order='F')
    sol_psi = xz_to_phys(ssol_psi, xg, zg)
    sol_w = xz_to_phys(ssol_w, xg, zg)
    sol_psi = transf.tocheb2d(sol_psi)
    sol_w = transf.tocheb2d(sol_w)
    err_psi = np.abs(lhs_psi - sol_psi)
    err_w = np.abs(lhs_w - sol_w)
    pl.subplot(1,2,1)
    pl.imshow(np.log10(np.abs(err_psi)))
    pl.title('Tau backward error Psi')
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(np.log10(np.abs(err_w)))
    pl.title('Tau backward error W')
    pl.colorbar()
    pl.show()

    print("lapl2he1_e1laplh:")
    print("\tbc = (40, 20), (20,20) (coupled)")
    import scipy.io as io
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A11 = c2d.lapl2h(nx, nz, k, 1, {'x':{0:40}, 'z':{0:10}})
    A12 = c2d.e1(nx, nz, 4, {'x':{0:0}, 'z':{0:10, 'c':-1.0}})
    A21 = c2d.e1(nx, nz, 2, {'x':{0:0}, 'z':{0:11}}, -1.0)
    A22 = c2d.laplh(nx, nz, k, 1, {'x':{0:20}, 'z':{0:11}})
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B11 = c2d.sid(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B12 = c2d.zblk(nx, nz, 4, 1, c2d.c2dbc.no_bc())
    B21 = c2d.zblk(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B22 = c2d.sid(nx, nz, 2, 1, c2d.c2dbc.no_bc())
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    io.mmwrite('matrix_A.mtx', A)
    io.mmwrite('matrix_B.mtx', B)

    ssol_psi = (1.0 - x**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-1,1)])
    ssol_w = ssol_psi*z
    sphys_psi = sy.expand(sy.diff(ssol_psi,x,x) - k**2*ssol_psi)
    sphys_psi = sy.expand(sy.diff(sphys_psi,x,x) - k**2*sphys_psi) + sy.expand(sy.diff(ssol_w,z))
    sphys_w = sy.expand(sy.diff(ssol_w,x,x) - k**2*ssol_w) - sy.expand(sy.diff(ssol_psi,z))

    rhs_psi = xz_to_phys(sphys_psi, xg, zg)
    rhs_w = xz_to_phys(sphys_w, xg, zg)
    rhs_psi = transf.tocheb2d(rhs_psi)
    rhs_w = transf.tocheb2d(rhs_w)
    rhs_psi = rhs_psi.reshape(nx*nz, order = 'F')
    rhs_w = rhs_w.reshape(nx*nz, order = 'F')
    rhs = np.append(rhs_psi, rhs_w)
    lhs = spsplin.spsolve(A,B*rhs)
    lhs_psi = lhs[0:nx*nz]
    lhs_w = lhs[nx*nz:]
    lhs_psi = lhs_psi.reshape(nx,nz, order='F')
    lhs_w = lhs_w.reshape(nx,nz, order='F')
    sol_psi = xz_to_phys(ssol_psi, xg, zg)
    sol_w = xz_to_phys(ssol_w, xg, zg)
    sol_psi = transf.tocheb2d(sol_psi)
    sol_w = transf.tocheb2d(sol_w)
    err_psi = np.abs(lhs_psi - sol_psi)
    err_w = np.abs(lhs_w - sol_w)
    pl.subplot(1,2,1)
    pl.imshow(np.log10(np.abs(err_psi)))
    pl.title('Tau backward error Psi')
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(np.log10(np.abs(err_w)))
    pl.title('Tau backward error W')
    pl.colorbar()
    pl.show()

def surfaceAvg(nx, nz, xg, zg):
    """Accuracy test for the surface average"""

    print("surfaceAvg:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.surfaceAvg(nx, nz)
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.integrate(sy.expand(sphys/4),(x,-1,1),(z,-1,1))
    test_value(A, sphys, ssol, xg, zg)


if __name__ == "__main__":
    # Set test parameters
    nx = 12
    nz = 12
    xg = transf.grid(nx)
    zg = transf.grid(nz)

    # run hardcoded operator tests
    print('Hard coded exact operators')
#    d1(nx, nz, xg, zg)
#    e1(nx, nz, xg, zg)
#    laplh(nx, nz, xg, zg)
#    lapl2h(nx, nz, xg, zg)
#    i2j1e1(nx, nz, xg, zg)
#    i2j2e2(nx, nz, xg, zg)
#    i2j2d2e2(nx, nz, xg, zg)
#    i2(nx, nz, xg, zg)
#    i2j1(nx, nz, xg, zg)
#    i2j2(nx, nz, xg, zg)
#    i2laplh(nx, nz, xg, zg)
#    i2j1laplh(nx, nz, xg, zg)
#    i2j2laplh(nx, nz, xg, zg)
#    i2j2lapl(nx, nz, xg, zg)
    i2j2lapl2D(nx, nz, xg, zg)
#    i4j1(nx, nz, xg, zg)
#    i4j2(nx, nz, xg, zg)
#    i4j4(nx, nz, xg, zg)
#    i4j1e1(nx, nz, xg, zg)
#    i4j2e2(nx, nz, xg, zg)
#    i4j1laplh(nx, nz, xg, zg)
#    i4j2laplh(nx, nz, xg, zg)
#    i4j4d1(nx, nz, xg, zg)
#    i4j4lapl(nx, nz, xg, zg)
#    i4j1lapl2h(nx, nz, xg, zg)
#    i4j2lapl2h(nx, nz, xg, zg)
#    i4j4lapl2(nx, nz, xg, zg)
#    lapl2he1_e1laplh(nx, nz, xg, zg)
#    i4j1_lapl2he1_e1laplh(nx, nz, xg, zg)

    # Run average operator tests
#    surfaceAvg(nx, nz, xg, zg)
