"""Check accuracy for cartesian 3D operators"""

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
import quicc.geometry.cartesian.cartesian_3d as c3d


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        if has_error_plot:
            pl.imshow(np.log10(np.abs(err.reshape((err.shape[0], err.shape[1]*err.shape[2])))))
            pl.title(title)
            pl.colorbar()
            pl.show()

def xyz_to_phys(expr, grid_x, grid_z, grid_y):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z, y), expr)
    vx, vz, vy = np.meshgrid(grid_x, grid_z, grid_y, indexing = 'ij')
    return func(vx, vz, vy)

def test_value(op, res_expr, sol, grid_x, grid_y, grid_z):
    """Perform a forward operation test"""

    print("\tValue test")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    nx = len(grid_x)
    ny = len(grid_y)
    nz = len(grid_x)
    mesh = xyz_to_phys(res_expr, grid_x, grid_z, grid_y)
    lhs = transf.tocheb3d(mesh)
    lhs = lhs.reshape(nx*ny*nz, order = 'F')
    rhs = op*lhs
    print(rhs)
    print(sol)
    err = np.abs(rhs - sol)
    print("\t\tValue error: " + str(err))

def test_forward(op, res_expr, sol_expr, grid_x, grid_y, grid_z, qx, qy, qz):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    nx = len(grid_x)
    ny = len(grid_y)
    nz = len(grid_x)
    mesh = xyz_to_phys(res_expr, grid_x, grid_z, grid_y)
    lhs = transf.tocheb3d(mesh)
    lhs = lhs.reshape(nx*ny*nz, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape((nx, nz, ny), order='F')
    mesh = xyz_to_phys(sol_expr, grid_x, grid_z, grid_y)
    sol = transf.tocheb3d(mesh)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid_x, grid_y, grid_z):
    """Perform a tau backward operation test"""

    print("\tBackward tau test")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    nx = len(grid_x)
    ny = len(grid_y)
    nz = len(grid_x)
    rhs = xyz_to_phys(res_expr, grid_x, grid_z, grid_y)
    rhs = transf.tocheb3d(rhs)
    rhs = rhs.reshape(nx*ny*nz, order='F')
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = lhs.reshape(nx,nz,ny, order = 'F')
    sol = xyz_to_phys(sol_expr, grid_x, grid_z, grid_y)
    sol = transf.tocheb3d(sol)
    err = np.abs(lhs - sol)
    vis_error(err, 'Tau backward error')
    print("\t\tMax tau backward error: " + str(np.max(err)))

def d1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.d1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    test_forward(A, sphys, ssol, xg, yg, zg, 1, 0, 0)

def e1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for e1 operator"""

    print("e1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.e1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,y))
    test_forward(A, sphys, ssol, xg, yg, zg, 0, 1, 0)

def f1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for f1 operator"""

    print("f1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.f1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    test_forward(A, sphys, ssol, xg, yg, zg, 0, 0, 1)

def i1j1k1d1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i1j1k1d1 operator"""

    print("i1j1k1d1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i1j1k1d1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i1j1k1e1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i1j1k1e1 operator"""

    print("i1j1k1e1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i1j1k1e1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,y))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i1j1k1f1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i1j1k1f1 operator"""

    print("i1j1k1f1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i1j1k1f1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2d2e2f2(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2d2e2f2 operator"""

    print("i2j2k2d2e2f2:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2d2e2f2(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sy.diff(sy.diff(sphys,x,x),y,y),z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2 operator"""

    print("i2j2k2:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.integrate(sphys,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2d1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2d1 operator"""

    print("i2j2k2d1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2d1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2e1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2e1 operator"""

    print("i2j2k2e1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2e1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,y))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2f1(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2f1 operator"""

    print("i2j2k2f1:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2f1(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

def i2j2k2lapl(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2lapl operator"""

#    print("i2j2k2lapl:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i2j2k2lapl(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 2, 2, 2)

    print("\tbc = 20, 20 ,20")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'xy'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - y**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 20, 20, 21")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:20}, 'z':{0:21}, 'priority':'xy'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - y**2)*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 20, 21, 20")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:21}, 'z':{0:20}, 'priority':'xz'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - y**2)**2*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 21, 20, 20")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:21}, 'y':{0:20}, 'z':{0:20}, 'priority':'zy'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - y**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 21, 21, 20")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:21}, 'y':{0:21}, 'z':{0:20}, 'priority':'zsy'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - y**2)**2*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 21, 20, 21")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:21}, 'y':{0:20}, 'z':{0:21}, 'priority':'ysx'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - y**2)*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 20, 21, 21")
    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:21}, 'z':{0:21}, 'priority':'xsz'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - y**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

#    print("\tbc = 21, 21, 21")
#    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:21}, 'y':{0:21}, 'z':{0:21}, 'priority':'sxsz'}).tocsr()
#    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**2*(1.0 - y**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-4,1)])
#    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)


def i4j4k4(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i4j4k4 operator"""

    print("i4j4k4:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 4, 4, 4)

def i4j4k4lapl(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i4j4k4lapl operator"""

    print("i4j4k4lapl:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.i4j4k4lapl(nx, ny, nz, c3d.c3dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,y,y,y,y)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    test_forward(A, sphys, ssol, xg, yg, zg, 4, 4, 4)

def i4j4k4lapl2(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i4j4k4lapl2 operator"""

    print("i4j4k4lapl2:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
#    A = c3d.i4j4k4lapl2(nx, ny, nz, c3d.c3dbc.no_bc())
#    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
#    ssol = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    ssol = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    ssol = sy.integrate(ssol,x,x,x,x)
#    ssol = sy.expand(ssol)
#    ssol = sy.integrate(ssol,y,y,y,y)
#    ssol = sy.expand(ssol)
#    ssol = sy.integrate(ssol,z,z,z,z)
#    test_forward(A, sphys, ssol, xg, yg, zg, 4, 4, 4)
#
#    print("\tbc = 40, 40 ,40")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:40}, 'y':{0:40}, 'z':{0:40}, 'priority':'xy'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**2*(1.0 - y**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-4,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)
#
#    print("\tbc = 40, 40 ,41")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:40}, 'y':{0:40}, 'z':{0:41}, 'priority':'xy'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**2*(1.0 - y**2)**2*(1.0 - z**2)**3*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-6,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)
#
#    print("\tbc = 40, 41 ,40")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:40}, 'y':{0:41}, 'z':{0:40}, 'priority':'xz'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**2*(1.0 - y**2)**3*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-6,1)]) for k in np.arange(0,nz-4,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)
#
#    print("\tbc = 41, 40 ,40")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:41}, 'y':{0:40}, 'z':{0:40}, 'priority':'yz'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**3*(1.0 - y**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-6,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-4,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

    print("\tbc = 41, 41 ,40")
    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:41}, 'y':{0:41}, 'z':{0:40}, 'priority':'zsx'}).tocsr()
    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**3*(1.0 - y**2)**3*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-6,1)]) for j in np.arange(0,ny-6,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)
#
#    print("\tbc = 41, 40 ,41")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:41}, 'y':{0:40}, 'z':{0:41}, 'priority':'ysz'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**3*(1.0 - y**2)**2*(1.0 - z**2)**3*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-6,1)]) for j in np.arange(0,ny-4,1)]) for k in np.arange(0,nz-6,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)
#
#    print("\tbc = 40, 41 ,41")
#    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:40}, 'y':{0:41}, 'z':{0:41}, 'priority':'xsy'}).tocsr()
#    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
#    ssol = (1.0 - x**2)**2*(1.0 - y**2)**3*(1.0 - z**2)**3*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,ny-6,1)]) for k in np.arange(0,nz-6,1)])
#    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
#    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
#    test_backward_tau(A, B, sphys, ssol, xg, yg, zg)

def volumeAvg(nx, ny, nz, xg, yg, zg):
    """Accuracy test for volume average"""

    print("volumeAvg:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.volumeAvg(nx, ny, nz)
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    ssol = sy.integrate(sy.expand(sphys/8),(x,-1,1),(y,-1,1), (z,-1,1))
    test_value(A, sphys, ssol, xg, yg, zg)

def avgFlux_z(nx, ny, nz, xg, yg, zg):
    """Accuracy test for volume average"""

    print("avgFlux_z:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    A = c3d.avgFlux_z(nx, ny, nz, zscale = 2.0)
    sphys = np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,ny,1)]) for k in np.arange(0,nz,1)])
    tmp = sy.integrate(sy.expand(sphys/4),(x,-1,1),(y,-1,1))
    func = sy.utilities.lambdify(z, sy.expand(sy.diff(tmp,z)))
    ssol = func(1.0)
    test_value(A, sphys, ssol, xg, yg, zg)

if __name__ == "__main__":
    # Set test parameters
    nx = 8
    ny = 8
    nz = 8
    xg = transf.grid(nx)
    yg = transf.grid(nx)
    zg = transf.grid(nz)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    d1(nx, ny, nz, xg, yg, zg)
    e1(nx, ny, nz, xg, yg, zg)
    f1(nx, ny, nz, xg, yg, zg)
    i1j1k1d1(nx, ny, nz, xg, yg, zg)
    i1j1k1e1(nx, ny, nz, xg, yg, zg)
    i1j1k1f1(nx, ny, nz, xg, yg, zg)
    i2j2k2d2e2f2(nx, ny, nz, xg, yg, zg)
    i2j2k2(nx, ny, nz, xg, yg, zg)
    i2j2k2d1(nx, ny, nz, xg, yg, zg)
    i2j2k2e1(nx, ny, nz, xg, yg, zg)
    i2j2k2f1(nx, ny, nz, xg, yg, zg)
    i2j2k2lapl(nx, ny, nz, xg, yg, zg)
    i4j4k4(nx, ny, nz, xg, yg, zg)
    i4j4k4lapl(nx, ny, nz, xg, yg, zg)
    i4j4k4lapl2(nx, ny, nz, xg, yg, zg)

    # Run average operator tests
    volumeAvg(nx, ny, nz, xg, yg, zg)
    avgFlux_z(nx, ny, nz, xg, yg, zg)
