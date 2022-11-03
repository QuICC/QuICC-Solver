"""Generate test problem matrices for sparse solver tests"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.io as io
import scipy.sparse as spsp
if True:
    import matplotlib.pylab as pl
    has_error_plot = True
else:
    has_error_plot = False

import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_2d as c2d
import quicc.geometry.cartesian.cartesian_3d as c3d
import quicc.transform.cartesian as transf
import quicc.base.utils as utils

save_path = ""


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

def xz_to_phys(expr, grid_x, grid_z):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z), expr)
    vx, vz = np.meshgrid(grid_x, grid_z, indexing = 'ij')
    return func(vx, vz)

def xyz_to_phys(expr, grid_x, grid_z, grid_y):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z, y), expr)
    vx, vz, vy = np.meshgrid(grid_x, grid_z, grid_y, indexing = 'ij')
    return func(vx, vz, vy)

def build1DSol(opB, res_expr, sol_expr, grid):
    """Compute spectral expansion of solution"""

    x = sy.Symbol('x')
    nx = len(grid)
    rhs = transf.tocheb(x_to_phys(res_expr,grid))
    rhs = opB*rhs
    rhs = rhs.reshape(nx, 1, order = 'F')
    sol = transf.tocheb(x_to_phys(sol_expr,grid))
    sol = sol.reshape(nx, 1, order = 'F')

    return (spsp.coo_matrix(rhs), spsp.coo_matrix(sol))

def build2DSol(opB, res_expr, sol_expr, grid_x, grid_z):
    """Compute spectral expansion of solution"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    rhs = xz_to_phys(res_expr, grid_x, grid_z)
    rhs = transf.tocheb2d(rhs)
    rhs = rhs.reshape(nx*nz, order = 'F')
    rhs = opB*rhs
    rhs = rhs.reshape(nx*nz, 1, order = 'F')
    sol = xz_to_phys(sol_expr, grid_x, grid_z)
    sol = transf.tocheb2d(sol)
    sol = sol.reshape(nx*nz, 1, order = 'F')

    return (spsp.coo_matrix(rhs), spsp.coo_matrix(sol))

def build3DSol(opB, res_expr, sol_expr, grid_x, grid_y, grid_z):
    """Compute spectral expansion of solution"""

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    nx = len(grid_x)
    ny = len(grid_y)
    nz = len(grid_x)
    rhs = xyz_to_phys(res_expr, grid_x, grid_z, grid_y)
    rhs = transf.tocheb3d(rhs)
    rhs = rhs.reshape(nx*ny*nz, order='F')
    rhs = opB*rhs
    rhs = rhs.reshape(nx*ny*nz, 1, order = 'F')
    sol = xyz_to_phys(sol_expr, grid_x, grid_z, grid_y)
    sol = transf.tocheb3d(sol)
    sol = sol.reshape(nx*ny*nz, 1, order = 'F')

    return (spsp.coo_matrix(rhs), spsp.coo_matrix(sol))

def writeTest(opA, rhs, sol, name):
    """Perform a tau backward operation test"""

    nn = opA.shape[0]
    io.mmwrite(save_path + name + '_A_' +  str(nn) + '.mtx', opA)
    io.mmwrite(save_path + name + '_rhs_' +  str(nn) + '.mtx', rhs)
    io.mmwrite(save_path + name + '_sol_' +  str(nn) + '.mtx', sol)

def laplacian1D(nx, ny, nz, restriction = None):
    """Accuracy test for i2lapl operator"""

    k = nx/2
    l = nx/3
    A = c1d.i2lapl(nx, k, l, {0:20}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()

    return (A, B)

def bilaplacian1D(nx, ny, nz, restriction = None):
    """Accuracy test for i4lapl2 operator"""

    k = nx/2
    l = nx/3
    A = c1d.i4lapl2(nx, k, l, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()

    return (A, B)

def laplacian2D(nx, ny, nz, restriction = None):
    """Accuracy test for i2j2lapl operator"""

    k = nx/3
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}}, restriction = restriction).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc(), restriction = restriction).tocsr()

    return (A, B)

def bilaplacian2D(nx, ny, nz, restriction = None):
    """Accuracy test for i4j4lapl2 operator"""

    k = nx/3
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}, restriction = restriction).tocsr()
    B = c2d.i4j4(nx, nz, c2d.c2dbc.no_bc(), restriction = restriction).tocsr()

    return (A, B)

def laplacian3D(nx, ny, nz, restriction = None):
    """Accuracy test for i2j2k2lapl operator"""

    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'xy'}, restriction = restriction).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc(), restriction = restriction).tocsr()

    return (A, B)

def bilaplacian3D(nx, ny, nz, restriction = None):
    """Accuracy test for i4j4k4lapl2 operator"""

    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:41}, 'y':{0:41}, 'z':{0:40}, 'priority':'zsx'}, restriction = restriction).tocsr()
    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc(), restriction = restriction).tocsr()

    return (A, B)

def laplacian1D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i2lapl operator"""

    A, B = laplacian1D(nx, ny, nx, restriction = restriction)

    xg = transf.grid(nx)

    x = sy.Symbol('x')
    k = nx/2
    l = nx/3
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    rhs, sol = build1DSol(B, sphys, ssol, xg)

    return (A, rhs, sol)

def bilaplacian1D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i4lapl2 operator"""

    A, B = bilaplacian1D(nx, ny, nz, restriction = restriction)

    xg = transf.grid(nx)

    x = sy.Symbol('x')
    k = nx/2
    l = nx/3
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + k**4*ssol + l**4*ssol - 2*k**2*sy.diff(ssol,x,x) - 2*l**2*sy.diff(ssol,x,x) + 2*k**2*l**2*ssol)
    rhs, sol = build1DSol(B, sphys, ssol, xg)

    return (A, rhs, sol)

def laplacian2D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i2j2lapl operator"""

    A, B = laplacian2D(nx, ny, nz, restriction = restriction)

    xg = transf.grid(nx)
    zg = transf.grid(nz)

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = nx/3
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    rhs, sol = build2DSol(B, sphys, ssol, xg, zg)

    return (A, rhs, sol)

def bilaplacian2D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i4j4lapl2 operator"""

    A, B = bilaplacian2D(nx, ny, nz, restriction = restriction)

    xg = transf.grid(nx)
    zg = transf.grid(nz)

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = nx/3
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    sphys = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    rhs, sol = build2DSol(B, sphys, ssol, xg, zg)

    return (A, rhs, sol)

def laplacian3D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i2j2k2lapl operator"""

    A, B = laplacian3D(nx, ny, nz, restriction = None)

    xg = transf.grid(nx)
    yg = transf.grid(ny)
    zg = transf.grid(nz)

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')

    np.random.seed(42)

    ssol = (1.0 - x**2)*(1.0 - y**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    rhs, sol = build3DSol(B, sphys, ssol, xg, yg, zg)

    A, B = laplacian3D(nx, ny, nz, restriction = restriction)
    R = c3d.qid(nx, ny, nz, 0, 0, 0, {'x':{0:0}, 'y':{0:0}, 'z':{0:0}}, restriction = restriction)
    rhs = R*rhs
#    sol = R*sol

    return (A, rhs, sol)

def bilaplacian3D_exact(nx, ny, nz, restriction = None):
    """Accuracy test for i4j4k4lapl2 operator"""

    A, B = bilaplacian3D(nx, ny, nz, restriction = None)

    xg = transf.grid(nx)
    yg = transf.grid(ny)
    zg = transf.grid(nz)

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')

    np.random.seed(42)

    ssol = (1.0 - x**2)**3*(1.0 - y**2)**3*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-6,1)]) for j in np.arange(0,ny-6,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
    rhs, sol = build3DSol(B, sphys, ssol, xg, yg, zg)
    A, B = bilaplacian3D(nx, ny, nz, restriction = restriction)

    R = c3d.qid(nx, ny, nz, 0, 0, 0, {'x':{0:0}, 'y':{0:0}, 'z':{0:0}}, restriction = restriction)
    rhs = R*rhs
#    sol = R*sol

    return (A, rhs, sol)

def rb1dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rb1dboxvc operator"""

    import quicc.model.boussinesq_rb1dbox_vc as mod
    model = mod.BoussinesqRB1DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, 0, 0]
    bc_vel = 1
    bc_temp = 0
    phi = 35
    kp = 2.221441469
    kx = kp*np.cos(phi*np.pi/180.0);
    ky = (kp**2-kx**2)**0.5;
    eq_params = {'prandtl':1, 'rayleigh':657.5113645, 'scale1d':1.0}
    eigs = [kx, ky]
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields)

    return (A, B)

def rrb1dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rrb1dboxvc operator"""

    import quicc.model.boussinesq_rrb1dbox_vc as mod
    model = mod.BoussinesqRRB1DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, 0, 0]
    bc_vel = 1
    bc_temp = 0
    kx = 1
    ky = 2.710
    eq_params = {'prandtl':1, 'rayleigh':1676.12, 'taylor':1e3, 'scale1d':1.0}
    eigs = [kx, ky]
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields)

    return (A, B)

def rb2dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rb2dboxvc operator"""

    import quicc.model.boussinesq_rb2dbox_vc as mod
    model = mod.BoussinesqRB2DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, 0, nz]
    bc_vel = 2
    bc_temp = 1
    eigs = [3]
    eq_params = {'prandtl':1, 'rayleigh':10823.2, 'scale1d':1.0, 'scale3d':1.0} # m = 3, n = 1, aspect ration 1:1
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields, restriction = restriction)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields, restriction = restriction)

    return (A, B)

def rrb2dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rrb2dboxvc operator"""

    import quicc.model.boussinesq_rrb2dbox_vc as mod
    model = mod.BoussinesqRRB2DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, 0, nz]
    bc_vel = 2
    bc_temp = 1
    eigs = [5]
    eq_params = {'prandtl':1, 'rayleigh':1e5, 'taylor':1e3, 'scale1d':1.0, 'scale3d':1.0} # m = 1, n = 1, aspect ration 1:1
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields, restriction = restriction)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields, restriction = restriction)

    return (A, B)

def rb3dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rb3dboxvc operator"""

    import quicc.model.boussinesq_rb3dbox_vc as mod
    model = mod.BoussinesqRB3DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, ny, nz]
    bc_vel = 6
    bc_temp = 4
    eq_params = {'prandtl':1, 'rayleigh':1315.022729, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1, m = 1, n = 1, aspect ration 1:1:1
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    eigs = []
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields, restriction = restriction)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields, restriction = restriction)

    return (A, B)

def rrb3dboxvc(nx, ny, nz, restriction = None):
    """Accuracy test for rrb3dboxvc operator"""

    import quicc.model.boussinesq_rrb3dbox_vc as mod
    model = mod.BoussinesqRRB3DBoxVC()
    model.linearize = True
    model.use_galerkin = False
    fields = model.stability_fields()
    res = [nx, ny, nz]
    bc_vel = 6
    bc_temp = 4
    eq_params = {'prandtl':1, 'rayleigh':779.2727283, 'taylor':1e3, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1|0, m = 0|1, n = 1, aspect ration 1:1:1
    bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}
    eigs = []
    A = model.implicit_linear(res, eq_params, eigs, bcs, fields, restriction = restriction)
    bcs['bcType'] = model.SOLVER_NO_TAU
    B = model.time(res, eq_params, eigs, bcs, fields, restriction = restriction)

    return (A, B)

if __name__ == "__main__":
    # Set test parameters
    ns = [16, 32, 64, 128, 256, 512]
    for nn in ns:
        A, rhs, sol = laplacian1D_exact(nn, 0, 0)
        writeTest(A, rhs, sol, 'laplacian1D')
        A, rhs, sol = bilaplacian1D_exact(nn, 0, 0)
        writeTest(A, rhs, sol, 'bilaplacian1D')

    ns = [8, 16, 24, 32, 48]
    for nn in ns:
        A, rhs, sol = laplacian2D_eaxt(nn, 0, nn)
        writeTest(A, rhs, sol, 'laplacian2D')
        A, rhs, sol = bilaplacian2D_exact(nn, 0, nn)
        writeTest(A, rhs, sol, 'bilaplacian2D')

    ns = [8, 10, 12, 14, 16]
    for nn in ns:
        A, rhs, sol = laplacian3D_exact(nn, nn, nn)
        writeTest(A, rhs, sol, 'laplacian3D')
        A, rhs, sol = bilaplacian3D_exact(nn, nn, nn)
        writeTest(A, rhs, sol, 'bilaplacian3D')
