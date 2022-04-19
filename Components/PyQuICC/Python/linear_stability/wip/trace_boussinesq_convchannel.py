"""Script to run a marginal curve trace for the Boussinesq rotating convection in a periodic channel model (Velocity-continuity)"""

import numpy as np

import quicc.model.boussinesq_convchannel as mod

# Create the model and activate linearization
model = mod.BoussinesqConvChannel()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [100, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':1707.8}
phi = 0
kp = 3.117
kx = kp*np.cos(phi*np.pi/180.0);
ky = (kp**2-kx**2)**0.5;
eigs = [kx, ky]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp, 'pressure':bc_vel}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
if False:
    import matplotlib.pylab as pl
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if True:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if True:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

    sol_u = evp_vec[0:res[0],-1]
    sol_v = evp_vec[res[0]:2*res[0],-1]
    sol_w = evp_vec[2*res[0]:3*res[0],-1]
    sol_t = evp_vec[3*res[0]:4*res[0],-1]
    sol_p = evp_vec[4*res[0]:5*res[0],-1]
    sol_c = 1j*(kx/2.0)*np.array(sol_u) + 1j*(ky/2.0)*np.array(sol_v) + mod.c1d.d1(res[0], mod.no_bc()).dot(sol_w)

    import matplotlib.pylab as pl
    pl.subplot(2,3,1)
    pl.semilogy(abs(sol_u))
    pl.title('u')
    pl.subplot(2,3,2)
    pl.semilogy(abs(sol_u))
    pl.title('v')
    pl.subplot(2,3,3)
    pl.semilogy(abs(sol_u))
    pl.title('w')
    pl.subplot(2,3,4)
    pl.semilogy(abs(sol_t))
    pl.title('T')
    pl.subplot(2,3,5)
    pl.semilogy(abs(sol_p))
    pl.title('p')
    pl.subplot(2,3,6)
    pl.semilogy(abs(sol_c))
    pl.title('Continuity')
    pl.show()

    import quicc.transform.cartesian as transf
    phys_u = transf.tophys(sol_u.real)
    phys_v = transf.tophys(sol_v.real)
    phys_w = transf.tophys(sol_w.real)
    phys_t = transf.tophys(sol_t.real)
    phys_cr = transf.tophys(np.real(sol_c))
    phys_ci = transf.tophys(np.imag(sol_c))

    grid_x = transf.grid(res[0])
    
    pl.subplot(2,3,1)
    pl.plot(grid_x, phys_u)
    pl.title('u')
    pl.subplot(2,3,2)
    pl.plot(grid_x, phys_v)
    pl.title('v')
    pl.subplot(2,3,3)
    pl.plot(grid_x, phys_w)
    pl.title('w')
    pl.subplot(2,3,4)
    pl.plot(grid_x, phys_t)
    pl.title('T')
    pl.subplot(2,3,5)
    pl.plot(grid_x, phys_cr)
    pl.title('Continuity (real)')
    pl.subplot(2,3,6)
    pl.plot(grid_x, phys_ci)
    pl.title('Continuity (imag)')
    pl.show()

