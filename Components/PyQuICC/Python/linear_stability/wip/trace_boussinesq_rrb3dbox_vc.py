"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)"""

import numpy as np

import quicc.model.wip.boussinesq_rrb3dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRRB3DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [6, 6, 6]
#eq_params = {'prandtl':1, 'rayleigh':2340.687, 'zxratio':1.0, 'yxratio':1.0}
eq_params = {'prandtl':1, 'rayleigh':5011.73, 'taylor':1e4, 'zxratio':1.0, 'yxratio':1.0}
eigs = []
bc_vel = 0 # 0: NS/NS/NS, 1: SF/SF/SF, 2: SF/NS/NS, 3: SF/NS/NS
bc_temp = 0 # 0: FT/FT/FT, 1: FF/FF/FF, 2: FF/FT/FT, 3: FT/FF/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import quicc.transform.cartesian as transf

# Show the "spy" of the two matrices
if show_spy:
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

if show_solution:
    mode = -1

    # Get solution vectors
    sol_u = evp_vec[0:res[0]*res[1]*res[2],mode]
    sol_v = evp_vec[res[0]*res[1]*res[2]:2*res[0]*res[1]*res[2],mode]
    sol_w = evp_vec[2*res[0]*res[1]*res[2]:3*res[0]*res[1]*res[2],mode]
    sol_t = evp_vec[3*res[0]*res[1]*res[2]:4*res[0]*res[1]*res[2],mode]
    sol_p = evp_vec[4*res[0]*res[1]*res[2]:5*res[0]*res[1]*res[2],mode]
    # Extract continuity from velocity 
    sol_c = mod.c3d.d1(res[0], res[1], res[2], mod.no_bc(), sy = 0, sz = 0)*sol_u + mod.c3d.e1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sz = 0)*sol_v + mod.c3d.f1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sy = 0)*sol_w
    
    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_u))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.semilogy(np.abs(sol_v))
    pl.title("v")
    pl.subplot(2,3,3)
    pl.semilogy(np.abs(sol_w))
    pl.title("w")
    pl.subplot(2,3,4)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.subplot(2,3,5)
    pl.semilogy(np.abs(sol_c))
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.semilogy(np.abs(sol_p))
    pl.title("p")
    pl.show()
    pl.close("all")
