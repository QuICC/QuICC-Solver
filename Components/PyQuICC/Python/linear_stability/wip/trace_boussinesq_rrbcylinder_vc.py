"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a cylinder model (velocity-continuity formulation)"""

import numpy as np

import quicc.model.wip.boussinesq_rrbcylinder_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRRBCylinderVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [6, 0, 6]
eq_params = {'prandtl':1, 'rayleigh':2901.55, 'taylor':1e4, 'ratio31':1.0}
#eq_params = {'prandtl':1, 'rayleigh':0., 'ratio31':1.0}
eigs = [1]
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':bc_vel, 'velocity_y':bc_vel, 'velocity_z':bc_vel, 'temperature':bc_temp}

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
    import quicc.transform.cylinder as transf

# Show the "spy" of the two matrices
if show_spy:
    import matplotlib.pylab as pl
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
    sol_u = evp_vec[0:res[0]*res[2],mode]
    sol_v = evp_vec[res[0]*res[2]:2*res[0]*res[2],mode]
    sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],mode]
    sol_t = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],mode]
    sol_p = evp_vec[4*res[0]*res[2]:5*res[0]*res[2],mode]
    # Extract continuity from velocity 
    sol_c = mod.cylinder.x1div(res[0], res[2], (eigs[0]+1)%2, mod.no_bc(), sz = 0)*sol_u + 1j*eigs[0]*sol_v + mod.cylinder.x1e1(res[0], res[2], eigs[0]%2, mod.no_bc(), sr = 0)*sol_w
    
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
    
    # Create solution matrices
    mat_u = sol_u.reshape(res[0], res[2], order = 'F')
    mat_v = sol_v.reshape(res[0], res[2], order = 'F')
    mat_w = sol_w.reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_p = sol_p.reshape(res[0], res[2], order = 'F')
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')

    # Visualize spectrum matrix
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(mat_v)))
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(mat_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.imshow(np.log10(np.abs(mat_t)))
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(mat_c)))
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(mat_p)))
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    # Compute physical space values
    grid_r = transf.rgrid(res[0])
    grid_z = transf.zgrid(res[2])
    phys_u = transf.tophys2d(mat_u, (eigs[0]+1)%2)
    phys_v = transf.tophys2d(mat_v, (eigs[0]+1)%2)
    phys_w = transf.tophys2d(mat_w, eigs[0]%2)
    phys_t = transf.tophys2d(mat_t, eigs[0]%2)
    phys_p = transf.tophys2d(mat_p, eigs[0]%2)
    phys_c = transf.tophys2d(mat_c, (eigs[0]+1)%2)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_z, grid_r, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_z, grid_r, phys_v, 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.contourf(grid_z, grid_r, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_z, grid_r, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_z, grid_r, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_z, grid_r, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
