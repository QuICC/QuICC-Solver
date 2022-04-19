"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a square (velocity-continuity formulation)"""

import numpy as np

import quicc.model.boussinesq_rbsquare_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRBSquareVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]

# SF/SF, FF/FT
bc_vel = 2
bc_temp = 1
# SF/SF, FF/FT, Aspect ratio 1:1
#eq_params = {'prandtl':1, 'rayleigh':779.273, 'scale1d':1.0, 'scale3d':1.0} # m = 1, n = 1, aspect ration 1:1
#eq_params = {'prandtl':1, 'rayleigh':3044.03, 'scale1d':1.0, 'scale3d':1.0} # m = 2, n = 1, aspect ration 1:1
#eq_params = {'prandtl':1, 'rayleigh':10823.2, 'scale1d':1.0, 'scale3d':1.0} # m = 3, n = 1, aspect ration 1:1
#eq_params = {'prandtl':1, 'rayleigh':12176.1, 'scale1d':1.0, 'scale3d':1.0} # m = 1, n = 2, aspect ration 1:1
eq_params = {'prandtl':1, 'rayleigh':12468.4, 'scale1d':1.0, 'scale3d':1.0} # m = 2, n = 2, aspect ration 1:1
# SF/SF, FF/FT, Aspect ratio 3:1
#eq_params = {'prandtl':1, 'rayleigh':660.518, 'scale1d':1.0/3.0, 'scale3d':1.0} # m = 2, n = 1, aspect ratio 3:1
#eq_params = {'prandtl':1, 'rayleigh':779.273, 'scale1d':1.0/3.0, 'scale3d':1.0} # m = 3, n = 1, aspect ratio 3:1
#eq_params = {'prandtl':1, 'rayleigh':1174.4, 'scale1d':1.0/3.0, 'scale3d':1.0} # m = 4, n = 1, aspect ratio 3:1
#eq_params = {'prandtl':1, 'rayleigh':1202.58, 'scale1d':1.0/3.0, 'scale3d':1.0} # m = 1, n = 1, aspect ratio 3:1
#eq_params = {'prandtl':1, 'rayleigh':10568.3, 'scale1d':1.0/3.0, 'scale3d':1.0} # m = 4, n = 2, aspect ratio 3:1
# SF/SF, FF/FT, Aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':133.62, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 1, n = 1, aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':293.563, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 1, n = 2, aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':779.273, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 1, n = 3, aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':1692.07, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 2, n = 1, aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':2087.81, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 1, n = 4, aspect ratio 1:3
#eq_params = {'prandtl':1, 'rayleigh':2137.92, 'scale1d':1.0, 'scale3d':1.0/3.0} # m = 2, n = 2, aspect ratio 1:3

# SF/NS, FF/FT
bc_vel = 1 
bc_temp = 1 
# SF/SF/NS, FF/FF/FT, Aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':1500.0, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)
eq_params = {'prandtl':1, 'rayleigh':2000.0, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp}

eigs = []

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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -2e0, np.inf)
    print("Found " + str(len(evp_lmb)) + " eigenvalues\n")

if show_solution:
    viz_mode = 1
    xscale = eq_params['scale1d']
    zscale = eq_params['scale3d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0]*res[2],mode]
        sol_w = evp_vec[res[0]*res[2]:2*res[0]*res[2],mode]

        # Extract continuity from velocity 
        sol_c = mod.c2d.d1(res[0], res[2], mod.no_bc(), xscale = xscale, sz = 0)*sol_u + mod.c2d.e1(res[0], res[2], mod.no_bc(), zscale = zscale, sx = 0)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0]*res[2],viz_mode]
    sol_w = evp_vec[res[0]*res[2]:2*res[0]*res[2],viz_mode]
    sol_t = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],viz_mode]
    sol_p = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],viz_mode]

    # Extract continuity from velocity 
    sol_c = mod.c2d.d1(res[0], res[2], mod.no_bc(), xscale = xscale, sz = 0)*sol_u + mod.c2d.e1(res[0], res[2], mod.no_bc(), zscale = zscale, sx = 0)*sol_w
    
    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_u))
    pl.title("u")
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
    mat_w = sol_w.reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_p = sol_p.reshape(res[0], res[2], order = 'F')
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')

    # Visualize spectrum matrix
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.colorbar()
    pl.title("u")
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
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])
    phys_u = transf.tophys2d(mat_u)
    phys_w = transf.tophys2d(mat_w)
    phys_t = transf.tophys2d(mat_t)
    phys_p = transf.tophys2d(mat_p)
    phys_c = transf.tophys2d(mat_c)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_z, grid_x, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,3)
    pl.contourf(grid_z, grid_x, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_z, grid_x, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_z, grid_x, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_z, grid_x, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.plot(grid_z, phys_u[res[0]//2,:], label = 'Z')
    pl.plot(grid_x, phys_u[:,res[1]//2], label = 'X')
    pl.legend(loc='upper left')
    pl.title("u")
    pl.subplot(2,3,3)
    pl.plot(grid_z, phys_w[res[0]//2,:], label = 'Z')
    pl.plot(grid_x, phys_w[:,res[1]//2], label = 'X')
    pl.legend(loc='upper left')
    pl.title("w")
    pl.subplot(2,3,4)
    pl.plot(grid_z, phys_t[res[0]//2,:], label = 'Z')
    pl.plot(grid_x, phys_t[:,res[1]//2], label = 'X')
    pl.legend(loc='upper left')
    pl.title("T")
    pl.subplot(2,3,6)
    pl.plot(grid_z, phys_p[res[0]//2,:], label = 'Z')
    pl.plot(grid_x, phys_p[:,res[1]//2], label = 'X')
    pl.legend(loc='upper left')
    pl.title("p")
    pl.show()
    pl.close("all")
