"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a 3D box (streamfunction formulation)"""

import numpy as np

import quicc.model.boussinesq_rbcbox_st as mod

# Create the model and activate linearization
model = mod.BoussinesqRBCBoxST()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [16, 16, 16]

# SF/SF/SF, FF/FF/FT
#bc_str = 6 
#bc_temp = 4 
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':779.2727283, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1|0, m = 0|1, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':1315.022729, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1, m = 1, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':3044.034095, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 2|0, m = 0|2, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':4208.072733, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 2|1, m = 1|2, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':8876.403420, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 2, m = 2, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':10520.18183, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1, m = 1, n = 2, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':10823.23234, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 3|0, m = 0|3, n = 1, aspect ration 1:1:1
#eq_params = {'prandtl':1, 'rayleigh':12965.15002, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # l = 3|1, m = 1|3, n = 1, aspect ration 1:1:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':660.5178179, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 2, m = 0, n = 1, aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':779.2727283, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 0|3, m = 1|0, n = 1, aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':824.8505622, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1, m = 0, n = 1, aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':985.0066489, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 2, m = 1, n = 1, aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':1174.395870, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 4, m = 0, n = 1, aspect ratio 3:1:1
#eq_params = {'prandtl':1, 'rayleigh':1202.581371, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # l = 1, m = 0, n = 1, aspect ratio 3:1:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':660.5178179, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 0, m = 2, n = 1, aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':779.2727283, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 1|0, m = 0|3, n = 1, aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':824.8505622, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 0, m = 1, n = 1, aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':985.0066489, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 1, m = 2, n = 1, aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':1174.395870, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 0, m = 4, n = 1, aspect ratio 1:3:1
#eq_params = {'prandtl':1, 'rayleigh':1202.581371, 'scale1d':1.0, 'scale2d':1.0/3.0, 'scale3d':1.0} # l = 0, m = 1, n = 1, aspect ratio 1:3:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':133.6201523, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 0|1, m = 1|0, n = 1, aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':293.5634746, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 0|1, m = 1|0, n = 2, aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':458.2503123, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 1, m = 1, n = 1, aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':711.3936909, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 1, m = 1, n = 2, aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':779.2727283, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 0|1, m = 1|0, n = 3, aspect ratio 1:1:3
#eq_params = {'prandtl':1, 'rayleigh':1315.022729, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0/3.0} # l = 1, m = 1, n = 3, aspect ratio 1:1:3

# SF/SF/NS, FF/FF/FT
bc_str = 4 
bc_temp = 4 
# SF/SF/NS, FF/FF/FT, Aspect ratio 3:1:1
eq_params = {'prandtl':1, 'rayleigh':1500.0, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)
#eq_params = {'prandtl':1, 'rayleigh':2000.0, 'scale1d':1.0/3.0, 'scale2d':1.0, 'scale3d':1.0} # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':bc_str, 'temperature':bc_temp}

eigs = []

# Generate the operator A for the generalized EVP Ax = sigm B x
print("Constructing matrix A")
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
print("Constructing matrix B")
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl
    #import mayavi.mbab as ml

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
    print("Writing matrices to file")
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    print("Solve EVP")
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, np.inf)
    print(evp_lmb)

if show_solution:
    viz_mode = -1
    xscale = eq_params['scale1d']
    yscale = eq_params['scale2d']
    zscale = eq_params['scale3d']

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_s = evp_vec[0:res[0]*res[1]*res[2],viz_mode]
    sol_t = evp_vec[res[0]*res[1]*res[2]:2*res[0]*res[1]*res[2],viz_mode]
    
    # Create spectrum plots
    pl.subplot(1,2,1)
    pl.semilogy(np.abs(sol_s))
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.show()
    pl.close("all")
    
    # Create solution matrices
    mat_s = sol_s.reshape(res[0], res[2], res[1], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], res[1], order = 'F')

    # Compute physical space values
    grid_x = transf.grid(res[0])
    grid_y = transf.grid(res[1])
    grid_z = transf.grid(res[2])
    phys_s = transf.tophys3d(mat_s)
    phys_t = transf.tophys3d(mat_t)

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_y, grid_z, phys_s[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_y, grid_z, phys_t[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_y, grid_x, phys_s[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_y, grid_x, phys_t[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_z, grid_x, phys_s[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_z, grid_x, phys_t[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")
