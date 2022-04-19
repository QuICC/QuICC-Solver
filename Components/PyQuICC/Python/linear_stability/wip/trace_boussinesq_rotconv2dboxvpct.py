"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (velocity-pressure-continuity-temperature)"""

import numpy as np

import quicc.model.boussinesq_rotconv2dboxvpct as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConv2DBoxVPCT()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':2340.687}
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.73}
eigs = [1]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bc_cont = 0

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityz':bc_vel, 'pressure':bc_vel, 'temperature':bc_temp, 'continuity':bc_cont}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
A = A.tolil()
A[3*res[0]*res[2],:] = 0
A[3*res[0]*res[2],3*res[0]*res[2]] = 1
A[3*res[0]*res[2]+1,:] = 0
A[3*res[0]*res[2]+1,3*res[0]*res[2]+res[0]-1] = 1
A[3*res[0]*res[2]+res[0],:] = 0
A[3*res[0]*res[2]+res[0],-res[0]] = 1
A[3*res[0]*res[2]+res[0]+1,:] = 0
A[3*res[0]*res[2]+res[0]+1,-1] = 1
#A[4*res[0]*res[2]-2*res[0],:] = 0
#A[4*res[0]*res[2]-2*res[0],3*res[0]*res[2]] = 1
#A[4*res[0]*res[2]-2*res[0] + 1,:] = 0
#A[4*res[0]*res[2]-2*res[0] + 1,3*res[0]*res[2]+res[0]-1] = 1
#A[4*res[0]*res[2]-res[0],:] = 0
#A[4*res[0]*res[2]-res[0],-res[0]] = 1
#A[4*res[0]*res[2]-res[0] + 1,:] = 0
#A[4*res[0]*res[2]-res[0] + 1,-1] = 1
A = A.tocsr()

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

    sol_u = evp_vec[0:res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_w = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_t = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_p = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_c = evp_vec[4*res[0]*res[2]:5*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_cc = (mod.c2d.d1d0(res[0], res[2], 0, mod.no_bc())*sol_u.reshape(res[0]*res[2], order ='F') + mod.c2d.d0d1(res[0], res[2], 0, mod.no_bc())*sol_w.reshape(res[0]*res[2], order ='F')).reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(sol_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(sol_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(sol_t)))
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,4)
    pl.imshow(np.log10(np.abs(sol_c)))
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(sol_cc)))
    pl.colorbar()
    pl.title("Continuity (calc)")
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(sol_p)))
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    import quicc.transform.cartesian as transf
    phys_u = transf.tophys2d(sol_u)
    phys_w = transf.tophys2d(sol_w)
    phys_t = transf.tophys2d(sol_t)
    phys_p = transf.tophys2d(sol_p)
    phys_c = transf.tophys2d(sol_c)
    phys_cc = transf.tophys2d(sol_cc)
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,4)
    pl.contourf(grid_x, grid_z, phys_c, 50)
    #pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, phys_cc, 50)
    #pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_cc)), 50)
    pl.colorbar()
    pl.title("Continuity(calc)")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
