"""Script to run a marginal curve trace for the Boussinesq rotating convection in a box model (Velocity-Pressure)"""

import numpy as np

import quicc.model.boussinesq_rotconvboxvp as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConvBoxVP()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [10, 0, 10]
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.7135}
eigs = [1]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp, 'pressure':bc_vel}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
A = A.tolil()
A[4*res[0]*res[2],4*res[0]*res[2]] = 1
A[4*res[0]*res[2]+1,4*res[0]*res[2]+res[0]-1] = 1
A[4*res[0]*res[2]+res[0],-res[0]] = 1
A[4*res[0]*res[2]+res[0]+1,-1] = 1

#A[4*res[0]*res[2],4*res[0]*res[2]] = 1
#A[4*res[0]*res[2]+1,4*res[0]*res[2]+1] = 1
#A[4*res[0]*res[2]+res[0],4*res[0]*res[2]+res[0]] = 1
#A[4*res[0]*res[2]+res[0]+1,4*res[0]*res[2]+res[0]+1] = 1
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
    sol_v = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_t = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_p = evp_vec[4*res[0]*res[2]:,-1].reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(1,5,1)
    pl.imshow(np.log10(np.abs(sol_u)))
    pl.colorbar()
    pl.subplot(1,5,2)
    pl.imshow(np.log10(np.abs(sol_v)))
    pl.colorbar()
    pl.subplot(1,5,3)
    pl.imshow(np.log10(np.abs(sol_w)))
    pl.colorbar()
    pl.subplot(1,5,4)
    pl.imshow(np.log10(np.abs(sol_t)))
    pl.colorbar()
    pl.subplot(1,5,5)
    pl.imshow(np.log10(np.abs(sol_p)))
    pl.colorbar()
    pl.show()
    pl.close("all")

    import quicc.transform.cartesian as transf
    phys_u = transf.tophys2d(sol_u)
    phys_v = transf.tophys2d(sol_v)
    phys_w = transf.tophys2d(sol_w)
    phys_t = transf.tophys2d(sol_t)
    phys_p = transf.tophys2d(sol_p)
    pl.subplot(1,5,1)
    pl.imshow(sol_u.real)
    pl.colorbar()
    pl.subplot(1,5,2)
    pl.imshow(sol_v.real)
    pl.colorbar()
    pl.subplot(1,5,3)
    pl.imshow(sol_w.real)
    pl.colorbar()
    pl.subplot(1,5,4)
    pl.imshow(sol_t.real)
    pl.colorbar()
    pl.subplot(1,5,5)
    pl.imshow(sol_p.real)
    pl.colorbar()
    pl.show()
    pl.close("all")

    pl.subplot(1,5,1)
    pl.imshow(sol_u.imag)
    pl.colorbar()
    pl.subplot(1,5,2)
    pl.imshow(sol_v.imag)
    pl.colorbar()
    pl.subplot(1,5,3)
    pl.imshow(sol_w.imag)
    pl.colorbar()
    pl.subplot(1,5,4)
    pl.imshow(sol_t.imag)
    pl.colorbar()
    pl.subplot(1,5,5)
    pl.imshow(sol_p.imag)
    pl.colorbar()
    pl.show()
    pl.close("all")
