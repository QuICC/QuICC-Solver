"""Script to run a marginal curve trace for the Boussinesq rotating convection in a cylindrical annulus model"""

import numpy as np

import quicc.model.boussinesq_rotconvannulus as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConvAnnulus()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [10, 0, 10]
eq_params = {'taylor':1e4, 'prandtl':1, 'rayleigh':5.5, 'ro':1, 'r_ratio':0.35}
eigs = [1]
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity_x':1, 'velocity_y':1, 'velocity_z':1, 'temperature':0, 'pressure':1}

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
