"""Script to run a marginal curve trace for the Boussinesq rotating convection in a disk model"""

import numpy as np

import quicc.model.boussinesq_rotconvdisk as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConvDisk()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [30, 0, 0]
eq_params = {'taylor':1e4, 'prandtl':1, 'rayleigh':5.5, 'ro':1, 'rratio':0.35}
eigs = [1]
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':0, 'velocityy':0, 'temperature':0}

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
if False:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)
