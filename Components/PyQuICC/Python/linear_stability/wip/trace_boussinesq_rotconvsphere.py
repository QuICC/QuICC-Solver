"""Script to run a marginal curve trace for the Boussinesq convection in a rotating sphere"""

import quicc.model.boussinesq_rotconvsphere as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConvSphere()
model.linearize = True
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
l = 0
m = 6
res = [30,m+20, 0]
eigs = [l, m]
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':4.33385e4}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':0, 'temperature':0}

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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, 1)
    print(evp_lmb)
