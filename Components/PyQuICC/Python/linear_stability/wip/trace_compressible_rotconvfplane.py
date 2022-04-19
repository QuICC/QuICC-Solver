"""Script to run a marginal curve trace for the compressible convection in a rotating F-plane model"""

import quicc.model.compressible_rotconvfplane as mod

# Create the model and activate linearization
model = mod.CompressibleRotConvFPlane()
model.linearize = True
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [50, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':4.761e6, 'taylor':1e4, 'density_scales':1, 'polytropic_index':1.49, 'theta':1, 'gamma':5.0/3.0}
eigs = [1, 1]
bcs = {'bcType':0, 'velocityx':0, 'velocityy':0, 'velocityz':0, 'pressure':0, 'density':0, 'temperature':0, 'entropy':0}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = 2
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
#import matplotlib.pylab as pl
#pl.spy(A, markersize=0.2)
#pl.show()
#pl.spy(B, markersize=0.2)
#pl.show()

# Export the two matrices to matrix market format
import scipy.io as io
io.mmwrite("matrix_A.mtx", A)
io.mmwrite("matrix_B.mtx", B)
