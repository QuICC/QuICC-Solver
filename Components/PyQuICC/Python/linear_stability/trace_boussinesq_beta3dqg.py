"""Script to run a marginal curve trace for the Boussinesq Beta 3DQG model"""

import numpy as np

import quicc.model.boussinesq_beta3dqg as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqBeta3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [16, 0, 16]
chi = 1
eq_params = {'prandtl':1, 'rayleigh':700.0, 'gamma':1, 'chi':chi}
eigs = [2.221403788]

# No-slip/No-slip
#bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'temperature':0, 'velocityz':1} 
# Stress-free/Stress-free
#bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':1, 'temperature':0, 'velocityz':1}
# Stress-free/No-slip (simple Beta)
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':2, 'temperature':0, 'velocityz':2, 'vorticityz':2}

# Wave number function from single "index"
def wave(k):
    return [k]

eigs = wave(1)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-16
marginal_options['geometry'] = '2d'
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['minimum_int'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, k-2), k+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
