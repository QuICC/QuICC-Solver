"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a cylindrical annulus model (velocity-continuity formulation)"""

import numpy as np

import quicc.model.boussinesq_rbcannulus_vc as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCAnnulusVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

# Create parameters
m = 7 
res = [64, 0, 64]
eq_params = {'prandtl':1, 'rayleigh':5e3, 'ro':1.0, 'r_ratio':0.35, 'heating':0, 'scale3d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-16
marginal_options['geometry'] = 'annulus'
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['minimum_int'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
