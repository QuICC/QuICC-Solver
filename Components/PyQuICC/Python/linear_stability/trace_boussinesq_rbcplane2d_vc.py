"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a plane layer (2D) (Velocity-continuity formulation)"""

import numpy as np
import functools

import quicc.model.boussinesq_rbcplane2d_vc as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCPlane2DVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [32, 0]

# SF, FT,
bc_vel = 0
bc_temp = 0
heating = 0

k = np.sqrt(2)*np.pi
Pr = 1; Ra = 108*np.pi**4

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'heating':heating, 'scale1d':2.0}

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def wave(k):
    return [k]

# Wave number function from single "index"
eigs = wave(k)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['geometry'] = 'c1d'
marginal_options['ellipse_radius'] = 1e5
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['solve_nev'] = 10
marginal_options['point_k'] = k
marginal_options['plot_point'] = True
marginal_options['plot_curve'] = True
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['write_mtx'] = True
marginal_options['viz_mode'] = 1
marginal_options['curve_points'] = np.arange(max(1, k-4), k+5, 0.3)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
