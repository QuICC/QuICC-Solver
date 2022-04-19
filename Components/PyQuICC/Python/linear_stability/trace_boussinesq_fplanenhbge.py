"""Script to run a marginal curve trace for the Boussinesq F-plane NHBGE model"""

import numpy as np
import functools

import quicc.model.boussinesq_fplanenhbge as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqFPlaneNHBGE()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [128, 0, 0]
#eq_params = {'eady':-1, 'prandtl':1.0, 'rayleigh':8.69, 'ekman':1e-9, 'mean_flow':0, 'scale1d':2.0} # No mean flow, unstable
eq_params = {'eady':1, 'prandtl':1.0, 'rayleigh':13.53, 'ekman':1e-9, 'gamma':1, 'scale1d':2.0} # Mean flow, stable
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}
phi = 0
kp = 1.304807403014

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(1.0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = kp
#marginal_options['target'] = 0.1
marginal_options['plot_point'] = True
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['write_mtx'] = True
marginal_options['curve_points'] = np.arange(kp-0.05, kp+0.05, 0.01)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
