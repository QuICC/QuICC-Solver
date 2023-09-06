"""Script to run a marginal curve trace for the Boussinesq tilted F-plane 3DQG model"""

import numpy as np
import functools

import quicc.model.boussinesq_tilted_fplane3dqg as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTiltedFPlane3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]
kp = 1.3048074030314.; Ra = 8.6956307154148
eq_params = {'prandtl':1, 'rayleigh':100, 'theta':0.0, 'scale1d':2.0}
#kp = 1.2898160619671; Ra = 8.3028296084564
#eq_params = {'prandtl':1, 'rayleigh':Ra, 'theta':0.0, 'scale1d':2.0}
#kp = 1.2437150207239; Ra = 7.178085000681
#eq_params = {'prandtl':1, 'rayleigh':Ra, 'theta':30.0, 'scale1d':2.0}
#kp = 1.1624855229293; Ra = 5.4779041217364
#eq_params = {'prandtl':1, 'rayleigh':Ra, 'theta':45.0, 'scale1d':2.0}

# Set wave number
phi = 90

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocity_z':0, 'temperature':0}

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
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['solve'] = False
marginal_options['point_k'] = kp
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = False
marginal_options['show_physical'] = False
marginal_options['write_mtx'] = True
marginal_options['curve_points'] = np.arange(kp-0.05, kp+0.06, 0.01)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
