"""Script to run a marginal curve trace for the rescaled Boussinesq rotating Rayleigh-Benard Boussinesq model"""

import numpy as np
import functools

import quicc.model.boussinesq_rescaledrrbcplane as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRescaledRRBCPlane()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [1024, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'ekman':1e-11, 'scale1d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':1, 'velocityz':0, 'temperature':0, 'pressure':1}
phi = 0
kp = 1.3048

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
#marginal_options['curve'] = False
marginal_options['ellipse_radius'] = 1e2
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['solve_nev'] = 5
marginal_options['point_k'] = kp
marginal_options['plot_curve'] = False
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = False
marginal_options['save_spectra'] = True
marginal_options['show_physical'] = False
marginal_options['save_physical'] = True
marginal_options['save_pdf'] = True
marginal_options['viz_mode'] = -1
marginal_options['curve_points'] = np.arange(max(0, kp-0.02), kp+0.021, 0.01)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
