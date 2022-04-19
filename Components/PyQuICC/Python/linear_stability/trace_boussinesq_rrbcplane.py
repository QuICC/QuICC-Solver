"""Script to run a marginal curve trace for the rescaled Boussinesq rotating Rayleigh-Benard Boussinesq model (Toroidal/Poloidal formulation)"""

import numpy as np
import functools

import quicc.model.boussinesq_rrbcplane as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRRBCPlane()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [512, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.69672, 'ekman':1e-7, 'scale1d':2.0, 'fast_mean':0, 'rescaled':1}
eq_params = {'prandtl':1, 'rayleigh':20, 'ekman':1e-8, 'scale1d':2.0, 'fast_mean':0, 'rescaled':1}
auto_params = model.automatic_parameters(eq_params)
for k,v in auto_params.items():
    eq_params[k] = v
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':2, 'temperature':0}
phi = 0
kp = 1.30476
kpm = 0.02
kpp = 0.03
dkp = 0.0035
if eq_params['rescaled'] == 0:
    kp *= eq_params['ekman']**(-1./3.)
    kpm *= eq_params['ekman']**(-1./3.)
    kpp *= eq_params['ekman']**(-1./3.)
    dkp *= eq_params['ekman']**(-1./3.)
    eq_params['rayleigh'] *= eq_params['ekman']**(-4./3.)

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
marginal_options['evp_tol'] = 1e-12
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['solve_nev'] = 3
marginal_options['point_k'] = kp
marginal_options['plot_curve'] = False
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['save_spectra'] = False
marginal_options['show_physical'] = True
marginal_options['save_physical'] = False
marginal_options['save_pdf'] = False
marginal_options['write_mtx'] = True
marginal_options['viz_mode'] = 0
marginal_options['curve_points'] = np.arange(max(0, kp-kpm), kp+kpp, dkp)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
