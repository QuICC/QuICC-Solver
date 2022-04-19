"""Script to run a marginal curve trace for the tilted Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np
import functools

import quicc.model.boussinesq_tilted_rrbcplane_vc as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTiltedRRBCPlaneVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
Rac = None
kc = None

# SF, FT,
#bc_vel = 1
#heating = 0
#bc_temp = 0
#phi = 35

#phi = 45
#kp = 3.710
#kp = 129
#kp = 60

# NS, FT,
bc_vel = 0
heating = 0
bc_temp = 0
phi = 0
theta = 15

res = [512, 0, 0]
Ta = 1e8;

# Create parameters (rescaling to proper nondimensionalisation)
if kc is None:
    kp = 1.3*Ta**(1./6.) # Asymptotic prediction for minimum
else:
    kp = kc
if Rac is None:
    Ra = 9*Ta**(2./3.) # Asymptotic prediction for critical Rayleigh number
else:
    Ra = Rac

eq_params = {'taylor':Ta, 'prandtl':1, 'rayleigh':Ra, 'theta':theta, 'heating':heating, 'scale1d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(kp)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-12
#marginal_options['euler'] = (int(2e0), 1e-2)
#marginal_options['conv_idx'] = int(0.1*res[0])
#marginal_options['spectrum_conv'] = 1e4
#marginal_options['ellipse_radius'] = 1e5
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['solve'] = True
#marginal_options['solve_nev'] = 5
marginal_options['root_tol'] = 1e-12
marginal_options['evp_tol'] = 1e-12
marginal_options['point_k'] = kp
marginal_options['plot_point'] = True
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['save_hdf5'] = False
marginal_options['save_spectra'] = False
marginal_options['save_physical'] = False
marginal_options['curve_points'] = np.arange(max(0, kp-10), kp+11, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
