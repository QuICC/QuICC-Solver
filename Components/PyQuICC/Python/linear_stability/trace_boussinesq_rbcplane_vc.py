"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a plane layer (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np
import functools

import quicc.model.boussinesq_rbcplane_vc as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCPlaneVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [7, 0, 0]

# SF, FT,
bc_vel = 0
bc_temp = 0
heating = 0
#phi = 0
#kp = 3.
#eq_params = {'prandtl':1, 'rayleigh':667.0098243, 'heating':0, 'scale1d':2.0}
#eq_params = {'prandtl':1, 'rayleigh':1285, 'heating':0, 'scale1d':2.0}
#eq_params = {'prandtl':1, 'rayleigh':675, 'heating':0, 'scale1d':2.0}
## kx = 0, ky = 2
#kx = 0
#ky = 2
#eq_params = {'prandtl':1, 'rayleigh':667.0098243, 'scale1d':2.0}
# kx = 0, ky = 1
#kx = 0
#ky = 1
#eq_params = {'prandtl':1, 'rayleigh':1284.225280, 'scale1d':2.0}
## kx = 0, ky = 5.35
#kx = 0
#ky = 5.35
#eq_params = {'prandtl':1, 'rayleigh':1992.541617, 'scale1d':2.0}
## kx = 2, ky = 0
#kx = 2
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':667.0098243, 'scale1d':2.0}
## kx = 1, ky = 0
#kx = 1
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':1284.225280, 'scale1d':2.0}
## kx = 5.35, ky = 0
#kx = 5.35
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':1992.541617, 'scale1d':2.0}
## Minimum n = 1: k_ = 2.221441469
#phi = 35
#kp = 2.221441469
#eq_params = {'prandtl':1, 'rayleigh':657.5113645, 'scale1d':2.0}
## Minimum n = 2: k_ = 4.442882938
#phi = 15
#kp = 4.442882938
#eq_params = {'prandtl':1, 'rayleigh':10520.18183, 'scale1d':2.0}
#kx = 1.31
#ky = 1.52
#eq_params = {'prandtl':1, 'rayleigh':833.12, 'scale1d':2.0}
# Minimum n = 2: k_ = 4.442882938
#phi = 35
#kp = np.sqrt(2)*np.pi
#eq_params = {'prandtl':1, 'rayleigh':108*np.pi**4, 'heating':0, 'scale1d':2.0}

#kp = 3
#eq_params = {'prandtl':1, 'rayleigh':2.19e4, 'heating':1, 'scale1d':2.0}
#kp = 4
#eq_params = {'prandtl':1, 'rayleigh':1.9137e4, 'heating':1, 'scale1d':2.0}
#kp = 4.05
#eq_params = {'prandtl':1, 'rayleigh':1.9119e4, 'heating':1, 'scale1d':2.0}
#kp = 4.09
#eq_params = {'prandtl':1, 'rayleigh':1.9111e4, 'heating':1, 'scale1d':2.0}
#kp = 4.1
#eq_params = {'prandtl':1, 'rayleigh':1.91099e4, 'heating':1, 'scale1d':2.0}
#kp = 4.11
#eq_params = {'prandtl':1, 'rayleigh':1.91091e4, 'heating':1, 'scale1d':2.0}
#kp = 4.12
#eq_params = {'prandtl':1, 'rayleigh':1.91086e4, 'heating':1, 'scale1d':2.0}
#kp = 4.1296
#eq_params = {'prandtl':1, 'rayleigh':1.910844559e4, 'heating':1, 'scale1d':2.0}
#kp = 4.13
#eq_params = {'prandtl':1, 'rayleigh':1.910844591e4, 'heating':1, 'scale1d':2.0}
#kp = 4.14
#eq_params = {'prandtl':1, 'rayleigh':1.910863e4, 'heating':1, 'scale1d':2.0}
#kp = 4.15
#eq_params = {'prandtl':1, 'rayleigh':1.91092e4, 'heating':1, 'scale1d':2.0}
#kp = 4.25
#eq_params = {'prandtl':1, 'rayleigh':1.9132e4, 'heating':1, 'scale1d':2.0}
#kp = 4.5
#eq_params = {'prandtl':1, 'rayleigh':1.932e4, 'heating':1, 'scale1d':2.0}
#kp = 5
#eq_params = {'prandtl':1, 'rayleigh':2.018e4, 'heating':1, 'scale1d':2.0}
#phi = 0

#kp = 4.1296
#phi = 0
#eq_params = {'prandtl':1, 'rayleigh':5e6, 'heating':1, 'scale1d':2.0}

kp = np.sqrt(2)*np.pi
Pr = 1; Ra = 108*np.pi**4; phi = 35

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'heating':heating, 'scale1d':2.0}

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    if phi == 90:
        kx = 0
        ky = kp
    elif phi == 0:
        kx = kp
        ky = 0
    else:
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
marginal_options['ellipse_radius'] = 1e5
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = kp
marginal_options['plot_point'] = True
marginal_options['plot_curve'] = True
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['write_mtx'] = True
marginal_options['viz_mode'] = 2
marginal_options['curve_points'] = np.arange(max(1, kp-2), kp+3, 0.1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
