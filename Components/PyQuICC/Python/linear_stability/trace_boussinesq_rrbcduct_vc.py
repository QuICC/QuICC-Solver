"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a infinite duct (1 periodic direction) (velocity-continuity formulation)"""

import numpy as np

import quicc.model.boussinesq_rrbcduct_vc2 as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRRBCDuctVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#res = [8, 0, 8]
#res = [16, 0, 16]
res = [32, 0, 32]
#res = [48, 0, 48]
#res = [64, 0, 64]

# SF/SF, FF/FT
bc_vel = 2
bc_temp = 1
heating = 0
k = 5.7300054934036
# SF/SF, FF/FT, Aspect ratio 1:1
Pr = 1; Ra = 275718.48221715; Ta = 1e8; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
Pr = 1; Ra = 2757.48221715; Ta = 1e8; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
#k = 10
#Pr = 1; Ra = 1e6; Ta = 1e8; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1

# NS/NS, FT/FT, k = 0
bc_vel = 2
bc_temp = 1
heating = 0
k = 5.7281101194
Pr = 1; Ra = 1152782.00348; Ta = 1e8; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
Pr = 1; Ra = 11.00348; Ta = 1e8; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'taylor':Ta, 'heating':heating, 'scale1d':2.0*A1d, 'scale3d':2.0*A3d} # m = 1, n = 1, aspect ration 1:1

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def wave(k):
    return [k]

# Wave number function from single "index" (k perpendicular)
eigs = wave(k)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['mode'] = 0
marginal_options['ellipse_radius'] = 1e5
marginal_options['geometry'] = 'c2d'
marginal_options['point'] = False
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = k
marginal_options['plot_point'] = True
marginal_options['plot_curve'] = True
marginal_options['plot_spy'] = True
marginal_options['write_mtx'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['viz_mode'] = 0
marginal_options['curve_points'] = np.arange(max(1, k-5), k+5, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)

