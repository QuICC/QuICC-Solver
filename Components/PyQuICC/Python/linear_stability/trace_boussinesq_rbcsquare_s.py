"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a square (2D) (vorticity-streamfunction formulation)"""

import numpy as np

import quicc.model.boussinesq_rbcsquare_s as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCSquareS()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#res = [8, 8]
#res = [10, 10]
#res = [24, 24]
#res = [32, 32]
#res = [64, 64]
res = [128, 128]

# NS/NS, FF/FT
bc_stream = 0
bc_temp = 1
heating = 0
k = 0
Pr = 7; Ra = 552.438061798; A3d = 0.25
#Pr = 7; Ra = 2585.01869147; A3d = 1.0
#Pr = 7; Ra = 463428.864042; A3d = 4.0

# SF/SF, FF/FT
#bc_stream = 2
#bc_temp = 1
#heating = 0
#k = 0
#Pr = 7; Ra = 116.838589902; A3d = 0.25
#Pr = 7; Ra = 779.272728275; A3d = 1.0
#Pr = 7; Ra = 169113.005268; A3d = 4.0

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'heating':heating, 'scale1d':2.0, 'scale2d':2.0*A3d} #

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':bc_stream, 'temperature':bc_temp}

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
#marginal_options['ellipse_radius'] = 1e3
marginal_options['geometry'] = 'c2d'
marginal_options['point'] = False
marginal_options['solve'] = True
marginal_options['solve_nev'] = 5
marginal_options['point_k'] = k
marginal_options['plot_point'] = True
marginal_options['plot_curve'] = False
marginal_options['plot_spy'] = True
marginal_options['write_mtx'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['viz_mode'] = 4
marginal_options['curve_points'] = np.arange(max(1, k-2), k+3, 0.1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)

