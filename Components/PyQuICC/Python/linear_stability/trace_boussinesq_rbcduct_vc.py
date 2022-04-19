"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a infinite duct(1 periodic direction) (velocity-continuity formulation)"""

import numpy as np

#import quicc.model.boussinesq_rbcduct_vc_diff as mod
import quicc.model.boussinesq_rbcduct_vc2 as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCDuctVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#res = [8, 0, 8]
#res = [6, 0, 6]
#res = [12, 0, 12]
#res = [16, 0, 16]
#res = [24, 0, 24]
#res = [32, 0, 32]
#res = [36, 0, 36]
#res = [48, 0, 48]
#res = [64, 0, 64]
res = [128, 0, 128]
#res = [256, 0, 256]

# SF/SF, FF/FT, k = 0
#bc_vel = 2
#bc_temp = 1
#heating = 0
#k = 0
# SF/SF, FF/FT, Aspect ratio 1:1
#Pr = 1; Ra = 779.273; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
#Pr = 1; Ra = 3044.03; A1d = 1.0; A3d = 1.0 # m = 2, n = 1, aspect ration 1:1
#Pr = 1; Ra = 10823.2; A1d = 1.0; A3d = 1.0 # m = 3, n = 1, aspect ration 1:1
#Pr = 1; Ra = 12176.1; A1d = 1.0; A3d = 1.0 # m = 1, n = 2, aspect ration 1:1
#Pr = 1; Ra = 12468.4; A1d = 1.0; A3d = 1.0 # m = 2, n = 2, aspect ration 1:1
# SF/SF, FF/FT, Aspect ratio 3:1
#Pr = 1; Ra = 660.518; A1d = 1.0/3.0; A3d = 1.0 # m = 2, n = 1, aspect ratio 3:1
#Pr = 1; Ra = 779.273; A1d = 1.0/3.0; A3d = 1.0 # m = 3, n = 1, aspect ratio 3:1
#Pr = 1; Ra = 1174.40; A1d = 1.0/3.0; A3d = 1.0 # m = 4, n = 1, aspect ratio 3:1
#Pr = 1; Ra = 1202.58; A1d = 1.0/3.0; A3d = 1.0 # m = 1, n = 1, aspect ratio 3:1
#Pr = 1; Ra = 10568.3; A1d = 1.0/3.0; A3d = 1.0 # m = 4, n = 2, aspect ratio 3:1
# SF/SF, FF/FT, Aspect ratio 1:3
#Pr = 1; Ra = 133.620; A1d = 1.0; A3d = 1.0/3.0 # m = 1, n = 1, aspect ratio 1:3
#Pr = 1; Ra = 293.563; A1d = 1.0; A3d = 1.0/3.0 # m = 1, n = 2, aspect ratio 1:3
#Pr = 1; Ra = 779.273; A1d = 1.0; A3d = 1.0/3.0 # m = 1, n = 3, aspect ratio 1:3
#Pr = 1; Ra = 1692.07; A1d = 1.0; A3d = 1.0/3.0 # m = 2, n = 1, aspect ratio 1:3
#Pr = 1; Ra = 2087.81; A1d = 1.0; A3d = 1.0/3.0 # m = 1, n = 4, aspect ratio 1:3
#Pr = 1; Ra = 2137.92; A1d = 1.0; A3d = 1.0/3.0 # m = 2, n = 2, aspect ratio 1:3

# PAPER
bc_vel = 2
bc_temp = 1
heating = 0
#k = np.pi
#Pr = 1; Ra = 108*np.pi**4; A1d = 1.0; A3d = 1.0/2.0
#Pr = 7; Ra = 116.83858990055; A3d = 1.0/4.0
#Pr = 7; Ra = 800; A3d = 1.0
#Pr = 7; Ra = 116.83858990055; A3d = 4.0
k = (np.sqrt(7.0)/2.0)*np.pi
#k = 0
Pr = 7; Ra = 36*np.pi**4; A1d = 1.0/2.0; A3d = 1.0
#Pr = 1; Ra = 500*np.pi**4; A1d = 2.0; A3d = 1.0
#k = 0
#Pr = 7; Ra = 8*np.pi**4; A1d = 1.0/8.0; A3d = 1.0
#k = 0
#Pr = 1; Ra = 8000*np.pi**4; A1d = 8.0; A3d = 1.0

# SF/SF, FF/FT, k = 1.0
#bc_vel = 2
#bc_temp = 1
#heating = 0
#k = 1.0
# SF/SF, FF/FT, Aspect ratio 1:1
#Pr = 1; Ra = 820.6591462; A1d = 1.0; A3d = 1.0 # l = 1, n = 1, aspect ration 1:1
#Pr = 1; Ra = 1284.225280; A1d = 1.0; A3d = 1.0 # l = 0, n = 1, aspect ration 1:1
#Pr = 1; Ra = 3152.998132; A1d = 1.0; A3d = 1.0 # l = 2, n = 1, aspect ration 1:1
#Pr = 1; Ra = 11031.37354; A1d = 1.0; A3d = 1.0 # l = 3, n = 1, aspect ration 1:1
#Pr = 1; Ra = 11741.76818; A1d = 1.0; A3d = 1.0 # l = 1, n = 2, aspect ration 1:1
#Pr = 1; Ra = 12628.25262; A1d = 1.0; A3d = 1.0 # l = 2, n = 2, aspect ration 1:1

# SF/SF, FF/FT, k = 5.35
#bc_vel = 2
#bc_temp = 1
#heating = 0
#k = 5.35
#k = 3.0
# SF/SF, FF/FT, Aspect ratio 1:1
#Pr = 1; Ra = 1992.541617; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
#Pr = 1; Ra = 2938.551173; A1d = 1.0; A3d = 1.0 # m = 1, n = 1, aspect ration 1:1
#Pr = 1; Ra = 6960.466725; A1d = 1.0; A3d = 1.0 # m = 2, n = 1, aspect ration 1:1
#Pr = 1; Ra = 17572.19000; A1d = 1.0; A3d = 1.0 # m = 3, n = 1, aspect ration 1:1
#Pr = 1; Ra = 12314.58187; A1d = 1.0; A3d = 1.0 # m = 1, n = 2, aspect ration 1:1
#Pr = 1; Ra = 18282.41677; A1d = 1.0; A3d = 1.0 # m = 2, n = 2, aspect ration 1:1

# SF/NS, FF/FT
#bc_vel = 1 
#bc_temp = 1 
#heating = 0
# SF/SF/NS, FF/FF/FT, Aspect ratio 3:1:1
#Pr = 1; Ra = 1500.0; A1d = 1.0/3.0; A3d = 1.0 # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)
#Pr = 1; Ra = 2000.0; A1d = 1.0/3.0; A3d = 1.0 # Burroughs, Romero, Lehoucq, Salinger, 2001 (WARNING different scaling!)

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'heating':heating, 'scale1d':2.0, 'scale3d':2.0*A3d} #

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
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['solve_nev'] = 5
marginal_options['point_k'] = k
marginal_options['plot_point'] = True
marginal_options['plot_curve'] = False
marginal_options['plot_spy'] = True
marginal_options['write_mtx'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['viz_mode'] = 2
marginal_options['curve_points'] = np.arange(max(1, k-2), k+2, 0.2)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)

