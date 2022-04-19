"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)"""

import numpy as np

import quicc.model.boussinesq_rbcbox_vc as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCBoxVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#res = [4, 4, 4]
res = [12, 12, 12]
#res = [14, 14, 14]
#res = [16, 16, 16]
#res = [18, 18, 18]
#res = [20, 20, 20]
#res = [20, 20, 18]
#res = [18, 18, 18]
#res = [32, 32, 32]

# SF/SF/SF, FF/FF/FT
#bc_vel = 6 
#bc_temp = 4 
#heating = 0
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:1:1
#Pr = 1; Ra = 779.2727283; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 1|0, m = 0|1, n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 1315.022729; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 1  , m = 1  , n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 3044.034095; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 2|0, m = 0|2, n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 4208.072733; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 2|1, m = 1|2, n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 8876.403420; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 2  , m = 2  , n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 10520.18183; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 1  , m = 1  , n = 2, aspect ration 1:1:1
#Pr = 1; Ra = 10823.23234; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 3|0, m = 0|3, n = 1, aspect ration 1:1:1
#Pr = 1; Ra = 12965.15002; A1d = 1.0; A2d = 1.0; A3d = 1.0 # l = 3|1, m = 1|3, n = 1, aspect ration 1:1:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 3:1:1
#Pr = 1; Ra = 660.5178179; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 2  , m = 0  , n = 1, aspect ratio 3:1:1
#Pr = 1; Ra = 779.2727283; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 0|3, m = 1|0, n = 1, aspect ratio 3:1:1
#Pr = 1; Ra = 824.8505622; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 1  , m = 0  , n = 1, aspect ratio 3:1:1
#Pr = 1; Ra = 985.0066489; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 2  , m = 1  , n = 1, aspect ratio 3:1:1
#Pr = 1; Ra = 1174.395870; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 4  , m = 0  , n = 1, aspect ratio 3:1:1
#Pr = 1; Ra = 1202.581371; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # l = 1  , m = 0  , n = 1, aspect ratio 3:1:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:3:1
#Pr = 1; Ra = 660.5178179; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 0  , m = 2  , n = 1, aspect ratio 1:3:1
#Pr = 1; Ra = 779.2727283; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 1|0, m = 0|3, n = 1, aspect ratio 1:3:1
#Pr = 1; Ra = 824.8505622; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 0  , m = 1  , n = 1, aspect ratio 1:3:1
#Pr = 1; Ra = 985.0066489; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 1  , m = 2  , n = 1, aspect ratio 1:3:1
#Pr = 1; Ra = 1174.395870; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 0  , m = 4  , n = 1, aspect ratio 1:3:1
#Pr = 1; Ra = 1202.581371; A1d = 1.0; A2d = 1.0/3.0; A3d = 1.0 # l = 0  , m = 1  , n = 1, aspect ratio 1:3:1
# SF/SF/SF, FF/FF/FT, Aspect ratio 1:1:3
#Pr = 1; Ra = 133.6201523; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 0|1, m = 1|0, n = 1, aspect ratio 1:1:3
#Pr = 1; Ra = 293.5634746; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 0|1, m = 1|0, n = 2, aspect ratio 1:1:3
#Pr = 1; Ra = 458.2503123; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 1  , m = 1  , n = 1, aspect ratio 1:1:3
#Pr = 1; Ra = 711.3936909; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 1  , m = 1  , n = 2, aspect ratio 1:1:3
#Pr = 1; Ra = 779.2727283; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 0|1, m = 1|0, n = 3, aspect ratio 1:1:3
#Pr = 1; Ra = 1315.022729; A1d = 1.0; A2d = 1.0; A3d = 1.0/3.0 # l = 1  , m = 1  , n = 3, aspect ratio 1:1:3

# PAPER
bc_vel = 6 
bc_temp = 4 
heating = 0
Pr = 1; Ra = (216./5.)*np.pi**4; A1d = 1.0; A2d = 1.0; A3d = 1.0 # Paper
#Pr = 7; Ra = (61./9.)*np.pi**4; A1d = 1.0/8.0; A2d = 1.0/8.0; A3d = 1.0 # Paper
#Pr = 1; Ra = 5989.0*np.pi**4; A1d = 8.0; A2d = 8.0; A3d = 1.0 # Paper

# SF/SF/NS, FF/FF/FT
#bc_vel = 4 
#bc_temp = 4 
#heating = 0
# SF/SF/NS, FF/FF/FT, Aspect ratio 3:1:1
#Pr = 1; Ra = 1.5e3; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # Paper
#Pr = 1; Ra = 2e3; A1d = 1.0/3.0; A2d = 1.0; A3d = 1.0 # Paper

# NS/NS/NS, FF/FF/FT
#bc_vel = 0 
#bc_temp = 4
#heating = 0
# NS/NS/NS, FF/FF/FT
#Pr = 1; Ra = 1755.2; A1d = 1.0/6.0; A2d = 1.0/6.0; A3d = 1.0 # Michael Watson's thesis
#Pr = 1; Ra = 1813.0; A1d = 1.0/4.0; A2d = 1.0/4.0; A3d = 1.0 # Michael Watson's thesis
#Pr = 1; Ra = 2084.9; A1d = 1.0/2.0; A2d = 1.0/2.0; A3d = 1.0 # Michael Watson's thesis
#Pr = 1; Ra = 1e4; A1d = 1.0; A2d = 1.0; A3d = 1.0 # Primary bifurcation???

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'heating':0, 'scale1d':2.0*A1d, 'scale2d':2.0*A2d, 'scale3d':2.0*A3d} # Paper
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function that does nothing ...
def wave(k = 0):
    return []

# Wave number function from single "index" (k perpendicular)
eigs = wave()

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['mode'] = 0
marginal_options['ellipse_radius'] = 1e5
marginal_options['geometry'] = 'c3d'
marginal_options['point'] = False
marginal_options['solve'] = True
marginal_options['plot_point'] = True
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
