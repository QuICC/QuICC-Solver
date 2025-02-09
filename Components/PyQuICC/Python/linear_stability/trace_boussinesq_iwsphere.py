"""Script to run a marginal curve trace for the Boussinesq inertial wave in a sphere with Worland expansion (Toroidal/Poloidal formulation)"""

import numpy as np

import quicc.model.boussinesq_iwsphere as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqIWSphere()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
bc_vel = 1 # 0: NS, 1: SF
#E = (1.0/1e6)**0.5
#res = [32, 32, 0]
#E = (1.0/1e7)**0.5
#res = [32, 32, 0]
E = (1.0/1e8)**0.5
res = [48, 48, 0]
#E = (1.0/1e9)**0.5
#res = [48, 48, 0]
#E = (1.0/1e10)**0.5
#res = [64, 64, 0]
#E = (1.0/1e11)**0.5
#res = [96, 96, 0]
#E = (1.0/1e12)**0.5
#res = [128, 128, 0]
#E = (1.0/1e13)**0.5
#res = [192, 192, 0]
#E = (1.0/1e14)**0.5
#res = [256, 256, 0]
#E = (1.0/1e15)**0.5
#res = [512, 512, 0]
#E = (1.0/1e16)**0.5
#res = [512, 768, 0]
#E = (1.0/1e17)**0.5
#res = [512, 1024, 0]
#E = (1.0/1e18)**0.5
#res = [784, 1536, 0]
#E = (1.0/1e19)**0.5
#res = [512, 512, 0]

# Create parameters (rescaling to proper nondimensionalisation)
m = 1 #
res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
eq_params = {'ekman':E, 'rayleigh':0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(1)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-11
#marginal_options['ellipse_radius'] = 1e5
marginal_options['geometry'] = 'sphere_worland'
marginal_options['target'] = 0.6547j
marginal_options['curve'] = False
marginal_options['minimum'] = False
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['solve_nev'] = 15
marginal_options['point_k'] = m
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = False
marginal_options['save_hdf5'] = True
marginal_options['viz_mode'] = 1
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
