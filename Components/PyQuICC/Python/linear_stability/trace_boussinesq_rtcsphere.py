"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a sphere with Worland expansion (Toroidal/Poloidal formulation)"""

import numpy as np

import quicc.model.boussinesq_rtcsphere as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCSphere()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
Rac = None
mc = None

# SF, FT, internal heating
bc_vel = 1; bc_temp = 0
#E = (1.0/1e6)**0.5
#res = [32, 32, 0]
#E = (1.0/1e7)**0.5
#res = [32, 32, 0]
#E = (1.0/1e8)**0.5
#res = [32, 32, 0]
#E = (1.0/1e9)**0.5
#res = [48, 48, 0]
E = (1.0/1e10)**0.5
res = [64, 64, 0]
#E = (1.0/1e11)**0.5
#res = [96, 96, 0]
#E = (1.0/1e12)**0.5
#res = [256, 128, 0]
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

# NS, FT, internal heating
#bc_vel = 0; bc_temp = 0
#E = (1.0/1e6)**0.5
#res = [32, 32, 0]
#E = (1.0/1e7)**0.5
#res = [32, 32, 0]
#E = (1.0/1e8)**0.5
#res = [32, 32, 0]
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

# Create parameters (rescaling to proper nondimensionalisation)
if mc is None:
    m = np.int(0.3029*E**(-2./6.)) # Asymptotic prediction for minimum
else:
    m = mc
if Rac is None:
    Ra = (4.1173*E**(-4./3.) + 17.7815*E**(-2./2.))*E**(2./2.) # Asymptotic prediction for critical Rayleigh number
else:
    Ra = Rac

res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
eq_params = {'ekman':E, 'prandtl':1, 'rayleigh':Ra}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(1)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-16
marginal_options['geometry'] = 'sphere_worland'
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = False
marginal_options['solve'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
