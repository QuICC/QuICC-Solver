"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a cylinder model with Worland expansion (Toroidal/Poloidal formulation)"""

import numpy as np

import quicc.model.boussinesq_rbccylinder_v2 as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRBCCylinder()
model.linearize = True
model.use_galerkin = False

# Set boundary conditions
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 1 # 0: FT/FT, 2: FF/FT, 3: FT/FF

# Create parameters
m = 2
res = [48, 0, 48]
eq_params = {'prandtl':1, 'rayleigh':2.66491e3, 'gamma':1.0, 'scale3d':2.0}
#eq_params = {'prandtl':1, 'rayleigh':1193, 'gamma':1.0, 'scale3d':2.0}
eq_params = {'prandtl':1, 'rayleigh':2000, 'gamma':1.0, 'scale3d':2.0}
auto_params = model.automatic_parameters(eq_params)
for k,v in auto_params.items():
    eq_params[k] = v
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-12
marginal_options['geometry'] = 'cylinder_worland'
#marginal_options['ellipse_radius'] = 1e3
#marginal_options['target'] = 4e4
marginal_options['curve'] = False
marginal_options['minimum'] = False
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = False
marginal_options['solve'] = True
marginal_options['solve_nev'] = 3
marginal_options['point_k'] = m
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = True
marginal_options['viz_mode'] = -1
marginal_options['show_spectra'] = True
marginal_options['save_spectra'] = False
marginal_options['show_physical'] = True
marginal_options['save_physical'] = False
marginal_options['write_mtx'] = True
marginal_options['save_pdf'] = False
marginal_options['curve_points'] = np.arange(m, m+1, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
