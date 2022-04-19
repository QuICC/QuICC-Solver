"""Script to run a marginal curve trace for the Boussinesq inviscid inertial wave in a sphere with Worland expansion (Toroidal/Poloidal formulation)"""

import numpy as np

import quicc.model.boussinesq_iiwsphere as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqIIWSphere()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
bc_vel = 0
res = [20, 20, 0]

# Create parameters (rescaling to proper nondimensionalisation)
m = 8 #
res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
eq_params = {'rayleigh':0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(m)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-11
#marginal_options['ellipse_radius'] = 1e5
marginal_options['write_mtx'] = True
marginal_options['geometry'] = 'sphere_worland'
marginal_options['target'] = -0.072j
marginal_options['curve'] = False
marginal_options['minimum'] = False
marginal_options['plot_curve'] = False
marginal_options['solve'] = True
marginal_options['solve_nev'] = 15
marginal_options['point_k'] = m
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['save_hdf5'] = False
marginal_options['viz_mode'] = 0
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
