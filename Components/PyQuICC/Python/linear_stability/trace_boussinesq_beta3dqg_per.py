"""Script to run a marginal curve trace for the periodic Boussinesq Beta 3DQG model"""

import numpy as np
import functools

import quicc.model.boussinesq_beta3dqg_per as mod
import quicc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqBeta3DQGPer()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]
chi = 5
Pr = 1
G = 1e-2
tanchi = np.tan(chi*np.pi/180.)
beta = tanchi/G

kxc = 0.0
kyc = (2.0**(-1.0/6.0)*(beta*Pr/(1.0 + Pr))**(1.0/3.0))
Rac = 16*((kxc**2 + kyc**2)**3/kyc**2 + beta**2*Pr**2/((kxc**2 + kyc**2)*(Pr + 1.0)**2))
fc = -beta*kyc/((kxc**2 + kyc**2)*(1.0 + Pr))
print("Ra_c = " + str(Rac))
print("k_c = " + str(kyc))
print("f_c = " + str(fc))

# Set wave number
kp = 14.4046937485972/2.0
kp = 3/2.0
Ra = 192407.5882
Ra = 290
phi = 90
Ra = 500
kp = 1
Ra = Rac
kp = kyc

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'gamma':G, 'chi':chi, 'scale1d':1, 'elevator':1}
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}

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
    #kx = 1
    #ky = kp
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(kp)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['ellipse_radius'] = 1e3
marginal_options['mode'] = 0
marginal_options['point'] = False
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
marginal_options['viz_mode'] = 0
marginal_options['curve_points'] = np.arange(max(0.3, kp-2), kp+.5, 0.05)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
