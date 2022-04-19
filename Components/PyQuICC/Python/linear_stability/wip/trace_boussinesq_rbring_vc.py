"""Script to run a marginal curve trace for the 2D Boussinesq Rayleigh-Benard convection in a ring (velocity-continuity formulation)"""

import numpy as np

import quicc.model.boussinesq_rbring_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRBRingVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [500, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':20412, 'ro':1, 'rratio':0.35}
eigs = [4]
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 2 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import quicc.transform.annulus as transf

# Show the "spy" of the two matrices
if show_spy:
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    import quicc.linear_stability.solver as solver
    #evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e3, -1e2)
    print(evp_lmb)

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0],mode]
        sol_v = evp_vec[res[0]:2*res[0],mode]
        # Extract continuity from velocity 
        a, b = mod.annulus.linear_r2x(eq_params['ro'], eq_params['rratio'])
        sol_c = mod.annulus.x1div(res[0], a, b, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
        intg_c = mod.annulus.i1x1div(res[0], a, b, mod.no_bc())*sol_u + mod.annulus.i1(res[0], a, b, mod.no_bc(), 1j*eigs[0])*sol_v
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))) + ", Max integrated continuity: " + str(np.max(np.abs(intg_c))))

if show_solution:
    viz_mode = 1
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0],viz_mode]
    sol_v = evp_vec[res[0]:2*res[0],viz_mode]
    sol_t = evp_vec[2*res[0]:3*res[0],viz_mode]
    sol_p = evp_vec[3*res[0]:4*res[0],viz_mode]
    # Extract continuity from velocity 
    a, b = mod.annulus.linear_r2x(eq_params['ro'], eq_params['rratio'])
    sol_c = mod.annulus.x1div(res[0], a, b, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
    intg_c = mod.annulus.i1x1div(res[0], a, b, mod.no_bc())*sol_u + mod.annulus.i1(res[0], a, b, mod.no_bc(), 1j*eigs[0])*sol_v
    
    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_u))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.semilogy(np.abs(sol_v))
    pl.title("v")
    pl.subplot(2,3,3)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.subplot(2,3,4)
    pl.semilogy(np.abs(sol_c))
    pl.title("Continuity")
    pl.subplot(2,3,5)
    pl.semilogy(np.abs(intg_c))
    pl.title("Integrated continuity")
    pl.subplot(2,3,6)
    pl.semilogy(np.abs(sol_p))
    pl.title("p")
    pl.show()
    pl.close("all")

    # Compute physical space values
    grid_r = transf.rgrid(res[0], a, b)
    phys_u = transf.torphys(sol_u.real)
    phys_v = transf.torphys(sol_v.real)
    phys_t = transf.torphys(sol_t.real)
    phys_p = transf.torphys(sol_p.real)
    phys_c = transf.torphys(sol_c.real)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.plot(grid_r, phys_u)
    pl.title("u")
    pl.subplot(2,3,2)
    pl.plot(grid_r, phys_v)
    pl.title("v")
#    pl.subplot(2,3,3)
#    pl.plot(grid_r, phys_w)
#    pl.title("w")
    pl.subplot(2,3,3)
    pl.plot(grid_r, phys_t)
    pl.title("T")
    pl.subplot(2,3,5)
    pl.plot(grid_r, np.log10(np.abs(phys_c)))
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.plot(grid_r, phys_p)
    pl.title("p")
    pl.show()
    pl.close("all")
