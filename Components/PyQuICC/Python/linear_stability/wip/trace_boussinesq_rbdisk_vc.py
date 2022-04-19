"""Script to run a marginal curve trace for the 2D Boussinesq Rayleigh-Benard convection in a disk (velocity-continuity formulation)"""

import numpy as np

#import quicc.model.boussinesq_rbdisk_vc as mod
#import quicc.model.boussinesq_rbdisk_vc2 as mod
import quicc.model.boussinesq_rbdisk_vc3 as mod

# Create the model and activate linearization
model = mod.BoussinesqRBDiskVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [500, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':5901.55}
#eq_params = {'prandtl':1, 'rayleigh':0.}
eigs = [7]
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

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
    import quicc.transform.cylinder as transf

# Show the "spy" of the two matrices
if show_spy:
    import matplotlib.pylab as pl
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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e3, -1e2)
    print(evp_lmb)

if show_solution:
    viz_mode = 1

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0],mode]
        sol_v = evp_vec[res[0]:2*res[0],mode]
        # Extract continuity from velocity 
        #sol_c = mod.cylinder.x1div(res[0], (eigs[0]+1)%2, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
        #intg_c = mod.cylinder.i1x1div(res[0], (eigs[0]+1)%2, mod.no_bc())*sol_u + mod.cylinder.i1(res[0], (eigs[0]+1)%2 , mod.no_bc(), 1j*eigs[0])*sol_v
        sol_c = mod.cylinder.x1d1(res[0], eigs[0]%2, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
        intg_c = mod.cylinder.i1x1d1(res[0], eigs[0]%2, mod.no_bc())*sol_u + mod.cylinder.i1(res[0], eigs[0]%2 , mod.no_bc(), 1j*eigs[0])*sol_v
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))) + ", Max integrated continuity: " + str(np.max(np.abs(intg_c))))

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0],viz_mode]
    sol_v = evp_vec[res[0]:2*res[0],viz_mode]
    sol_t = evp_vec[2*res[0]:3*res[0],viz_mode]
    sol_p = evp_vec[3*res[0]:4*res[0],viz_mode]
    # Extract continuity from velocity 
    #sol_c = mod.cylinder.x1div(res[0], (eigs[0]+1)%2, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
    #intg_c = mod.cylinder.i1x1div(res[0], (eigs[0]+1)%2, mod.no_bc())*sol_u + mod.cylinder.i1(res[0], (eigs[0]+1)%2 , mod.no_bc(), 1j*eigs[0])*sol_v
    sol_c = mod.cylinder.x1d1(res[0], eigs[0]%2, mod.no_bc(), zr = 0)*sol_u + 1j*eigs[0]*sol_v
    intg_c = mod.cylinder.i1x1d1(res[0], eigs[0]%2, mod.no_bc())*sol_u + mod.cylinder.i1(res[0], eigs[0]%2 , mod.no_bc(), 1j*eigs[0])*sol_v
    
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
    pl.subplot(2,3,5)
    pl.semilogy(np.abs(sol_c))
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.semilogy(np.abs(sol_p))
    pl.title("p")
    pl.show()
    pl.close("all")

    # Compute physical space values
    grid_r = transf.rgrid(res[0])
    phys_u = transf.torphys(sol_u.real, eigs[0]%2)
    phys_v = transf.torphys(sol_v.real, eigs[0]%2)
    phys_t = transf.torphys(sol_t.real, eigs[0]%2)
    phys_p = transf.torphys(sol_p.real, eigs[0]%2)
    phys_c = transf.torphys(sol_c.real, eigs[0]%2)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.plot(grid_r, phys_u)
    pl.title("u")
    pl.subplot(2,3,2)
    pl.plot(grid_r, phys_v)
    pl.title("v")
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
