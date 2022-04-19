"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (streamfunction-temperature)"""

import numpy as np
import scipy.sparse.linalg as splin

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':2340.687}
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.73}
eigs = [1]
bc_str = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

# STREAMFUNCTION FORMULATION

import quicc.model.boussinesq_rotconv2dboxst as mod_st

# Create the model and activate linearization
model = mod_st.BoussinesqRotConv2DBoxST()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':bc_str, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
if False:
    import matplotlib.pylab as pl
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if True:
    import scipy.io as io
    io.mmwrite("matrix_A_st.mtx", A)
    io.mmwrite("matrix_B_st.mtx", B)

# Solve EVP with sptarn
if True:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

    sol_st_s = evp_vec[0:res[0]*res[2],-1]
    sol_st_t = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1]
    sol_st_u = mod_st.c2d.d0d1(res[0], res[2], 0, mod_st.no_bc())*sol_st_s
    sol_st_w = -mod_st.c2d.d1d0(res[0], res[2], 0, mod_st.no_bc())*sol_st_s
    sol_st_c = mod_st.c2d.d1d0(res[0], res[2], 0, mod_st.no_bc())*sol_st_u + mod_st.c2d.d0d1(res[0], res[2], 0, mod_st.no_bc())*sol_st_w

    rhs = (eq_params['rayleigh']/16.)*mod_st.c2d.i2j2d0d1(res[0], res[2], mod_st.no_bc())*sol_st_t
    bc = mod_st.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:22}, 'z':{0:0}, 'priority':'sx'})*sol_st_u + mod_st.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:22}, 'priority':'sx'})*sol_st_w
    poisson = mod_st.c2d.i2j2lapl(res[0], res[2], 0, {'x':{0:21},'z':{0:21}, 'priority':'sx'})
    poisson[0,:] = 0
    poisson[0,0] = 1
    sol_st_p = splin.spsolve(poisson,rhs + bc)

    sol_st_s = sol_st_s/np.max(np.abs(sol_st_s))
    sol_st_t = sol_st_t/np.max(np.abs(sol_st_t))
    sol_st_u = sol_st_u/np.max(np.abs(sol_st_u))
    sol_st_w = sol_st_w/np.max(np.abs(sol_st_w))
    sol_st_p = sol_st_p/np.max(np.abs(sol_st_p))

    mat_st_s = sol_st_s.reshape(res[0], res[2], order = 'F')
    mat_st_t = sol_st_t.reshape(res[0], res[2], order = 'F')
    mat_st_u = sol_st_u.reshape(res[0], res[2], order = 'F')
    mat_st_w = sol_st_w.reshape(res[0], res[2], order = 'F')
    mat_st_c = sol_st_c.reshape(res[0], res[2], order = 'F')
    mat_st_p = sol_st_p.reshape(res[0], res[2], order = 'F')

# VELOCITY_PRESSURE

import quicc.model.boussinesq_rotconv2dboxvpt as mod_vpt

# Create the model and activate linearization
model = mod_vpt.BoussinesqRotConv2DBoxVPT()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityz':bc_vel, 'pressure':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
A = A.tolil()
## TOP CLEARED OUT
#A[3*res[0]*res[2],:] = 0
#A[3*res[0]*res[2],3*res[0]*res[2]] = 1
#A[3*res[0]*res[2]+1,:] = 0
#A[3*res[0]*res[2]+1,3*res[0]*res[2]+res[0]-1] = 1
#
#A[3*res[0]*res[2]+res[0],:] = 0
#A[3*res[0]*res[2]+res[0],-res[0]] = 1
#A[3*res[0]*res[2]+res[0]+1,:] = 0
#A[3*res[0]*res[2]+res[0]+1,-1] = 1
#
#A[-(res[0]+2),:] = 0
#A[-(res[0]+2),-(res[0]+2)] = 1

# BOTTOM CLEARED OUT
A[-2*res[0],:] = 0
A[-2*res[0],3*res[0]*res[2]] = 1
A[-2*res[0]+1,:] = 0
A[-2*res[0]+1,3*res[0]*res[2]+res[0]-1] = 1

A[-res[0],:] = 0
A[-res[0],-res[0]] = 1
A[-res[0]+1,:] = 0
A[-res[0]+1,-1] = 1

A[-(res[0]+2),:] = 0
A[-(res[0]+2),4*res[0]*res[2] - res[0] - 2] = 1

A = A.tocsr()

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
if False:
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if True:
    io.mmwrite("matrix_A_vpt.mtx", A)
    io.mmwrite("matrix_B_vpt.mtx", B)

# Solve EVP with sptarn
if True:
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    evp_sel = evp_vec[:,-1]
    print(evp_lmb)

    sol_vpt_u = evp_sel[0:res[0]*res[2]]
    sol_vpt_w = evp_sel[res[0]*res[2]:2*res[0]*res[2]]
    sol_vpt_t = evp_sel[2*res[0]*res[2]:3*res[0]*res[2]]
    sol_vpt_p = evp_sel[3*res[0]*res[2]:]
    sol_vpt_c = mod_vpt.c2d.d1d0(res[0], res[2], 0, mod_vpt.no_bc())*sol_vpt_u + mod_vpt.c2d.d0d1(res[0], res[2], 0, mod_vpt.no_bc())*sol_vpt_w

    sol_vpt_u = sol_vpt_u/np.max(np.abs(sol_vpt_u))
    sol_vpt_t = sol_vpt_t/np.max(np.abs(sol_vpt_t))
    sol_vpt_w = sol_vpt_w/np.max(np.abs(sol_vpt_w))
    sol_vpt_p = sol_vpt_p/np.max(np.abs(sol_vpt_p))
        
    mat_vpt_u = sol_vpt_u.reshape(res[0], res[2], order = 'F')
    mat_vpt_w = sol_vpt_w.reshape(res[0], res[2], order = 'F')
    mat_vpt_t = sol_vpt_t.reshape(res[0], res[2], order = 'F')
    mat_vpt_p = sol_vpt_p.reshape(res[0], res[2], order = 'F')
    mat_vpt_c = sol_vpt_c.reshape(res[0], res[2], order = 'F')



    import matplotlib.pylab as pl
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_st_u - mat_vpt_u)))
    pl.colorbar()
    pl.title('u')
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(mat_st_w - mat_vpt_w)))
    pl.colorbar()
    pl.title('w')
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(mat_st_t - mat_vpt_t)))
    pl.colorbar()
    pl.title('T')
    pl.subplot(2,3,4)
    pl.imshow(np.log10(np.abs(mat_st_s)))
    pl.colorbar()
    pl.title('Streamfunction')
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(mat_st_c - mat_vpt_c)))
    pl.colorbar()
    pl.title('Continuity')
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(mat_st_p - mat_vpt_p)))
    pl.colorbar()
    pl.title('p')
    pl.show()
    pl.close("all")

    import quicc.transform.cartesian as transf
    phys_s = transf.tophys2d(mat_st_s)
    phys_t = transf.tophys2d(mat_st_t - mat_vpt_t)
    phys_u = transf.tophys2d(mat_st_u - mat_vpt_u)
    phys_w = transf.tophys2d(mat_st_w - mat_vpt_w)
    phys_c = transf.tophys2d(mat_st_c - mat_vpt_c)
    phys_p = transf.tophys2d(mat_st_p - mat_vpt_p)
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_u)), 50)
    pl.colorbar()
    pl.title('u')
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_w)), 50)
    pl.colorbar()
    pl.title('w')
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_t)), 50)
    pl.colorbar()
    pl.title('T')
    pl.subplot(2,3,4)
    pl.contourf(grid_x, grid_z, phys_s, 50)
    pl.colorbar()
    pl.title('Streamfunction')
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title('Continuity')
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p, 50)
    pl.colorbar()
    pl.title('p')
    pl.show()
    pl.close("all")
