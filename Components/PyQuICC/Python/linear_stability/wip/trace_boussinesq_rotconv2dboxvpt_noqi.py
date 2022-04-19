"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (velocity-pressure-temperature) without quasi-inverse"""

import numpy as np

import quicc.model.boussinesq_rotconv2dboxvpt_noqi as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConv2DBoxVPTNoQI()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]
pN = 0
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':2340.687}
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.73}
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':0.0}
eigs = [1]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

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
A[-res[0]-2,:] = 0
A[:,-res[0]-2] = 0
A[-res[0]-2,-res[0]-2] = 1

A[-res[0]-1,:] = 0
A[:,-res[0]-1] = 0
A[-res[0]-1,-res[0]-1] = 1

A[-2,:] = 0
A[:,-2] = 0
A[-2,-2] = 1

A[-1,:] = 0
A[:,-1] = 0
A[-1,-1] = 1

A[3*res[0]*res[2],:] = 0
A[:,3*res[0]*res[2]] = 0
A[3*res[0]*res[2],3*res[0]*res[2]] = 1

A = A.tocsr()

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
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if True:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, np.inf)
    evp_sel = evp_vec[:,-1]
    print(evp_lmb)

    sol_u = evp_sel[0:res[0]*res[2]]
    print('U value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:20}, 'z':{0:0}, 'priority':'x'})*sol_u
    print(t[np.nonzero(t)])
    print('U value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:20}, 'priority':'z'})*sol_u
    print(t[np.nonzero(t)])
    print('Dx U value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:21}, 'z':{0:0}, 'priority':'x'})*sol_u
    print(t[np.nonzero(t)])
    print('Dz U value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:21}, 'priority':'z'})*sol_u
    print(t[np.nonzero(t)])
    sol_w = evp_sel[res[0]*res[2]:2*res[0]*res[2]]
    print('W value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:20}, 'z':{0:0}, 'priority':'x'})*sol_w
    print(t[np.nonzero(t)])
    print('W value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:20}, 'priority':'z'})*sol_w
    print(t[np.nonzero(t)])
    print('Dx W value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:21}, 'z':{0:0}, 'priority':'x'})*sol_w
    print(t[np.nonzero(t)])
    print('Dz W value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:21}, 'priority':'z'})*sol_w
    print(t[np.nonzero(t)])
    sol_t = evp_sel[2*res[0]*res[2]:3*res[0]*res[2]]
    sol_p = evp_sel[3*res[0]*res[2]:]
    sol_c = mod.c2d.d1d0(res[0], res[2], mod.no_bc(), sz = 0)*sol_u + mod.c2d.d0d1(res[0], res[2], mod.no_bc(), sx = 0)*sol_w
    print('C value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:20}, 'z':{0:0}, 'priority':'x'})*sol_c
    print(t[np.nonzero(t)])
    print('C value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:20}, 'priority':'z'})*sol_c
    print(t[np.nonzero(t)])
    sol_cu = mod.c2d.d1d0(res[0], res[2], mod.no_bc(), sz = 0)*sol_u
    sol_cw = mod.c2d.d0d1(res[0], res[2], mod.no_bc(), sx = 0)*sol_w

    print('T value on X boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:20}, 'z':{0:0}, 'priority':'x'})*sol_t
    print(t[np.nonzero(t)])
    print('T value on Z boundary')
    t = mod.c2d.zblk(res[0], res[0], 2, 2, {'x':{0:0}, 'z':{0:20}, 'priority':'z'})*sol_t
    print(t[np.nonzero(t)])

#    sol_u = sol_u/np.max(sol_u)
#    sol_t = sol_t/np.max(sol_t)
#    sol_w = sol_w/np.max(sol_w)
#    sol_p = sol_p/np.max(sol_p)

    mat_u = sol_u.reshape(res[0], res[2], order = 'F')
    mat_w = sol_w.reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_p = sol_p.reshape(res[0]-pN, res[2]-pN, order = 'F')
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')
    mat_cu = sol_cu.reshape(res[0], res[2], order = 'F')
    mat_cw = sol_cw.reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.semilogy(np.abs(sol_c))
    pl.show()

    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(mat_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(mat_t)))
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(mat_c)))
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(mat_p)))
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    import quicc.transform.cartesian as transf
    phys_u = transf.tophys2d(mat_u)
    phys_w = transf.tophys2d(mat_w)
    phys_t = transf.tophys2d(mat_t)
    phys_p = transf.tophys2d(mat_p)
    phys_c = transf.tophys2d(mat_c)
    phys_cu = transf.tophys2d(mat_cu)
    phys_cw = transf.tophys2d(mat_cw)
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(2,5,1)
    pl.semilogy(grid_x, np.abs(phys_c[res[2]//2,:]))
    pl.title('d_x u + d_z w')
    pl.subplot(2,5,2)
    pl.semilogy(grid_x, np.abs(phys_cu[res[2]//2,:]))
    pl.title('d_x u')
    pl.subplot(2,5,3)
    pl.semilogy(grid_x, np.abs(phys_cw[res[2]//2,:]))
    pl.title('d_z w')
    pl.subplot(2,5,4)
    pl.semilogy(grid_x, np.abs(phys_u[res[2]//2,:]))
    print((phys_u[res[2]//2,0], phys_u[res[2]//2,-1]))
    pl.title('u')
    pl.subplot(2,5,5)
    pl.semilogy(grid_x, np.abs(phys_w[res[2]//2,:]))
    print((phys_w[res[2]//2,0], phys_w[res[2]//2,-1]))
    pl.title('w')
    pl.subplot(2,5,6)
    pl.semilogy(grid_x, np.abs(phys_c[:,res[2]//2]))
    pl.title('d_x u + d_z w')
    pl.subplot(2,5,7)
    pl.semilogy(grid_x, np.abs(phys_cu[:,res[2]//2]))
    pl.title('d_x u')
    pl.subplot(2,5,8)
    pl.semilogy(grid_x, np.abs(phys_cw[:,res[2]//2]))
    pl.title('d_z w')
    pl.subplot(2,5,9)
    pl.semilogy(grid_x, np.abs(phys_u[:,res[2]//2]))
    print((phys_u[0,res[2]//2], phys_u[-1,res[2]//2]))
    pl.title('u')
    pl.subplot(2,5,10)
    pl.semilogy(grid_x, np.abs(phys_w[:,res[2]//2]))
    print((phys_w[0,res[2]//2], phys_w[-1,res[2]//2]))
    pl.title('w')
    pl.show()

    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    #pl.contourf(grid_x, grid_z, phys_c, 50)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
