import numpy as np
import quicc.geometry.spherical.sphere_radius_worland as geo
import scipy.sparse as sp
import scipy.sparse.linalg as SLA
import scipy.linalg as DLA
import scipy.io as io

def evalG(u, a, b, f, t):
    if f is None:
        g = np.zeros(u.shape)
    else:
        n = f.shape[0]
        g = np.zeros(n)
        for i in range(0, n):
            g[i] = (f[i,0]*np.cos(a*t) + f[i,1]*np.sin(a*t) + f[i,2]*np.cos(b*t) + f[i,3]*np.sin(b*t))
        #    g[i] = 2*u[i]

    return g

def tauOps(ops):

    A, B, Bc, Q, Qc, C = ops

    # Count bC
    nbc = 1
    if np.max(np.abs(C.tolil()[1,:])) != 0:
        nbc += 1
    
    # Get linear operator
    tQ = Q.tolil()
    tQ[0,-1] = 1
    if nbc == 2:
        tQ[1,-2] = 1
    tQ = tQ.tocsc()
    L = SLA.spsolve(tQ, A.tocsc())
    # Get mass matrix
    T = SLA.spsolve(tQ, Bc.tocsc())
    # Get "QI"
    Qc = geo.tid(nN, l, nbc, {0:0})
    Id = geo.tid(nN, l, 0, {0:0})

    # Move boundary condition to bottom
    C = C.tolil()
    C[-1,:] = C[0,:]
    C[-2,:] = C[1,:]
    C[0:2,:] = 0
    C = C.tocoo()

    return (L, Id, T, Id, Qc, C)


def pc2Step(t, dt, u0, nN, a, b, nl, ops, ncorr = 1):

    #A, B, _, Q, _, C = ops
    A, _, B, _, Q, C = ops

    c = 1./2.
    LHS = B - dt*c*A + C
    RHS = B + dt*c*A

    g0 = evalG(u0, a, b, nl, t)
    u = SLA.spsolve(LHS, RHS*u0 + dt*Q*g0)
    for i in range(0,ncorr):
        g1 = evalG(u, a, b, nl, t+dt)
        u += SLA.spsolve(LHS, c*dt*(Q*g1 - Q*g0))

    return u


def pc2TauStep(t, dt, u0, nN, a, b, nl, ops):

    ops_tau = tauOps(ops)

    return pc2Step(t, dt, u0, nN, a, b, nl, ops_tau)

def cnrkw3():
    # stages
    s = 4 

    # DIRK
    aIm = np.zeros((s,s))
    aIm[1,0] = 4./15.; aIm[1,1] = 4./15.
    aIm[2,0] = 4./15.; aIm[2,1] = 1./3.; aIm[2,2] = 1./15.
    aIm[3,0] = 4./15.; aIm[3,1] = 1./3.; aIm[3,2] = 7./30.; aIm[3,3] = 1./6.
    bIm = np.zeros(s)
    bIm[0] = 4./15.; bIm[1] = 1./3.; bIm[2] = 7./30.; bIm[3] = 1./6.
    cIm = np.zeros(s)
    cIm[1] = 8./15.; cIm[2] = 2./3.; cIm[3] = 1.

    # ERK
    aEx = np.zeros((s,s))
    aEx[1,0] = 8./15.
    aEx[2,0] = 1./4.; aEx[2,1] = 5./12.
    aEx[3,0] = 1./4.; aEx[3,2] = 3./4.
    bEx = np.zeros(s)
    bEx[0] = 1./4.; bEx[2] = 3./4.
    cEx = cIm

    return (aIm, bIm, cIm, aEx, bEx, cEx)

def cb2():
    # stages
    s = 3

    # DIRK
    aIm = np.zeros((s,s))
    aIm[1,1] = 2./5.;
    aIm[2,1] = 5./6.; aIm[2,2] = 1./6.
    bIm = np.zeros(s)
    bIm[1] = 5./6.; bIm[2] = 1./6.
    cIm = np.zeros(s)
    cIm[1] = 2./5.; cIm[2] = 1.

    # ERK
    aEx = np.zeros((s,s))
    aEx[1,0] = 2./5.
    aEx[2,1] = 1.
    bEx = bIm
    cEx = cIm

    return (aIm, bIm, cIm, aEx, bEx, cEx)

def rkLinStep(t, dt, u0, nN, a, b, nl, ops, scheme):

    A, B, Bc, Q, Qc, C = ops

    aIm, bIm, cEx, _, _, _ = scheme()
    ns = aIm.shape[0]

    x = u0

    for s in range(0,ns):
        if s == 0:
            y = x
            #x = Bc*x
        else:
            #y = x + (aIm[s,s-1] - bIm[s-1])*dt*Bc*z
            y = SLA.spsolve(LHS, RHS*x + (aIm[s,s-1] - bIm[s-1])*dt*A*y)
            #y = SLA.spsolve(Bc + C, y)
        LHS = Bc - aIm[s,s]*dt*A + C
        RHS = Bc - aIm[s,s]*dt*A
        #z = SLA.spsolve(LHS, A*y)
        #x = x + bIm[s]*dt*Bc*z
        x = SLA.spsolve(LHS, RHS*x + bIm[s]*dt*A*y)

    #u = SLA.spsolve(Bc+C, x)
    u = x

    return u

def rkStep(t, dt, u0, nN, a, b, nl, ops, scheme):

    A, _, Bc, _, Qc, C = ops

    dense = False
    if dense:
        A = A.todense()
        Bc = Bc.todense()
        Qc = Qc.todense()
        C = C.todense()

    aIm, bIm, cEx, aEx, bEx, cEx = scheme()
    ns = aIm.shape[0]

    if dense:
        x = u0.reshape((u0.shape[0],1))
    else:
        x = u0

    for s in range(0,ns):
        if s == 0:
            y = x
            x = Bc*x
        else:
            y = x + (aIm[s,s-1] - bIm[s-1])*dt*Bc*z + (aEx[s,s-1] - bEx[s-1])*dt*y
            if dense:
                y = DLA.solve(Bc + C, y)
            else:
                y = SLA.spsolve(Bc + C, y)
        LHS = Bc - aIm[s,s]*dt*A + C
        if dense:
            z = DLA.solve(LHS, A*y)
        else:
            z = SLA.spsolve(LHS, A*y)
        y = Qc*evalG(y + aIm[s,s]*dt*z, a, b, nl, t + cEx[s]*dt)
        x = x + bIm[s]*dt*Bc*z + bEx[s]*dt*y

    u = SLA.spsolve(Bc+C, x)

    return u

def checkSolution(u, ref, tag, showSol = False):
    if showSol:
        print(f'Solution {tag}: {u}')
    err = np.max(np.abs(u - ref))
    print(f'{tag} error: {err}')

    return err
#
# 2nd order system
#
print(100*'-')
print(100*'=')
print(100*'-')
for id in range(0,7):

    fbase = f'bessel_id{id}'
    basedir = '../_refdata/Framework/Timestep/Worland/'
    meta = np.genfromtxt(basedir + f'{fbase}_meta.dat')
    nN = int(meta[1])
    l = meta[2]
    a,b = meta[3:5]
    dts = 10**meta[5:]

    u0 = np.genfromtxt(basedir + f'{fbase}_in.dat')
    try:
        nl = np.genfromtxt(basedir + f'{fbase}_forcing.dat')
    except:
        nl = None
    uRef = np.genfromtxt(basedir + f'{fbase}_ref.dat')

    A = geo.i2lapl(nN, l, {0:0})
    B = geo.i2(nN, l, {0:0})
    Bcorr = B - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
    Q = geo.i2(nN, l, {0:0})
    Qcorr = Q - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
    C = geo.zblk(nN, l, {0:10})
    ops = (A, B, Bcorr, Q, Qcorr, C)

    #import scipy.io as io
    #AA = io.mmread('I2Lapl_id99_ref.dat')

    print(f'BC Error u0: {np.max(np.abs(C*u0))}')

    showSol = True
    nt = dts.shape[0]
    errors = np.zeros((nt, 8))
    for idt in range(0, nt):
        t = 0
        dt = dts[idt]
        print(f'nN = {nN}')
        print(f'l = {l}')
        print(f'a = {a}')
        print(f'b = {b}')
        print(f'dt = {dt}')

        if showSol:
            print(f'Reference: {uRef[:,idt]}')

        errors[idt, 0] = dt

        u_pc2 = pc2Step(t, dt, u0, nN, a, b, nl, ops)
        errors[idt, 1] = checkSolution(u_pc2, uRef[:,idt], 'PC2', showSol)

        u_pc2_tau = pc2TauStep(t, dt, u0, nN, a, b, nl, ops)
        errors[idt, 2] = checkSolution(u_pc2_tau, uRef[:,idt], 'Tau PC2', showSol)

        u_rk_cb2 = rkStep(t, dt, u0, nN, a, b, nl, ops, cb2)
        errors[idt, 3] = checkSolution(u_rk_cb2, uRef[:,idt], 'RK CB2', showSol)

        u_rk_cnrkw3 = rkStep(t, dt, u0, nN, a, b, nl, ops, cnrkw3)
        errors[idt, 4] = checkSolution(u_rk_cnrkw3, uRef[:,idt], 'RK CN/RKW3', showSol)

        print(100*'+')
        print(f'Diff PC2 vs Ref: {np.abs(u_pc2 - uRef[:,idt])}')
        print(f'Diff PC2 vs RK CB2: {np.abs(u_pc2 - u_rk_cb2)}')
        print(f'Diff PC2 vs RK CN/RKW3: {np.abs(u_pc2 - u_rk_cnrkw3)}')

    import matplotlib.pylab as pl

    pl.title(f' l = {l}')
    pl.xscale('log')
    pl.yscale('log')
    pl.plot(errors[:,0], errors[:,1], label='PC2')
    pl.plot(errors[:,0], errors[:,2], label='Tau PC2')
    pl.plot(errors[:,0], errors[:,3], label='RK CB2')
    pl.plot(errors[:,0], errors[:,4], label='RK CN/RKW3')
    pl.plot(errors[:,0], errors[:,0], 'k--', label='X')
    pl.plot(errors[:,0], errors[:,0]**2, 'k-.', label='X^2')
    pl.plot(errors[:,0], errors[:,0]**3, 'k:', label='X^3')
    pl.legend()
    pl.show()

##
## 4th order system
##
#print(100*'-')
#print(100*'=')
#print(100*'-')
#id = 6
#
#fbase = f'bidiffusion_id{id}'
#basedir = '../_refdata/Framework/Timestep/Worland/'
#meta = np.genfromtxt(basedir + f'{fbase}_meta.dat')
#nN = int(meta[1])
#l = meta[2]
#a,b = meta[3:5]
#dts = 10**meta[5:]
#
#u0 = np.genfromtxt(basedir + f'{fbase}_in.dat')
#nl = np.genfromtxt(basedir + f'{fbase}_forcing.dat')
#uRef = np.genfromtxt(basedir + f'{fbase}_ref.dat')
#
#showSol = False
#idt = 4
#t = 0
#dt = dts[idt]
#print(f'nN = {nN}')
#print(f'l = {l}')
#print(f'a = {a}')
#print(f'b = {b}')
#print(f'dt = {dt}')
#
#A = geo.i4lapl2(nN, l, {0:0})
#corr_mat = sp.lil_matrix((nN,nN))
#corr_mat[-2,-1] = geo.i4lapl(nN+1, l, {0:0}).tolil()[-1,-2]/geo.i4(nN+1, l, {0:0}).tolil()[-1,-3]
#B = geo.i4lapl(nN, l, {0:0})
#Bcorr = B - geo.i4(nN, l, {0:0})*corr_mat
#Q = geo.i4(nN, l, {0:0})
#Qcorr = Q - geo.i4(nN, l, {0:0})*geo.qid(nN, l, nN-2, {0:0})
#C = geo.zblk(nN, l, {0:21})
#ops = (A, B, Bcorr, Q, Qcorr, C)
#
#print(f'BC Error u0: {np.max(np.abs(C*u0))}')
#
#u_pc2 = pc2Step(t, dt, u0, nN, a, b, nl, ops)
#u = u_pc2
#if showSol:
#    print(f'Solution 4th order PC2: {u}')
#err = np.max(np.abs(u - uRef[:,idt]))
#print(f'PC2 error for dt = {dt}: {err}')
#
#u_cb2 = rkcb2Step(t, dt, u0, nN, a, b, nl, ops)
#u = u_cb2
#if showSol:
#    print(f'Solution 4th order RKCB2: {u}')
##print(f'Diff: {np.abs(u - uRef[:,idt])}')
#err = np.max(np.abs(u - uRef[:,idt]))
#print(f'RKCB2 error for dt = {dt}: {err}')
#
#u_cb2B = rkcb2StepB(t, dt, u0, nN, a, b, nl, ops)
#u = u_cb2B
#if showSol:
#    print(f'Solution 4th order RKCB2B: {u}')
##print(f'Diff: {np.abs(u - uRef[:,idt])}')
#err = np.max(np.abs(u - uRef[:,idt]))
#print(f'RKCB2B error for dt = {dt}: {err}')
#
#u_cb2C = rkcb2StepC(t, dt, u0, nN, a, b, nl, ops)
#u = u_cb2C
#if showSol:
#    print(f'Solution 4th order RKCB2C: {u}')
##print(f'Diff: {np.abs(u - uRef[:,idt])}')
#err = np.max(np.abs(u - uRef[:,idt]))
#print(f'RKCB2C error for dt = {dt}: {err}')
#print(100*'+')
#
#print(f'Diff PC2 vs RKCB2: {np.abs(u_pc2 - u_cb2)}')
#print(f'Diff PC2 vs RKCB2B: {np.abs(u_pc2 - u_cb2B)}')
#print(f'Diff PC2 vs RKCB2C: {np.abs(u_pc2 - u_cb2C)}')
#print(f'Diff RKCB2 vs RKCB2B: {np.abs(u_cb2 - u_cb2B)}')


## Transport 2nd order
#print(100*'/')
#print("Transport equation")
#print(100*'/')
#
#tsp_T = geo.tid(nN, l, 1, {0:0})
#tsp_C = geo.zblk(nN, l, {0:10}).tolil()
#tsp_C[-1,:] = tsp_C[0,:]
#tsp_C[0,:] = 0
#tsp_C = tsp_C.tocoo()
#tsp_L = sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#tsp_LHS = tsp_T - 0.5*dt*tsp_L + tsp_C
#tsp_RHS = tsp_T + 0.5*dt*tsp_L
#
#tsp_q_T = geo.i2(nN, l, {0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#tsp_q_L = geo.i2lapl(nN, l, {0:0})
#tsp_q_C = geo.zblk(nN, l, {0:10})
#tsp_q_NL = geo.i2(nN, l, {0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#tsp_q_LHS = tsp_q_T - 0.5*dt*tsp_q_L + tsp_q_C
#tsp_q_RHS = tsp_q_T + 0.5*dt*tsp_q_L
#
#old_q = np.genfromtxt(f'tsp_old_q_l{l}.dat')
#old_e = np.genfromtxt(f'tsp_old_e_l{l}.dat')
#print(f'Error OLD: {np.max(np.abs(old_q - old_e))}')
#old = old_e
#
#nl_e = np.genfromtxt(f'tsp_nl_e_l{l}.dat')
#nl_q = -np.genfromtxt(f'tsp_nl_q_l{l}.dat')
#nl = nl_e
#print(tsp_q_NL*nl_e - nl_q)
#
#tsp_sol = SLA.spsolve(tsp_LHS, tsp_RHS*old + dt*geo.tid(nN, l, 1, {0:0})*nl)
#tsp_q_sol = SLA.spsolve(tsp_q_LHS, tsp_q_RHS*old + dt*tsp_q_NL*nl)
#q_sol = np.genfromtxt(f'tsp_sol_q_l{l}.dat')
#e_sol = np.genfromtxt(f'tsp_sol_e_l{l}.dat')
#print(f'Error SOL Q vs E: {np.max(np.abs(q_sol - e_sol))}')
#print(tsp_sol)
#print(f'Error sol QuICC: {np.abs(q_sol - tsp_sol)}')
#print(f'Error sol QI: {np.abs(tsp_q_sol - tsp_sol)}')
#print(f'Error sol EPM: {np.abs(e_sol - tsp_sol)}')
#print(f'Error BC tau: {np.max(np.abs(tsp_C*tsp_sol))}')
#print(f'Error BC QuICC: {np.max(np.abs(tsp_C*q_sol))}')
#print(f'Error BC EPM: {np.max(np.abs(tsp_C*e_sol))}')
#
## Transport 2nd order
#print(100*'/')
#print("Toroidal equation")
#print(100*'/')
#
#tor_T = geo.tid(nN, l, 1, {0:0})
#tor_C = geo.zblk(nN, l, {0:12}).tolil()
#tor_C[-1,:] = tor_C[0,:]
#tor_C[0,:] = 0
#tor_C = tor_C.tocoo()
#tor_L = sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#tor_LHS = tor_T - 0.5*dt*tor_L + tor_C
#tor_RHS = tor_T + 0.5*dt*tor_L
#tor_q_NL = geo.i2(nN, l, {0:0})
#
#tor_q_T = geo.i2(nN, l, {0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#tor_q_L = geo.i2lapl(nN, l, {0:0})
#tor_q_C = geo.zblk(nN, l, {0:12})
#tor_q_NL = geo.i2(nN, l, {0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#tor_q_LHS = tor_q_T - 0.5*dt*tor_q_L + tor_q_C
#tor_q_RHS = tor_q_T + 0.5*dt*tor_q_L
#
#old_q = np.genfromtxt(f'tor_old_q_l{l}.dat')
#old_e = np.genfromtxt(f'tor_old_e_l{l}.dat')
#print(f'Error OLD: {np.max(np.abs(old_q - old_e))}')
#old = old_e
#
#nl_e = np.genfromtxt(f'tor_nl_e_l{l}.dat')
#nl_q = -np.genfromtxt(f'tor_nl_q_l{l}.dat')
#nl = nl_e
#print(tor_q_NL*nl_e - nl_q)
#
#tor_sol = SLA.spsolve(tor_LHS, tor_RHS*old + dt*geo.tid(nN, l, 1, {0:0})*nl)
#tor_q_sol = SLA.spsolve(tor_q_LHS, tor_q_RHS*old + dt*tor_q_NL*nl)
#q_sol = np.genfromtxt(f'tor_sol_q_l{l}.dat')
#e_sol = np.genfromtxt(f'tor_sol_e_l{l}.dat')
#print(f'Error SOL Q vs E: {np.max(np.abs(q_sol - e_sol))}')
#print(tor_sol)
#print(f'Error sol QuICC: {np.abs(q_sol - tor_sol)}')
#print(f'Error sol QI: {np.abs(tor_q_sol - tor_sol)}')
#print(f'Error sol EPM: {np.abs(e_sol - tor_sol)}')
#print(f'Error BC tau: {np.max(np.abs(tor_C*tor_sol))}')
#print(f'Error BC QuICC: {np.max(np.abs(tor_C*q_sol))}')
#print(f'Error BC EPM: {np.max(np.abs(tor_C*e_sol))}')
#
#
# Poloidal 4th order
#print(100*'/')
#print("Poloidal equation")
#print(100*'/')
#
#pol_T = geo.tid(nN, l, 2, {0:0})*sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#pol_C = geo.zblk(nN, l, {0:21}).tolil()
#pol_C[-2:,:] = pol_C[0:2,:]
#pol_C[0:2,:] = 0
#pol_C = pol_C.tocoo()
#pol_L = sp.csr_matrix(io.mmread(f'slapl2_l{l}.mtx'))
#pol_LHS = pol_T - 0.5*dt*pol_L + pol_C
#pol_RHS = pol_T + 0.5*dt*pol_L
#
#pol_q_T = geo.i4lapl(nN, l, {0:0}) - geo.i4(nN, l, {0:0})*geo.qid(nN, l, nN-2, {0:0})*sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#pol_q_L = geo.i4lapl2(nN, l, {0:0})
#pol_q_C = geo.zblk(nN, l, {0:21})
#pol_q_NL = geo.i4(nN, l, {0:0}) - geo.i4(nN, l, {0:0})*geo.qid(nN, l, nN-2, {0:0})
#pol_q_LHS = pol_q_T - 0.5*dt*pol_q_L + pol_q_C
#pol_q_RHS = pol_q_T + 0.5*dt*pol_q_L
#
#old_q = np.genfromtxt(f'pol_old_q_l{l}.dat')
#old_e = np.genfromtxt(f'pol_old_e_l{l}.dat')
#print(f'Error OLD: {np.max(np.abs(old_q - old_e))}')
#old = old_e
#
#nl_e = np.genfromtxt(f'pol_nl_e_l{l}.dat')
#nl_q = -np.genfromtxt(f'pol_nl_q_l{l}.dat')
#nl = nl_e
#print(f'Error NL: {pol_q_NL*nl_e - nl_q}')
#
#pol_sol = SLA.spsolve(pol_LHS, pol_RHS*old + dt*geo.tid(nN, l, 2, {0:0})*nl)
#pol_sol_4th = pol_sol
#pol_q_sol = SLA.spsolve(pol_q_LHS, pol_q_RHS*old + dt*pol_q_NL*nl)
#q_sol = np.genfromtxt(f'pol_sol_q_l{l}.dat')
#e_sol = np.genfromtxt(f'pol_sol_e_l{l}.dat')
#print(f'Solution: {pol_sol}')
#print(f'Error SOL Q vs E: {np.max(np.abs(q_sol - e_sol))}')
#print(f'Error sol QuICC: {np.abs(q_sol - pol_sol)}')
#print(f'Error sol QI: {np.abs(pol_q_sol - pol_sol)}')
#print(f'Error sol EPM: {np.abs(e_sol - pol_sol)}')
#print(f'Error BC tau: {np.max(np.abs(pol_C*pol_sol))}')
#print(f'Error BC QuICC: {np.max(np.abs(pol_C*q_sol))}')
#print(f'Error BC EPM: {np.max(np.abs(pol_C*e_sol))}')
#
## Poloidal influence matrix
#print(100*'/')
#print("Poloidal equation")
#print(100*'/')
#
#pol_Q = sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#pol_Q_C = geo.zblk(nN, l, {0:10}).tolil()
#pol_Q_C[-1,:] = pol_Q_C[0,:]
#pol_Q_C[0,:] = 0
#pol_Q_C = pol_Q_C.tocoo()
#pol_LHS_Q = pol_Q + pol_Q_C
#
#pol_T = geo.tid(nN, l, 1, {0:0})
#pol_C = geo.zblk(nN, l, {0:14}).tolil()
#pol_C[-1,:] = pol_C[0,:]
#pol_C[0,:] = 0
#pol_C = pol_C.tocoo()
#pol_L = sp.csr_matrix(io.mmread(f'slapl_l{l}.mtx'))
#pol_LHS = pol_T - 0.5*dt*pol_L + pol_C
#pol_RHS = pol_T + 0.5*dt*pol_L
#
#pol_C_4 = geo.zblk(nN, l, {0:21})
#
#pol_rhs_G = np.zeros(nN);pol_rhs_G[-1] = 1.0
#pol_sol_g_G = SLA.spsolve(pol_LHS_Q, pol_rhs_G)
#pol_sol_G = SLA.spsolve(pol_LHS, -dt*geo.tid(nN, l, 1, {0:0})*pol_sol_g_G)
#pol_sol_G /= (pol_Q_C*pol_sol_G)[-1]
#print(f'Green\'s function: {pol_sol_G}')
#
#pol_q_Q =  geo.i2lapl(nN, l, {0:0})
#pol_q_Q_C = geo.zblk(nN, l, {0:10})
#pol_q_LHS_Q = pol_q_Q + pol_q_Q_C
#
#pol_q_T = geo.i2(nN, l,{0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#pol_q_C = geo.zblk(nN, l, {0:14})
#pol_q_L = geo.i2lapl(nN, l, {0:0})
#pol_q_LHS = pol_q_T - 0.5*dt*pol_q_L + pol_q_C
#pol_q_RHS = pol_q_T + 0.5*dt*pol_q_L
#pol_q_NL = geo.i2(nN, l, {0:0}) - geo.i2(nN, l, {0:0})*geo.qid(nN, l, nN-1, {0:0})
#
#pol_q_rhs_G = np.zeros(nN);pol_q_rhs_G[0] = 1.0
#pol_q_sol_g_G = SLA.spsolve(pol_q_LHS_Q, pol_q_rhs_G)
#pol_q_sol_G = SLA.spsolve(pol_q_LHS, -dt*pol_q_NL*pol_q_sol_g_G)
#pol_q_sol_G /= (pol_q_Q_C*pol_q_sol_G)[0]
#print(f'Green\'s QI function: {pol_q_sol_G}')
#
#old_q = np.genfromtxt(f'pol_old_q_l{l}.dat')
#old_e = np.genfromtxt(f'pol_old_e_l{l}.dat')
#print(f'Error OLD: {np.max(np.abs(old_q - old_e))}')
#old = old_e
#
#nl_e = np.genfromtxt(f'pol_nl_e_l{l}.dat')
#nl_q = -np.genfromtxt(f'pol_nl_q_l{l}.dat')
#nl = nl_e
#print(f'Error NL: {pol_q_NL*nl_e - nl_q}')
#
#pol_sol_g = SLA.spsolve(pol_LHS_Q, geo.tid(nN, l, 1, {0:0})*nl)
#pol_sol = SLA.spsolve(pol_LHS, pol_RHS*old - dt*geo.tid(nN, l, 1, {0:0})*pol_sol_g)
#bc_err = (pol_Q_C*pol_sol)[-1]
#pol_sol -= bc_err*pol_sol_g_G
#print(f'Solution: {pol_sol}')
#
#pol_q_sol_g = SLA.spsolve(pol_q_LHS_Q, pol_q_NL*nl)
#pol_q_sol = SLA.spsolve(pol_q_LHS, pol_q_RHS*old - dt*pol_q_NL*pol_q_sol_g)
#bc_err = (pol_q_Q_C*pol_q_sol)[0]
#pol_q_sol -= bc_err*pol_q_sol_g_G
#print(f'Solution QI: {pol_q_sol}')
#
#q_sol = np.genfromtxt(f'pol_sol_q_l{l}.dat')
#e_sol = np.genfromtxt(f'pol_sol_e_l{l}.dat')
#print(f'Error SOL Q vs E: {np.max(np.abs(q_sol - e_sol))}')
#print(f'Error sol vs 4th order: {np.abs(pol_sol_4th - pol_sol)}')
##print(f'Error sol QuICC: {np.abs(q_sol - pol_sol)}')
#print(f'Error sol QI: {np.abs(pol_q_sol - pol_sol)}')
#print(f'Error sol EPM: {np.abs(e_sol - pol_sol)}')
#print(f'Error BC tau: {np.max(np.abs(pol_C_4*pol_sol))}')
##print(f'Error BC QuICC: {np.max(np.abs(pol_C_4*q_sol))}')
#print(f'Error BC EPM: {np.max(np.abs(pol_C_4*e_sol))}')
