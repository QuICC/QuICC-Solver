"""Check accuracy for radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import scipy.special as special

import quicc.transform.cylinder_worland as transf
import quicc.geometry.cylindrical.cylinder_radius_worland as geo
import quicc.geometry.worland.wnl as wnl
np.set_printoptions(precision=15)


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    print(expr)
    return func(grid)

def test_forward(op, m, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.torspec(x_to_phys(res_expr,grid), m, op.shape[0])
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torspec(t, m, op.shape[0])
    #print(rhs.T)
    #print("=================================================================")
    #print(sol.T)
    err = np.abs(rhs[0:-(1+q)] - sol[q:-1])
    relerr = err/(1.0 + np.abs(sol[q:-1]))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err.T)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr.T)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, m, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    rhs = transf.torspec(x_to_phys(res_expr,grid), m, opA.shape[0])
    rhs = rhs[0:opA.shape[0]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = np.reshape(lhs, (lhs.shape[0],1))
    sol = transf.torspec(x_to_phys(sol_expr,grid), m, opA.shape[0])
    sol = sol[0:opA.shape[0]]
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err.T)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr.T)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

def zblk(nr, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.zblk(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = 0
        test_forward(A, m, sphys, ssol, rg, 0)

def i2(nr, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,2*nr,2)])
        ssol = 4**2*sy.integrate(sy.integrate(sphys*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 1)

def i2laplh(nr, rg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')

    print("\tForward:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(j) for j in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2)*x**(-m)
        ssol = 4**2*sy.integrate(sy.integrate(stmp*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 1)

    print("\tbc = 10:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, {0:10}).tocsr()
        B = geo.i2(nr, m, {0:0}).tocsr()
        ssol = sy.expand((x**2-1)*np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)]))
        sphys = sy.expand(sy.diff(ssol,x,x) + sy.diff(ssol,x)/x - m**2*ssol/x**2)
        test_backward_tau(A, B, m, sphys, ssol, rg)

    print("\tbc = 11:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, {0:11}).tocsr()
        B = geo.i2(nr, m, {0:0}).tocsr()
        ssol = sy.expand((x**2-1)**2*np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)]))
        sphys = sy.expand(sy.diff(ssol,x,x) + sy.diff(ssol,x)/x - m**2*ssol/x**2)
        test_backward_tau(A, B, m, sphys, ssol, rg)

def i4(nr, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        ssol = 4**4*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sphys*x,x)*x,x)*x)*x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 2)

def i4laplh(nr, rg):
    """Accuracy test for i4laplh operator"""

    print("i4laplh:")
    x = sy.Symbol('x')
    print("\tForward:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2)*x**(-m)
        ssol = 4**4*sy.integrate(sy.integrate(sy.integrate(sy.integrate(stmp*x,x)*x,x)*x)*x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 2)

def i4lapl2h(nr, rg):
    """Accuracy test for i4lapl2h operator"""

    print("i4lapl2h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4lapl2h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2)
        stmp = sy.expand(sy.diff(stmp,x,x) + sy.diff(stmp,x)/x - m**2*stmp/x**2)*x**(-m)
        ssol = 4**4*sy.integrate(sy.integrate(sy.integrate(sy.integrate(stmp*x,x)*x,x)*x)*x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 2)

def i6(nr, rg):
    """Accuracy test for i6 operator"""

    print("i6:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        ssol = 4**6*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sphys*x,x)*x,x)*x,x)*x,x)*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 3)

def i6laplh(nr, rg):
    """Accuracy test for i6laplh operator"""

    print("i6laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2)*x**(-m)
        ssol = 4**6*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(stmp*x,x)*x,x)*x,x)*x,x)*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 3)

def i6lapl2h(nr, rg):
    """Accuracy test for i6lapl2h operator"""

    print("i6lapl2h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6lapl2h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2)
        stmp = sy.expand(sy.diff(stmp,x,x) + sy.diff(stmp,x)/x - m**2*stmp/x**2)*x**(-m)
        ssol = 4**6*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(stmp*x,x)*x,x)*x,x)*x,x)*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 3)

def i6lapl3h(nr, rg):
    """Accuracy test for i6lapl3h operator"""

    print("i6lapl3h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6lapl3h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([x**i for i in np.arange(0,2*nr,2)])
        stmp = sy.expand(sy.simplify(sy.diff(sphys*x**m,x,x) + sy.diff(sphys*x**m,x)/x - m**2*sphys*x**m/x**2))
        stmp = sy.expand(sy.simplify(sy.diff(stmp,x,x) + sy.diff(stmp,x)/x - m**2*stmp/x**2))
        stmp = sy.expand(sy.simplify(sy.diff(stmp,x,x) + sy.diff(stmp,x)/x - m**2*stmp/x**2))*x**(-m)
        ssol = 4**6*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(stmp*x,x)*x,x)*x,x)*x,x)*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 3)

def test_worland(nr, m):
    """Test Worland transform"""

    # Create physical space test function
    t = wnl.eval_poly(0, m, nr)
    for i in range(1, nr-m):
        t += np.random.ranf()*wnl.eval_poly(i, m, nr)
    t = np.matrix(t).T

    # Compute spectral expansion
    s = transf.torspec(t, m, nr-m)

    # Project spectral expansion to physical space
    st = transf.torphys(s, m, nr)

    # Print error for transform loop
    print(np.abs(st.T-t.T))

def test_fft(n, m):
    """Test Worland transform"""

    nr = np.ceil(3.0*n/2. + 3.0*m/4.0 + 1.0)
    # Create physical space test function
    t = wnl.eval_poly(0, m, nr)
    for i in range(1, n):
        t += np.random.ranf()*wnl.eval_poly(i, m, nr)
    t = np.matrix(t).T

    # Compute spectral expansion
    s = transf.torspec(t, m, n)

    # Project spectral expansion to physical space
    st = transf.torphys(s, m, nr)

    # Print error for transform loop
    print(np.max(np.abs(st.T-t.T)))

    import scipy.fftpack as fft
    import scipy.linalg as linalg
    cheb2worl  = []
    for i in range(0, n):
        cheb2worl.append(fft.dct(wnl.eval_poly(i, m, nr))/(2*nr))
    cheb2worl = np.matrix(cheb2worl)
    f = fft.dct(t.T).T/(2*nr)
    mat = cheb2worl[:,m//2:m//2+n].T
    print(cheb2worl[:,16].T)
    print(cheb2worl[:,17].T)
    fs = linalg.solve_triangular(mat, f[m//2:m//2+n])

    # Print error for FFT transform loop
    print(np.max(np.abs(s.T-fs.T)))

def test_endpoints(n, m):
    """Test endpoint values"""

    t = []
    t.append(wnl.eval_bc_poly(n, m))
    t.append(wnl.eval_bc_diff(n, m))
    t.append(wnl.eval_bc_diff2(n, m))
    t.append(wnl.eval_bc_diff3(n, m))
    t.append(wnl.eval_bc_diff4(n, m))
    t.append(wnl.eval_bc_rdiffdivr(n, m))
    t.append(wnl.eval_bc_divrdiffr(n, m))
    t.append(wnl.eval_bc_laplh_cyl(n, m))
    t.append(wnl.eval_bc_dlaplh_cyl(n, m))
    t.append(wnl.eval_bc_lapl2h_cyl(n, m))
    t = np.array(t).T
    np.savetxt("cylinder_worland_endpoints_l"+str(m)+".dat", t)
    print("Endpoint values:")
    print("Value, Diff, Diff2, Diff3, Diff4, rdiffdivr, divrdiff, laplh, dlaplh, lapl2h")
    print(t)

def test_stencils(n, m):
    """Test endpoint values"""

    print("Stencil errors: m = " + str(m))
    print("\t Value:")
    err = geo.radbc.stencil_value(n,m).T*wnl.eval_bc_poly(n, m)
    print(np.max(np.abs(err)))
    print("\t Diff:")
    val = geo.radbc.stencil_diff(n,m).T*wnl.eval_bc_diff(n, m)
    print(np.max(np.abs(err)))
    print("\t Value + Diff:")
    t = np.array([wnl.eval_bc_poly(n, m), wnl.eval_bc_diff(n, m)])
    err = geo.radbc.stencil_value_diff(n,m).T*t.T
    print(np.max(np.abs(err),0))
    # print("\t Value + Diff2:")
    #t = np.array([wnl.eval_bc_poly(n, m), wnl.eval_bc_diff2(n, m)])
    #err = geo.radbc.stencil_value_diff2(n,m).T*t.T
    #print(np.max(np.abs(err),0))
    print("\t Value + Laplh:")
    t = np.array([wnl.eval_bc_poly(n, m), wnl.eval_bc_laplh_cyl(n, m)])
    err = geo.radbc.stencil_value_laplh(n,m).T*t.T
    print(np.max(np.abs(err),0))

def divr(n):
    """Test 1/r operator"""

    for i in range(0,1):
        #m = np.random.randint(1, n-1)
        #m = m + (m+i)%2
        m = 32
        print("1/r operator: m = " + str(m))
#        A = geo.divr(n, m, {0:0}).tolil()
#        #A[0,:] = 0
#        #A[0,-1] = 1
#        A = A.tocsr()

        # Build unit W_n^{l-1} vector
        #v = np.array([1.0/i**4 for i in range(1,n+1)])
        v = np.array([1.0 for i in range(1,n+1)])
        v[-1:] = 0
        v = v.reshape((n,1))

        # Remove divergent part (multiply by r^2)
        pv = transf.torphys(v, m-1, 2*n+m)
        g = wnl.get_grid(pv.shape[0]/2)
        g = g.reshape((g.shape[0],1))
        v = transf.torspec(np.multiply(g**2,pv), m-1, n)
        v = v.reshape((n,1))
        pv = transf.torphys(v, m-1, 2*n+m)

        # Make bad v
        v_bad = v.copy()
        v_bad[0] += 1.0
        pv_bad = transf.torphys(v_bad, m-1, 2*n+m)

        print(np.asscalar(v[0]))
        #
        for j in range(4, 24):
            vt = transf.torphys(v_bad, m-1, 2*n+m)
            t = 0.5*transf.torspec(vt, m-1, j)
            mid = np.array([(-1.0)**i*np.exp(special.gammaln(i+m-0.5)-special.gammaln(i+1.0)-special.gammaln(m-0.5))*wnl.get_invnorm(i,m-1)/wnl.get_invnorm(0,m-1) for i in range(0,j)])
            mid = mid.reshape((mid.shape[0],1))
            mid[0] = 0
            corr = -np.asscalar(np.dot(mid[:,0],t[:,0]))
            print(corr)
            #v_bad[0] = corr
            print(v_bad[0] - v[0])

        print("---------------------------------------------------")
        mid = np.array([(-1.0)**i*np.exp(special.gammaln(i+m-0.5)-special.gammaln(i+1.0)-special.gammaln(m-0.5))*wnl.get_invnorm(i,m-1)/wnl.get_invnorm(0,m-1) for i in range(0,n)])
        mid = mid.reshape((mid.shape[0],1))
        mid[0] = 0
#        print(np.multiply(mid, v))
#        print(np.dot(mid[:,0],v_bad[:,0]))
#        print(np.dot(mid[:,0],v_bad[:,0]) + v[0])
#        #print(np.dot(mid[:,0],v_bad[:,0]) + 2.008339172142914)
#        print(np.dot(mid[:,0],v_bad[:,0]) + 2.001601876723355)
#        #print(v_bad)
#        #print(mid[1]*v_bad[1,0])
#        #print((np.multiply(mid,v_bad)).T)
#        print(np.cumsum(np.multiply(mid,v_bad)))
#        print(v)
        #print(np.cumsum(np.multiply(mid,v_bad[:,0])))

#        # Build clean RHS (shift up and zero last)
#        rhs = np.zeros(v.shape)
#        rhs = v_bad.copy()
#        #rhs[0:-2] = v[1:-1]
#        #rhs[0:-2] = v[0:-2]
#
#        import scipy.io as io
#        io.mmwrite('matrix_divr.mtx', A)
#
#        # Solve for clean LHS
#        #print(rhs)
#        lhs = spsplin.spsolve(A, rhs)
#        import numpy.linalg as LA
#        #print(lhs)
#        lhs = lhs.reshape((n,1))
#        #lhs[-2] = 0
#        #print(rhs - A*lhs)
#
#        # Plot solution
#        import matplotlib.pylab as pl
#        from matplotlib import rc
#        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#        ## for Palatino and other serif fonts use:
#        #rc('font',**{'family':'serif','serif':['Palatino']})
#        rc('text', usetex=True)
#
#        pl.subplot(221)
#        pl.title("Spectrum")
#        pl.semilogy(np.abs(v), 'r<', label = '$W_^{m-1}$', linewidth=2, markersize=10)
#        pl.semilogy(np.abs(v_bad), 'c+', label = '$W_^{m-1} + \\alpha W_0^{m-1}$', linewidth=2, markersize=10)
#        pl.semilogy(np.abs(lhs), 'ko', label = '$W_n^m = \\frac{W_n^{m-1}}{r}$', linewidth=2, markersize=10)
#        r = 0.5*transf.torspec(pv/g, m, n)
#        r = r.reshape((n,1))
#        print(r)
#        pl.semilogy(np.abs(r), 'mx', label = 'Manual clean', linewidth=2, markersize=10)
#        r_bad = 0.5*transf.torspec(pv_bad/g, m, n)
#        r_bad = r_bad.reshape((n,1))
#        pl.semilogy(np.abs(r_bad), 'cd', label = 'Manual bad', linewidth=2, markersize=10)
#        pl.legend()
#
#        pl.subplot(222)
#        pl.title("Physical")
#        pl.plot(g, pv, 'r-', label = '$W_^{m-1}$', linewidth = 2)
#        pl.plot(g, pv/g, 'bx',  label = '$W_n^{m-1}/r$',linewidth = 1)
#        pl.plot(g, pv_bad, 'c-', label = '$W_n^{m-1} + \\alpha W_0^{m-1}$', linewidth = 2)
#        pl.plot(g, pv_bad/g, 'm-',  label = '$(W_n^{m-1} + \\alpha W_0^{m-1})/r$',linewidth = 1)
#        plhs = transf.torphys(lhs, m, 2*n+m)
#        pl.plot(g, plhs, 'k-',  label = '$W_n^m = W_n^{m-1}/r$',linewidth = 2)
#        pr = transf.torphys(r, m, 2*n+m)
#        pl.plot(g, pr, 'k:',  label = '$W_n^m = W_n^{m-1}/r$',linewidth = 2)
#        pr_bad = transf.torphys(r_bad, m, 2*n+m)
#        pl.plot(g, pr_bad, 'm:',  label = '$W_n^m = W_n^{m-1}/r$',linewidth = 2)
#        pl.legend(loc = 'upper left')
#
#        pl.subplot(223)
#        pl.title('Error')
#        pl.semilogy(g, np.abs(pv/g-plhs), 'k+-', label = 'Clean - Op', linewidth = 2)
#        pl.semilogy(g, np.abs(pv/g-pr), 'kx:', label = 'Clean - Man clean', linewidth = 2)
#        pl.semilogy(g, np.abs(pv/g-pr_bad), 'mo:', label = 'Clean - Man bad', linewidth = 2)
#        pl.legend(loc = 'lower right')
#
#        pl.subplot(224)
#        pl.semilogy(np.abs(lhs - r), 'k+-', label = 'Clean - Op', linewidth = 2)
#        pl.semilogy(np.abs(lhs - r_bad), 'kx:', label = 'Clean - Man clean', linewidth = 2)
#        pl.legend()
#        pl.show()


if __name__ == "__main__":
    # Set test parameters
    n = 64
    nr = int(np.ceil(3.0*n/2.0 + 3.0*n/4.0 + 1))
    print("Grid: " + str((n, nr)))
    rg = wnl.get_grid(nr)

    # run tests
#    test_worland(nr, 110)
#    for i in range(0,4):
#        test_endpoints(n, i)
#    test_stencils(n, 0)
#    test_fft(nr, 32)
    divr(n)
#    zblk(nr, rg)
#    i2(n, rg)
#    i2laplh(nr, rg)
#    i4(n, rg)
#    i4laplh(nr, rg)
#    i4lapl2h(nr, rg)
#    i6(nr, rg)
#    i6laplh(nr, rg)
#    i6lapl2h(nr, rg)
#    i6lapl3h(nr, rg)
