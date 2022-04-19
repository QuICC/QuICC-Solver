"""Check accuracy for the theta direction in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl

import quicc.transform.shell as transf
import quicc.geometry.spherical.shell_sh as shell


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        pl.semilogy(np.abs(err))
        pl.title(title)
        pl.show()

def sh_to_phys(expr, tg):
    """Convert sympy expression to grid values"""

    th = sy.Symbol('th')
    ph = sy.Symbol('ph')
    func = sy.utilities.lambdify((th, ph), expr)
    t = func(0.1,0)
    phys = np.zeros((len(tg),1))
    for i in range(0,len(tg)):
        phys[i] = func(tg[i],0).evalf()

    return phys

def test_forward(op, maxl, m, slhs, ssol, tg):
    """Perform a forward operation test"""

    print("\tForward test")
    rhs = op*slhs
    mesh = sh_to_phys(ssol, tg)
    sol = transf.totleg(mesh, maxl, m)
    err = np.abs(sol - rhs)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def coriolisdr(maxl, m, tg):
    """Accuracy test for coriolisdr operator"""

    print("coriolisdr:")
    th = sy.Symbol('th')
    ph = sy.Symbol('ph')
    A = shell.coriolisdr(maxl, m)
    sphys = [np.random.ranf() for i in np.arange(m,maxl+1)]
    ssol = 0
    for i, c in enumerate(sphys):
        ssol = ssol + c*sy.cos(th)*sy.Ynm(m+i, m, th, ph)

    test_forward(A, maxl, m, sphys, ssol, tg)

def coriolis_r(maxl, m, tg):
    """Accuracy test for coriolis_r operator"""

    print("coriolis_r:")
    th = sy.Symbol('th')
    ph = sy.Symbol('ph')
    A = shell.coriolis_r(maxl, m)
    sphys = [np.random.ranf() for i in np.arange(m,maxl+1)]
    ssol = 0
    for i, c in enumerate(sphys):
        ssol = ssol + c*sy.cos(th)*sy.Ynm(m+i, m, th, ph)
    test_forward(A, sphys, ssol, tg)


if __name__ == "__main__":
    # Set test parameters
    maxl = 4
    m = np.random.randint(1, maxl//2)
    tg = transf.tgrid(maxl, m)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    #zblk(nx, xg)
    coriolisdr(maxl, m, tg)
    #coriolis_r(maxl, m, tg)
