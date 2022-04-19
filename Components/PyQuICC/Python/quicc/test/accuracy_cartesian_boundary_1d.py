"""Check accuracy for cartesian 1D operators"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import mpmath
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
#import matplotlib.pylab as pl
import scipy.io as io

import quicc.transform.cartesian as transf
transf.min_grid = 1
import quicc.geometry.cartesian.cartesian_1d as c1d
import quicc.geometry.cartesian.cartesian_boundary_1d as bc1d


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.lambdify(x, expr)
    return func(grid)

def test_value(op, res_expr, sol, grid):
    """Perform a value operation test"""

    print("\tValue test")
    x = sy.Symbol('x')
    lhs = transf.tocheb(x_to_phys(res_expr,grid))
    rhs = op*lhs
    print("Matrix integral: " + str(rhs))
    print("Sympy integral: " + str(sol))
    err = np.abs(rhs - sol)
    print("\t\tValue error: " + str(err))

def test_tau_insulating(nx, grid, cscale = 2.0):
    """Test insulating boundary condition"""
    print("Insulating boundary condition (c = "+str(cscale)+")")
    k_perp = 7.0
    bc = {'cscale':cscale, 'k_perp':k_perp}
    x = sy.Symbol('x')
    xexpr = x
    xexpr = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    chebx = transf.tocheb(x_to_phys(xexpr, xg))[0:nx]
    tau_line_top = bc1d.tau_insulating(nx, 1, bc)
    tau_line_bot = bc1d.tau_insulating(nx, -1, bc)
    tau_line = bc1d.tau_insulating(nx, 0, bc)
    bcVal_top = np.dot(tau_line_top, chebx)[0]
    bcVal_bot = np.dot(tau_line_bot, chebx)[0]
    bcVal = np.dot(tau_line, chebx)
    ref_top = x_to_phys(cscale*sy.diff(xexpr,x) + k_perp*xexpr,np.array([1.0]))[0]
    ref_bot = x_to_phys(cscale*sy.diff(xexpr,x) - k_perp*xexpr,np.array([-1.0]))[0]
    ref = np.array([(ref_top+ref_bot)/2.0, (ref_top-ref_bot)/2.0])
    # Absolute error
    print("Absolute error:")
    print(bcVal_top - ref_top)
    print(bcVal_bot - ref_bot)
    print(bcVal - ref)
    # Relative error
    print("Relative error:")
    print((bcVal_top - ref_top)/ref_top)
    print((bcVal_bot - ref_bot)/ref_bot)
    print((bcVal - ref)/ref)
    print("")

if __name__ == "__main__":
    # Set test parameters
    nx = 32
    xg = transf.grid(nx)

    test_tau_insulating(nx, xg, cscale = 1.0)
    test_tau_insulating(nx, xg, cscale = 2.0)
