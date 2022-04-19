"""Module to compute the stencils for Galerkin basis for the spectral Jacobi operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
from sympy.abc import a,b,c,d,e,l,m,n,q
import numpy as np

cs = [a, b, c, d, e]

w_alpha_chebyshev = -sy.Rational(1,2)
w_alpha_legendre = 0
w_beta_minhalf = l - sy.Rational(1,2)
w_beta_zero = l
w_beta_plushalf = l + sy.Rational(1,2)

w_a = w_alpha_chebyshev
w_b = w_beta_minhalf

def w(n, l, a, b):
    expr = sy.gamma(n + a + 1)/(sy.gamma(n + 1)*sy.gamma(a + 1))
    return expr

def dw(n, l, a, b):
    expr = 2*(n + l)*w(n-1,l,a+1,b+1) + l*w(n,l,a,b)
    return expr

def rdr_1w(n, l, a, b):
    expr = 2*(n + l)*w(n-1,l,a+1,b+1) + (l - 1)*w(n,l,a,b)
    return expr

def r_1drw(n, l, a, b):
    expr = 2*(n + l)*w(n-1,l,a+1,b+1) + (l + 1)*w(n,l,a,b)
    return expr

def insulating(n, l, a, b):
    expr = 2*(n + l)*w(n-1,l,a+1,b+1) + (2*l + 1)*w(n,l,a,b)
    return expr

def ddw(n, l, a, b):
    expr = 4*(n + l + 1)*(n + l)*w(n-2,l,a+2,b+2) + 2*(2*l + 1)*(n + l)*w(n-1,l,a+1,b+1) + l*(l - 1)*w(n,l,a,b)
    return expr

def laplh(n, l, a, b):
    expr = 4*(n + l + 1)*(n + l)*w(n-2,l,a+2,b+2) + 4*(l + 1)*(n + l)*w(n-1,l,a+1,b+1)
    return expr

def lapl2h(n, l, a, b):
    expr = 16*(n + l + 3)*(n + l + 2)*(n + l + 1)*(n + l)*w(n-4,l,a+4,b+4) + 32*(l + 2)*(n + l + 2)*(n + l + 1)*(n + l)*w(n-3,l,a+3,b+3) + 16*(l + 2)*(l + 1)*(n + l + 1)*(n + l)*w(n-2,l,a+2,b+2)
    return expr

def stencil(bN, cfunc, l, alpha, beta, stride = 1, show_col = False, show_cond = False, show_sol = False):
    C = sy.IndexedBase('C')
    rN = range(0,bN)
    coeffs = [cs[i+1] for i in rN]
    try:
        cond = list([])
        for cf in cfunc:
            cond.append(np.sum([coeffs[i]*cf(q+i, l, alpha, beta) for i in rN]))
    except:
        cond = np.sum([coeffs[i]*cfunc(q+i, l, alpha, beta) for i in rN])
    if show_cond:
        print('+'*5 + ' Condition(s) ' + '+'*5)
        print(cond)
    sol = sy.solvers.solve(cond, *coeffs)
    if type(sol) is dict:
        sol = [sol[k] for k in coeffs[:-1]]
        sol.append(coeffs[-1])
    else:
        sol = sol[0]
    sol = [sy.simplify(s) for s in sol]
    if show_sol:
        print('+'*5 + ' Solution(s) ' + '+'*5)
        print(sol)
    # Column solutions
    colSol = [sy.simplify(col.subs(coeffs[-1], C[q]).subs(q,m)) for col in sol]
    colSol = [sy.simplify(colSol[i]/(colSol[0]/C[m])) for i in rN]
    colZeroSol = [sy.simplify(col.subs(m,0)) for col in colSol]
    if show_col:
        print('+'*5 + ' Column solutions ' + '+'*5)
        print(colSol)
        print('+'*5 + ' Column zero solution ' + '+'*5)
        print(colZeroSol)
    # Row solutions
    rowSol = [sy.simplify(sol[-(i+1)].subs(coeffs[-1], C[q+1]).subs(q,n + stride*(i-bN))) for i in rN]
    rowSol = [sy.simplify(rowSol[i]/(rowSol[-1]/C[n])) for i in rN]
    rowZeroSol = [row/(colZeroSol[0]/C[0]) for row in colZeroSol[::-1]]
    subs = [(C[n-i],1) for i in rN]
    rowUnitSol = [sy.simplify(row.subs(subs)) for row in rowSol]
    print('+'*5 + ' Row solution ' + '+'*5 )
    print(rowSol)
    print('+'*5 + ' Row zero solution ' + '+'*5)
    print(rowZeroSol)
    print('+'*5 + ' Row unit solution ' + '+'*5)
    print(rowUnitSol)

print('='*30)
print('u(1) = 0')
print('-'*30)
stencil(2, w, l, w_a, w_b)
print('='*30)

print('='*30)
print('Du(1) = 0')
print('-'*30)
stencil(2, dw, l, w_a, w_b)
print('='*30)

print('='*30)
print('1/r D r u(1) = 0')
print('-'*30)
stencil(2, r_1drw, l, w_a, w_b)
print('='*30)

print('='*30)
print('r D 1/r u(1) = 0')
print('-'*30)
stencil(2, rdr_1w, l, w_a, w_b)
print('='*30)

print('='*30)
print('Spherical insulating boundary')
print('-'*30)
stencil(2, insulating, l, w_a, w_b)
print('='*30)

print('='*30)
print('No-slip, no penetration')
print('-'*30)
stencil(3, (dw, w), l, w_a, w_b)
print('='*30)

print('='*30)
print('Stress-free, no penetration')
print('-'*30)
stencil(3, (ddw, w), l, w_a, w_b)
print('='*30)
