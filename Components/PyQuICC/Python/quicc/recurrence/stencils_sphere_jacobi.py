"""Module to compute the stencils for Galerkin basis for the spectral Jacobi operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
from sympy.abc import a,b,c,d,e,l,m,n,q
from sympy import linsolve
import numpy as np

cs = [a, b, c, d, e]

w_alpha_chebyshev = -sy.Rational(1,2)
w_alpha_legendre = 0
w_beta_minhalf = l - sy.Rational(1,2)
w_beta_zero = l
w_beta_plushalf = l + sy.Rational(1,2)

w_type  = "SphEnergy"

# Setup for Chebyshev type
if w_type == "Chebyshev":
    w_a = w_alpha_chebyshev
    w_b = w_beta_minhalf
# Setup for SphEnergy type
elif w_type == "SphEnergy":
    w_a = w_alpha_legendre
    w_b = w_beta_plushalf

def w(n, l, a, b):
    expr = sy.gamma(n + a + 1)/(sy.gamma(n + 1)*sy.gamma(a + 1))
    return expr

def dw(n, l, a, b):
    expr = 2*(n + a + b + 1)*w(n-1,l,a+1,b+1) + l*w(n,l,a,b)
    return expr

def rdr_1w(n, l, a, b):
    expr = 2*(n + a + b + 1)*w(n-1,l,a+1,b+1) + (l - 1)*w(n,l,a,b)
    return expr

def r_1drw(n, l, a, b):
    expr = 2*(n + a + b + 1)*w(n-1,l,a+1,b+1) + (l + 1)*w(n,l,a,b)
    return expr

def insulating(n, l, a, b):
    expr = 2*(n + a + b + 1)*w(n-1,l,a+1,b+1) + (2*l + 1)*w(n,l,a,b)
    return expr

def ddw(n, l, a, b):
    expr = 4*(n + a + b + 1)*(n + a + b + 2)*w(n-2,l,a+2,b+2) + 2*(2*l + 1)*(n + a + b + 1)*w(n-1,l,a+1,b+1) + l*(l - 1)*w(n,l,a,b)
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
        cond = [np.sum([coeffs[i]*cfunc(q+i, l, alpha, beta) for i in rN])]
    if show_cond:
        print('+'*5 + ' Condition(s) ' + '+'*5)
        print(cond)
    print(coeffs)
    print(tuple(coeffs[1:]))
    sol = linsolve(cond, tuple(coeffs[1:]))
    sol = [sy.expand_func(s).expand().factor() for s in sol.args[0]]
    sol.insert(0,coeffs[0])
    if show_sol:
        print('+'*5 + ' Solution(s) ' + '+'*5)
        print(sol)
    # Column solutions
    colSol = [sy.simplify(col.subs(coeffs[0], C[q]).subs(q,m)) for col in sol]
    colSol = [sy.simplify(colSol[i]/(colSol[0]/C[m])) for i in rN]
    colSol = [sy.factor(sy.expand(s)) for s in colSol]
    colZeroSol = [sy.simplify(col.subs(m,0)) for col in colSol]
    colZeroSol = [sy.factor(sy.expand(s)) for s in colZeroSol]
    if show_col:
        print('+'*5 + ' Column solutions ' + '+'*5)
        print(colSol)
        print('+'*5 + ' Column zero solution ' + '+'*5)
        print(colZeroSol)
    # Row solutions
    rowSol = [sy.simplify(sol[i].subs(coeffs[0], C[q]).subs(q,n - stride*i)) for i in rN]
    #rowSol = [sy.simplify(rowSol[i]/(rowSol[0]/C[n-stride*i])) for i in rN]
    rowSol = [sy.factor(sy.expand(s)) for s in rowSol]
    rowZeroSol = [row/(colZeroSol[0]/C[0]) for row in colZeroSol[::-1]]
    rowZeroSol = [sy.factor(sy.expand(s)) for s in rowZeroSol]
    subs = [(C[n-i],1) for i in rN]
    rowUnitSol = [sy.simplify(row.subs(subs)).cancel() for row in rowSol]
    rowUnitSol = [sy.factor(sy.expand(s)) for s in rowUnitSol]
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
stencil(3, (dw, w), l, w_a, w_b, show_sol = True)
print('='*30)

print('='*30)
print('Stress-free, no penetration')
print('-'*30)
stencil(3, (ddw, w), l, w_a, w_b, show_sol = True)
print('='*30)
