"""Module to compute the recurrence relations for the chebyshev spectral operators of the linear map y = ax + b"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def yp(p):
    """y^p operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':p, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iq(q):
    """I^q operator"""

    # Setup terms in recurrence
    terms = [{'q':q, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqyp(q,p):
    """I^q y^p operator"""

    # Setup terms in recurrence
    terms = [{'q':q, 'p':p, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqds(q, s):
    """I^q D^s operator"""

    # Setup terms in recurrence
    terms = [{'q':q, 'p':0, 'd':s, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqypds(q,p,s):
    """I^q y^p D^s operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':q, 'p':p, 'd':s, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqypd1y1(q,p):
    """I^q y^p of 1/y D y operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':q, 'p':p+1, 'd':1, 'c':1}, {'q':q, 'p':p, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqlapl_3d(q):
    """I^q of cartesian laplacian 3D (1 chebyshev)"""

    # Setup terms in recurrence
    k, l = sympy.symbols('k l')
    terms = [{'q':q, 'p':0, 'd':2, 'c':1},{'q':q, 'p':0, 'd':0, 'c':-k**2},{'q':q, 'p':0, 'd':0, 'c':-l**2}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqlapl_2d(q):
    """I^q of cartesian laplacian 2D (1 chebyshev)"""

    # Setup terms in recurrence
    k = sympy.symbols('k')
    terms = [{'q':q, 'p':0, 'd':2, 'c':1},{'q':q, 'p':0, 'd':0, 'c':-k**2}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqypsphlapl(q,p):
    """I^q y^p of spherical laplacian"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':q, 'p':p, 'd':2, 'c':1}, {'q':q, 'p':p-1, 'd':1, 'c':2}, {'q':q, 'p':p-2, 'd':0, 'c':-l*(l+1)}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqrpsphlapl2(q,p):
    """I^q y^p of spherical bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':q, 'p':p, 'd':4, 'c':1}, {'q':q, 'p':p-1, 'd':3, 'c':4}, {'q':q, 'p':p-2, 'd':2, 'c':-2*l*(l+1)}, {'q':q, 'p':p-4, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4y4sphlaplyd1y1():
    """Spherical shell 4th integral of r^4 laplacian D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':3, 'c':1}, {'q':4, 'p':3, 'd':2, 'c':3}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqlapl2_3d(q):
    """I^q cartesian laplacian in 3D (1 Chebyshev)"""

    # Setup terms in recurrence
    k, l = sympy.symbols('k l')
    terms = [{'q':q, 'p':0, 'd':4, 'c':1},{'q':q, 'p':0, 'd':0, 'c':k**4},{'q':q, 'p':0, 'd':0, 'c':l**4},{'q':q, 'p':0, 'd':2, 'c':-2*k**2},{'q':q, 'p':0, 'd':2, 'c':-2*l**2},{'q':q, 'p':0, 'd':0, 'c':2*k**2*l**2}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def iqlapl2_2d(q):
    """I^q cartesian laplacian in 2D (1 Chebyshev)"""

    # Setup terms in recurrence
    k = sympy.symbols('k')
    terms = [{'q':q, 'p':0, 'd':4, 'c':1},{'q':q, 'p':0, 'd':0, 'c':k**4},{'q':q, 'p':0, 'd':2, 'c':-2*k**2}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4y2sphlapl2_l1():
    """Spherical shell 4th integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':2, 'd':4, 'c':1}, {'q':4, 'p':1, 'd':3, 'c':4}, {'q':4, 'p':0, 'd':2, 'c':-4}]
    terms = symbolic.change_variable(terms, 'linear')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
