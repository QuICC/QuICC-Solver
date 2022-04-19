"""Module to compute the recurrence relations for the spectral operators of the radial direction in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def x1():
    """Cylindrical annulus x multiplication operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def x2():
    """Cylindrical annulus x^2 multiplication operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Cylindrical annulus first integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1d1():
    """Cylindrical annulus first integral of x of first derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1div():
    """Cylindrical annulus 1st integral of x radial divergence operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':1, 'p':1, 'd':1, 'c':1}, {'q':1, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1():
    """Cylindrical annulus first integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x2():
    """Cylindrical annulus first integral of x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cylindrical annulus second integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x1():
    """Cylindrical annulus second integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2d1():
    """Cylindrical annulus second integral of x^2 of first derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3d1():
    """Cylindrical annulus second integral of x^3 of first derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':3, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3d1x_2():
    """Cylindrical annulus second integral of x^3 of first derivative 1/x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-2}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2d2():
    """Cylindrical annulus second integral of x^2 of second derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2():
    """Cylindrical annulus second integral of x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3():
    """Cylindrical annulus second integral of x^3 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':3, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2div():
    """Cylindrical annulus second integral of x^2 radial divergence operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}, {'q':2, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2laplh():
    """Cylindrical annulus second integral of x^2 horizontal laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-m**2}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2vlaplh():
    """Cylindrical annulus second integral of x^2 horizontal vector laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-(m**2 + 1)}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3vlaplhx_1():
    """Cylindrical annulus second integral of x^3 horizontal vector laplacian operator of f/x"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':-1}, {'q':2, 'p':0, 'd':0, 'c':-m**2}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cylindrical annulus fourth integral of x^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4():
    """Cylindrical annulus fourth integral of x^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4laplh():
    """Cylindrical annulus fourth integral of x^4 horizontal laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':1}, {'q':4, 'p':2, 'd':0, 'c':-m**2}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4lapl2h():
    """Cylindrical annulus fourth integral of x^4 horizontal bilaplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':2}, {'q':4, 'p':2, 'd':2, 'c':-(1+2*m**2)}, {'q':4, 'p':1, 'd':1, 'c':(1+2*m**2)}, {'q':4, 'p':0, 'd':0, 'c':(m**2 - 4)*m**2}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
