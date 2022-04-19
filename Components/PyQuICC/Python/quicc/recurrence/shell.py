"""Module to compute the recurrence relations for the spectral operators of the radial direction in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def r1():
    """Spherical shell r operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r2():
    """Spherical shell r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r4():
    """Spherical shell r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':4, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1r1():
    """Spherical shell 1st integral r operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Spherical shell 1st integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Spherical shell 2nd integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2d1():
    """Spherical shell 2nd integral operator of 1 derivative"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r1():
    """Spherical shell 2nd integral of r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2():
    """Spherical shell 2nd integral of r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r3():
    """Spherical shell 2nd integral of r^3 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':3, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r1d1r1():
    """Spherical shell 2nd integral of r D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}, {'q':2, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2d1():
    """Spherical shell 2nd integral of r^2 D operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2lapl():
    """Spherical shell 2nd integral of r^2 laplacianoperator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':2}, {'q':2, 'p':0, 'd':0, 'c':-l*(l+1)}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r3lapl():
    """Spherical shell 2nd integral of r^3 laplacianoperator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':3, 'd':2, 'c':1}, {'q':2, 'p':2, 'd':1, 'c':2}, {'q':2, 'p':1, 'd':0, 'c':-l*(l+1)}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i3():
    """Spherical shell 3rd integral operator"""

    # Setup terms in recurrence
    terms = [{'q':3, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Spherical shell 4th integral operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r1():
    """Spherical shell 4th integral of r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':1, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r3d2():
    """Spherical shell 4th integral of r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':3, 'd':2, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r1d1r1():
    """Spherical shell 4th integral of r D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':1, 'd':0, 'c':1}, {'q':4, 'p':2, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4laplrd1r1():
    """Spherical shell 4th integral of r^4 laplacian D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':3, 'c':1}, {'q':4, 'p':3, 'd':2, 'c':3}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r3d1r1():
    """Spherical shell 4th integral of r^3 D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':1, 'c':1}, {'q':4, 'p':3, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r2():
    """Spherical shell 4th integral of r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':2, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r3():
    """Spherical shell 4th integral of r^3 operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':3, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4():
    """Spherical shell 4th integral of r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4d1():
    """Spherical shell 4th integral of r^4 D operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4lapl():
    """Spherical shell 4th integral of r^4 laplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':2}, {'q':4, 'p':2, 'd':0, 'c':-l*(l+1)}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r2lapl2_l1():
    """Spherical shell 4th integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':2, 'd':4, 'c':1}, {'q':4, 'p':1, 'd':3, 'c':4}, {'q':4, 'p':0, 'd':2, 'c':-4}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4lapl2():
    """Spherical shell 4th integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':4}, {'q':4, 'p':2, 'd':2, 'c':-2*l*(l+1)}, {'q':4, 'p':0, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4d1():
    """Spherical shell 4th integral operator of 1 derivative"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':1, 'c':1}]
    terms = symbolic.change_variable(terms, 'linear_r2x')
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
