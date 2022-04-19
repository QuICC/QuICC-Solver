"""Module to compute the recurrence relations for the spectral operators of the radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def x1():
    """Cylinder multiplication by x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def x2():
    """Cylinder multiplication by x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Cylinder 1st integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1d1():
    """Cylinder 1st integral of x 1st derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1div():
    """Cylinder 1st integral of x radial divergence operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':1, 'c':1}, {'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1():
    """Cylinder 1st integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x2():
    """Cylinder 1st integral of x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cylinder 2nd integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x1():
    """Cylinder 2nd integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2d1():
    """Cylinder 2nd integral of x^2 1st derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3d1():
    """Cylinder 2nd integral of x^3 1st derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':3, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3d1x_2():
    """Cylinder 2nd integral of x^3 1st derivative 1/x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-2.0}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2d2():
    """Cylinder 2nd integral of x^2 2nd derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2():
    """Cylinder 2nd integral of x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3():
    """Cylinder 2nd integral of x^3 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':3, 'd':0, 'c':1}]
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
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2laplh():
    """Cylinder 2nd integral of x^2 horizontal laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-m**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x2vlaplh():
    """Cylinder 2nd integral of x^2 horizontal vector laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':1}, {'q':2, 'p':0, 'd':0, 'c':-(m**2 + 1)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x3vlaplhx_1():
    """Cylinder 2nd integral of x^3 horizontal laplacian 1/x operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':-1}, {'q':2, 'p':0, 'd':0, 'c':-m**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cylinder 4th integral operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4():
    """Cylinder 4th integral of x^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4laplh():
    """Cylinder 4th integral of x^4 horizontal laplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':1}, {'q':4, 'p':2, 'd':0, 'c':-m**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4x4lapl2h():
    """Cylinder 4th integral of x^4 horizontal bilaplacian operator"""

    # Setup terms in recurrence
    m = sympy.Symbol('m')
    terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':2}, {'q':4, 'p':2, 'd':2, 'c':-(1+2*m**2)}, {'q':4, 'p':1, 'd':1, 'c':(1+2*m**2)}, {'q':4, 'p':0, 'd':0, 'c':(m**2 - 4)*m**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
