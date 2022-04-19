"""Module to compute the recurrence relations for the spectral operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def r1():
    """Sphere r operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r2():
    """Sphere r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Sphere 1st integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Sphere 2nd integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2():
    """Sphere second integral of r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r1d1r1():
    """Sphere 2nd integral of r D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}, {'q':2, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2lapl():
    """Sphere second integral of r^2 laplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':2, 'c':1}, {'q':2, 'p':1, 'd':1, 'c':2}, {'q':2, 'p':0, 'd':0, 'c':-l*(l+1)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r2():
    """Sphere fourth integral of r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4():
    """Sphere fourth integral of r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':4, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r3d1r1():
    """Sphere 4th integral of r^3 D r operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':1, 'c':1}, {'q':4, 'p':3, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4lapl():
    """Sphere fourth integral of r^4 laplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':2, 'c':1}, {'q':4, 'p':3, 'd':1, 'c':2}, {'q':4, 'p':2, 'd':0, 'c':-l*(l+1)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r2lapl2_l1():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':2, 'd':4, 'c':1}, {'q':4, 'p':1, 'd':3, 'c':4}, {'q':4, 'p':0, 'd':2, 'c':-4}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4lapl2():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':4, 'c':1}, {'q':4, 'p':3, 'd':3, 'c':4}, {'q':4, 'p':2, 'd':2, 'c':-2*l*(l+1)}, {'q':4, 'p':0, 'd':0, 'c':(l-1)*l*(l+1)*(l+2)}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r1():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r3():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':3, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r2d1():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':2, 'p':2, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4r4d1():
    """Sphere fourth integral of r^4 bilaplacian operator"""

    # Setup terms in recurrence
    l = sympy.Symbol('l')
    terms = [{'q':4, 'p':4, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
