"""Module to compute the recurrence relations for the spectral operators of the a cartesian direction"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_chebyshev as mod

symbolic = mod.SymbolicChebyshev()

def z1():
    """Cartesian multiplication by z operator"""

    # Setup terms in recurrence
    cs = sympy.symbols('cs')
    terms = [{'q':0, 'p':0, 'd':0, 'c':1/cs},{'q':0, 'p':1, 'd':0, 'c':1/cs}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def x1():
    """Cartesian multiplication by x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def x2():
    """Cartesian multiplication by x^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def x4():
    """Cartesian multiplication by x^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':4, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Cartesian first integral operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1x1():
    """Cartesian first integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1z1():
    """Cartesian first integral of z operator"""

    # Setup terms in recurrence
    cs = sympy.symbols('cs')
    terms = [{'q':1, 'p':1, 'd':0, 'c':1/cs},{'q':1, 'p':0, 'd':0, 'c':1/cs}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cartesian second integral operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2z1():
    """Cartesian 2nd integral of z operator"""

    # Setup terms in recurrence
    cs = sympy.symbols('cs')
    terms = [{'q':2, 'p':1, 'd':0, 'c':1/cs},{'q':2, 'p':0, 'd':0, 'c':1/cs}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2x1():
    """Cartesian second integral of x operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2d1():
    """Cartesian second integral of first derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2lapl():
    """Cartesian second integral of laplacian operator"""

    # Setup terms in recurrence
    k, l, cs = sympy.symbols('k l cs')
    terms = [{'q':2, 'p':0, 'd':2, 'c':cs**2},{'q':2, 'p':0, 'd':0, 'c':-k**2},{'q':2, 'p':0, 'd':0, 'c':-l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2laplh():
    """Cartesian second integral of horizontal laplacian operator"""

    # Setup terms in recurrence
    k, cs = sympy.symbols('k cs')
    terms = [{'q':2, 'p':0, 'd':2, 'c':cs**2},{'q':2, 'p':0, 'd':0, 'c':-k**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i3():
    """Cartesian third integral operator"""

    # Setup terms in recurrence
    terms = [{'q':3, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cartesian fourth integral operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4d1():
    """Cartesian 4th integral of 1st derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':1, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4d2():
    """Cartesian 4th integral of 2nd derivative operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':2, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl():
    """Cartesian fourth integral of laplacian operator"""

    # Setup terms in recurrence
    k, l, cs = sympy.symbols('k l cs')
    terms = [{'q':4, 'p':0, 'd':2, 'c':cs**2},{'q':4, 'p':0, 'd':0, 'c':-k**2},{'q':4, 'p':0, 'd':0, 'c':-l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4laplh():
    """Cartesian fourth integral of horizontal laplacian operator"""

    # Setup terms in recurrence
    k, cs = sympy.symbols('k cs')
    terms = [{'q':4, 'p':0, 'd':2, 'c':cs**2},{'q':4, 'p':0, 'd':0, 'c':-k**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2():
    """Cartesian fourth integral of bilaplacian operator"""

    # Setup terms in recurrence
    k, l, cs = sympy.symbols('k l cs')
    terms = [{'q':4, 'p':0, 'd':4, 'c':cs**4},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':0, 'c':l**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2*cs**2},{'q':4, 'p':0, 'd':2, 'c':-2*l**2*cs**2},{'q':4, 'p':0, 'd':0, 'c':2*k**2*l**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2h():
    """Cartesian fourth integral of horizontal bilaplacian operator"""

    # Setup terms in recurrence
    k, cs = sympy.symbols('k cs')
    terms = [{'q':4, 'p':0, 'd':4, 'c':cs**4},{'q':4, 'p':0, 'd':0, 'c':k**4},{'q':4, 'p':0, 'd':2, 'c':-2*k**2*cs**2}]
    r = symbolic.build_recurrence(terms, {0:1})
    n = sympy.Symbol('n')

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
