"""Module to compute the recurrence relations for the spectral Jacobi operators of the radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_jacobi as mod

m = sympy.Symbol('m')
n = sympy.Symbol('n')

symbolic = mod.SymbolicJacobi(a = -sympy.Rational(1,2), b = m - sympy.Rational(1,2))

def chooseParameters(alpha, beta):
    """Choose the Jacobi polynomial to use"""

    global symbolic
    symbolic = mod.SymbolicJacobi(a = alpha, b = beta)

def useLegendreJacobi():
    """Setup Legendre type Jacobi polynomial"""

    chooseParameters(0, m)

def useChebyshevZeroJacobi():
    """Setup Chebyshev type Jacobi polynomial for m = 0"""

    chooseParameters(-sympy.Rational(1,2), -sympy.Rational(1,2))

def useChebyshevJacobi():
    """Setup Chebyshev type Jacobi polynomial"""

    chooseParameters(-sympy.Rational(1,2), m - sympy.Rational(1,2))

def x1():
    """Cylinder x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r2():
    """Cylinder r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,2)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r4():
    """Cylinder r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':sympy.Rational(1,4)},{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def r6():
    """Cylinder r^6 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':3, 'd':0, 'c':sympy.Rational(1,8)},{'q':0, 'p':2, 'd':0, 'c':sympy.Rational(3,8)},{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(3,8)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,8)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1():
    """Cylinder i1 (i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2():
    """Cylinder i2 ((i1r1)^2) operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2laplh():
    """Cylinder i2laplh ((i1r1)^2 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':2, 'p':1, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4():
    """Cylinder i4 ((i1r1)^4) operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4divrdiff():
    """Cylinder i4 ((i1r1)^4 1/r D) operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':1, 'c':4.0}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4laplh():
    """Cylinder i4laplh ((i1r1)^4 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':1, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4lapl2h():
    """Cylinder i4lapl2h ((i1r1)^4 bilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':2, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':4, 'c':128},
            {'q':4, 'p':0, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':3, 'c':128*(m+2)},
            {'q':4, 'p':0, 'd':3, 'c':128*(m+2)},
            {'q':4, 'p':0, 'd':2, 'c':64*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6():
    """Cylinder i6 ((i1r1)^6) operator"""

    # Setup terms in recurrence
    terms = [{'q':6, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6divrdiff():
    """Cylinder i6 ((i1r1)^6 1/r D) operator"""

    # Setup terms in recurrence
    terms = [{'q':6, 'p':0, 'd':1, 'c':4.0}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6laplh():
    """Cylinder i6laplh ((i1r1)^4 laplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':1, 'd':2, 'c':8},
            {'q':6, 'p':0, 'd':2, 'c':8},
            {'q':6, 'p':0, 'd':1, 'c':8*(m + 1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6lapl2h():
    """Cylinder i6lapl2 ((i1r1)^6 bilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':2, 'd':4, 'c':64},
            {'q':6, 'p':1, 'd':4, 'c':128},
            {'q':6, 'p':0, 'd':4, 'c':64},
            {'q':6, 'p':1, 'd':3, 'c':128*(m+2)},
            {'q':6, 'p':0, 'd':3, 'c':128*(m+2)},
            {'q':6, 'p':0, 'd':2, 'c':64*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6lapl3h():
    """Cylinder i6lapl3 ((i1r1)^6 trilaplh) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':6, 'p':3, 'd':6, 'c':512},
            {'q':6, 'p':2, 'd':6, 'c':1536},
            {'q':6, 'p':1, 'd':6, 'c':1536},
            {'q':6, 'p':0, 'd':6, 'c':512},
            {'q':6, 'p':2, 'd':5, 'c':1536*(m+3)},
            {'q':6, 'p':1, 'd':5, 'c':3072*(m+3)},
            {'q':6, 'p':0, 'd':5, 'c':1536*(m+3)},
            {'q':6, 'p':1, 'd':4, 'c':1536*(m+3)*(m+2)},
            {'q':6, 'p':0, 'd':4, 'c':1536*(m+3)*(m+2)},
            {'q':6, 'p':0, 'd':3, 'c':512*(m+3)*(m+2)*(m+1)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1dr():
    """Cylinder 1st integral operator of radial derivative (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), 1:(2*n+1)/(n+1)})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2dr():
    """Cylinder 2nd integral operator of radial derivative (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), -1:(2*n-1)/n})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i4dr():
    """Cylinder 4th integral operator of radial derivative (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':3, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), -1:(2*n-1)/n})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i6dr():
    """Cylinder 6th integral operator of radial derivative (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':5, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), -1:(2*n-1)/n})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i1r_1dr():
    """Cylinder 1st integral operator of 1/r D r (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), -1:4*n/(2*n-1)})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")

def i2r_1dr():
    """Cylinder 1st integral operator of 1/r D r (projects accross beta)"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    n = sympy.Symbol('n')
    r = symbolic.build_recurrence(terms, {0:sympy.Rational(2,1), 1:4*(n+1)/(2*n+1)})

    # Print recurrence relation per diagonals
    for k,rec in sorted(r.items()):
        print("\t" + str(k) + ": \t" + str(rec))
    print("\n")
