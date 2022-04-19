"""Module to compute the recurrence relations for the spectral Jacobi operators of the radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy

import quicc.recurrence.symbolic_jacobi as mod

l = sympy.Symbol('l')
n = sympy.Symbol('n')

w_alpha_any = sympy.Symbol('a')
w_alpha_legendre = 0
w_alpha_chebyshev = -sympy.Rational(1,2)

w_beta_any = sympy.Symbol('b')
w_beta_minhalf = l - sympy.Rational(1,2)
w_beta_zero = l
w_beta_plushalf = l + sympy.Rational(1,2)

w_layout_list = 0
w_layout_py = 1
w_layout_cpp = 2

w_alpha = w_alpha_chebyshev
w_beta = w_beta_minhalf
w_layout = w_layout_list

symbolic = mod.SymbolicJacobi(a = w_alpha, b = w_beta)

def showDiags(r, layout, lshift = 0, name = ""):

    print("="*30 + " " + name + " " + "="*30 + "\n")

    # Simple list of diagonals
    if layout == w_layout_list:
        # Print recurrence relation per diagonals
        for k,rec in sorted(r.items()):
            print("\t" + str(k) + ": \t" + str(rec))
    # Python functions
    elif layout == w_layout_py:
        import re
        t = 4
        # Generate Python function for diagonals
        for k,rec in sorted(r.items()):
            if k < 0:
                sk = "_"+str(abs(k))
                sd = str(abs(k))+". sub"
            elif k > 0:
                sk = str(k)
                sd = str(k)+". super"
            else:
                sk = str(k)
                sd = "main "
            print(" "*t + "# Generate " + sd + "diagonal")
            print(" "*t + "def d" + sk + "(n):")
            code_rec = re.sub(r'(\d+)', r'\1.0', str(rec))
            code_rec = re.sub(r'\*\*(\d+)\.0', r'**\1', code_rec)
            print(" "*2*t + "val = " + str(code_rec))
            if lshift == 0:
                sl = ""
            else:
                sl = ", " + str(lshift)
            print(" "*2*t + "return normalize_row(n, l, " + str(k) + sl +")*val"+"\n")
    # C++ functions
    elif layout == w_layout_cpp:
        import re
        t = 3
        # Generate C++ function for diagonals
        for k,rec in sorted(r.items()):
            if k < 0:
                sk = "_"+str(abs(k))
                sd = str(abs(k))+". sub"
            elif k > 0:
                sk = str(k)
                sd = str(k)+". super"
            else:
                sk = str(k)
                sd = "main "

            code_rec = re.sub(r'(\d+)', r'\1.0', str(rec))
            code_rec = re.sub(r'\*\*(\d+)\.0', r'**\1', code_rec)
            code_rec = re.sub(r'\bl\b\*\*(\d+)', r'l\1', code_rec)
            powers = list([])
            for p in range(2, 10):
                if re.search(r'\bn\b\*\*'+str(p), code_rec) is not None:
                    powers.append(p)
            code_rec = re.sub(r'\bn\b\*\*(\d+)', r'n.pow(\1)', code_rec)
            code_rec = re.sub(r'\bl\b', r'l1', code_rec)
            print(" "*t + "ACoeff " + name + "Diags::d" + sk +"(const ACoeff& n) const")
            print(" "*t + "{")
            print(" "*2*t + "MHDFloat l1 = this->l();")
            for p in powers:
                print(" "*2*t + "MHDFloat l" + str(p) + " = std::pow(l1, " + str(p) + ");")
            print(" "*2*t + "ACoeff val;" + "\n")
            print(" "*2*t + "val = " + code_rec + ";" + "\n")
            print(" "*2*t + "return this->normalizeDiag(n, " + str(k) +")*val;")
            print(" "*t + "}" + "\n")
    print("-"*70)
    print("\n")


def chooseParameters(alpha, beta, layout = 0):
    """Choose the Jacobi polynomial to use"""

    global symbolic, w_alpha, w_beta, w_layout
    w_alpha = alpha
    w_beta = beta
    w_layout = layout
    symbolic = mod.SymbolicJacobi(a = w_alpha, b = w_beta)

def x1():
    """Sphere x operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "X1")

def r2():
    """Sphere r^2 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,2)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "R2")

def r4():
    """Sphere r^4 operator"""

    # Setup terms in recurrence
    terms = [{'q':0, 'p':2, 'd':0, 'c':sympy.Rational(1,4)},{'q':0, 'p':1, 'd':0, 'c':sympy.Rational(1,2)},{'q':0, 'p':0, 'd':0, 'c':sympy.Rational(1,4)}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "R4")

def i1():
    """Sphere i1 (i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':1, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I1")

def i2():
    """Sphere i2 (i1r1i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I2")

def i2divrdiff():
    """Sphere i2 ((i1r1)^2 1/r D) operator"""

    # Setup terms in recurrence
    terms = [{'q':2, 'p':0, 'd':1, 'c':4.0}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I2R_1D1")

def i2lapl():
    """Sphere i2lapl (i1r1i1r1 lapl) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':2, 'p':1, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':2, 'c':8},
            {'q':2, 'p':0, 'd':1, 'c':4*(2*l + 3)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I2Lapl")

def i4():
    """Sphere i4 (i1r1i1r1i1r1i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I4")

def i4divrdiff():
    """Sphere i4 ((i1r1)^4 1/r D) operator"""

    # Setup terms in recurrence
    terms = [{'q':4, 'p':0, 'd':1, 'c':4.0}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I4R_1D1")

def i4lapl():
    """Sphere i4lapl (i1r1i1r1i1r1i1r1 lapl) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':1, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':2, 'c':8},
            {'q':4, 'p':0, 'd':1, 'c':4*(2*l + 3)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I4Lapl")

def i4lapl2():
    """Sphere i4lapl2 (i1r1i1r1i1r1i1r1 bilapl) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':4, 'p':2, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':4, 'c':128},
            {'q':4, 'p':0, 'd':4, 'c':64},
            {'q':4, 'p':1, 'd':3, 'c':64*(2*l+5)},
            {'q':4, 'p':0, 'd':3, 'c':64*(2*l+5)},
            {'q':4, 'p':0, 'd':2, 'c':16*(2*l+3)*(2*l+5)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I4Lapl2")

def i6():
    """Sphere i6 (i1r1i1r1i1r1i1r1i1r1i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':6, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I6")

def i1qm():
    """Sphere i1qm (i1r1 coriolis Q(l-1)) operator"""

    # Compute starting terms
    r = symbolic.spectral_increase({0:-4}, True)

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = -1, name = "I1Qm")

def i1qp():
    """Sphere i2qp (i1r1 coriolis Q(l+1)) operator"""

    partA = symbolic.spectral_integral_decrease({0:(2*l+1)}, True)
    partB = symbolic.spectral_decrease({0:2}, True)

    r = partA
    for k,v in partB.items():
        r[k] = r[k] + v
        r[k] = r[k].simplify().factor()

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = 1, name = "I1Qp")

def i2qm():
    """Sphere i2qm (i1r1i1r1 coriolis Q(l-1)) operator"""

    # Compute starting terms
    fs = symbolic.spectral_increase({0:-4}, False)

    # Setup terms in recurrence
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':1},
            ]
    r = symbolic.build_recurrence(terms, fs)

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = -1, name = "I2Qm")

def i2qp():
    """Sphere i2qp (i1r1i1r1 coriolis Q(l+1)) operator"""

    above = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)

    # Setup terms in recurrence
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':(2*l + 1)},
            ]
    tmp = above.build_recurrence(terms, {0:1}, False)
    partA = symbolic.spectral_decrease(tmp, True)

    # Compute starting terms
    fs = symbolic.spectral_decrease({0:-(2*l-1)}, False)

    # Setup terms in recurrence
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':1},
            ]
    partB = symbolic.build_recurrence(terms, fs, True)
    r = partA
    for k,v in partB.items():
        r[k] = r[k] + v
        r[k] = r[k].simplify().factor()

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = 1, name = "I2Qp")

def i4qm():
    """Sphere i4qm (i1r1i1r1i1r1i1r1 coriolis Q(l-1)) operator"""

    # Compute starting terms
    fs = symbolic.spectral_increase({0:-4}, False)

    # Setup terms in recurrence
    terms = [
            {'q':3, 'p':0, 'd':0, 'c':1},
            ]
    r = symbolic.build_recurrence(terms, fs)

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = -1, name = "I4Qm")

def i4qp():
    """Sphere i4qp (i1r1i1r1i1r1i1r1 coriolis Q(l+1)) operator"""

    above = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)

    # Setup terms in recurrence
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':(2*l + 1)},
            ]
    tmp = above.build_recurrence(terms, {0:1}, False)
    tmp2 = symbolic.spectral_decrease(tmp, False)
    terms = [
            {'q':2, 'p':0, 'd':0, 'c':1},
            ]
    partA = symbolic.build_recurrence(terms, tmp2, True)

    # Compute starting terms
    fs = symbolic.spectral_decrease({0:-(2*l-1)}, False)

    # Setup terms in recurrence
    terms = [
            {'q':3, 'p':0, 'd':0, 'c':1},
            ]
    partB = symbolic.build_recurrence(terms, fs, True)
    r = partA
    for k,v in partB.items():
        r[k] = r[k] + v
        r[k] = r[k].simplify().factor()

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = 1, name = "I4Qp")
