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

def i3():
    """Sphere i3 (i1r1i1r1i1r1) operator"""

    # Setup terms in recurrence
    terms = [{'q':3, 'p':0, 'd':0, 'c':1}]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I3")

def i3lapl():
    """Sphere i3lapl (i1r1i1r1i1r1 lapl) operator"""

    # Setup terms in recurrence
    terms = [
            {'q':3, 'p':1, 'd':2, 'c':8},
            {'q':3, 'p':0, 'd':2, 'c':8},
            {'q':3, 'p':0, 'd':1, 'c':4*(2*l + 3)}
            ]
    r = symbolic.build_recurrence(terms, {0:1})

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, name = "I4Lapl")

def i3qm():
    """Sphere i3qm (i1r1i1r1i1r1 coriolis Q(l-1)) operator"""

    # Compute starting terms
    fs = symbolic.spectral_increase({0:-4}, False)

    # Setup terms in recurrence
    terms = [
            {'q':2, 'p':0, 'd':0, 'c':1},
            ]
    r = symbolic.build_recurrence(terms, fs)

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = -1, name = "I3Qm")

def i3qp():
    """Sphere i3qp (i1r1i1r1i2r1 coriolis Q(l+1)) operator"""

    above = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)

    # Setup terms in recurrence
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':(2*l + 1)},
            ]
    tmp = above.build_recurrence(terms, {0:1}, False)
    tmp2 = symbolic.spectral_decrease(tmp, False)
    terms = [
            {'q':1, 'p':0, 'd':0, 'c':1},
            ]
    partA = symbolic.build_recurrence(terms, tmp2, True)

    # Compute starting terms
    fs = symbolic.spectral_decrease({0:-(2*l-1)}, False)

    # Setup terms in recurrence
    terms = [
            {'q':2, 'p':0, 'd':0, 'c':1},
            ]
    partB = symbolic.build_recurrence(terms, fs, True)
    r = partA
    for k,v in partB.items():
        r[k] = r[k] + v
        r[k] = r[k].simplify().factor()

    # Print recurrence relation per diagonals
    showDiags(r, w_layout, lshift = 1, name = "I3Qp")

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
    """Sphere i1qp (i1r1 coriolis Q(l+1)) operator"""

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





############################################################
###    EXTRA MATERIAL FOR SPHEROID-RELATED OPERATORS     ###
############################################################


## HIGHER RAISE/LOWER OPERATORS: FIRST LAYER OF COEFFICIENTS
## MUST STILL WORK OUT HARMONISING WITH PM'S SIGN CONVENTION


# Some helpers
from sympy.abc import x
oneHalf = sympy.Rational(1, 2)
lpno2 = lambda n : l + oneHalf*n
lp1h  = lpno2(1) 


# Collate operator info

def get_dnp_info(mu):
    op_terms  = [(x + lp1h + m) for m in range(1, mu+1)]
    prefactor = 2**mu
    return op_terms, prefactor

def get_lapldnp_info(mu):
    op_terms  = [(x + lp1h + m) for m in range(0, mu+1)]
    prefactor = 2**(mu+3)
    return op_terms, prefactor

def get_lapldnm_info(mu):
    s = abs(mu)
    op_terms  = [(x + lp1h - int(s))]
    prefactor = 2**(2*s+3)
    return op_terms, prefactor

# For now this is fine...
def get_op_info(mu, op_type):
    if op_type == 'dnp':
        return get_dnp_info(mu)
    if op_type == 'lapldnp':
        return get_lapldnp_info(mu)
    if op_type == 'lapldnm':
        return get_lapldnm_info(mu)
    raise TypeError(f'op_type {op_type} is wrong')



# Compute gamma for getting 
# 'YD' in terms of 'DY'

gamma_vec_zero = np.array([1]).astype(np.int32)
# gamma_vec_one  = np.array([-1, 1]).astype(np.int32)

def gamma_vec_j_from_jm1(g_vec_jm1):
    """
    g^{(j)} components from those of g^{(j-1)}
    Just as written in eq. (XX) of paper...
    """

    j = len(g_vec_jm1)  # g^{(j-1)} has length j

    # Initialise g^{(j)} (length j+1 !)
    # and set components
    g_vec_j     = np.empty(j+1)
    g_vec_j[ 0] = -1 * g_vec_jm1[ 0]
    g_vec_j[-1] =  1 * g_vec_jm1[-1]
    for p in range(1, j):
        g_vec_j[p] = g_vec_jm1[p-1] - (p+1)*g_vec_jm1[p]
    
    return g_vec_j.astype(np.int32)

def gamma_vec(j: int):
    """
    Recursively compute the coefficients
    Start from j=1:
        g_0^(1) = -1
        g_1^(1) = +1
    as in eq. (XX) of paper
    """
    if j<0: 
        raise ValueError(f"j too low in gamma-matrix computation: j={j}")
    
    # NOTE CHECK!!!!!
    if j == 0:
        return gamma_vec_zero
    # if j == 1: 
    #     return gamma_vec_one

    return gamma_vec_j_from_jm1(gamma_vec(j-1))




# Make the corresponding polynomial by:
# making a tuple containing its 
# constituents; taking the sympy product 
# of that tuple; then getting the coefficients
def get_op_coeffs(op_terms, prefac):
    # get the corresponding betas WITH the
    # op-dependent prefactor multplied in
    expr = sympy.Poly(sympy.prod(op_terms), x)
    beta = sympy.factor(expr.coeffs())
    beta = tuple(reversed(beta)) # to match my indexing
    beta = [prefac * b for b in beta]
    return sympy.factor(beta) # factor again just in case


def a_from_gb(beta, bigN):
    """
    get alpha from gamma and beta
    return it as a dict in the right form
    """
    alpha = []
    for i in range(bigN+1):
        tmp = 0
        for j in range(i, bigN+1):
            tmp += gamma_vec(j)[i] * beta[j]
        alpha.append(tmp)
    return dict(enumerate(sympy.factor(alpha)))


def op_coeff(mu, op_type):
    terms, prefac = get_op_info(mu, op_type); bigN = len(terms)
    beta = get_op_coeffs(terms, prefac)
    alpha = a_from_gb(beta, bigN)
    return alpha



# PLUS

def D2pCoeffs():
    return op_coeff(2, 'dnp')
def D3pCoeffs():
    return op_coeff(3, 'dnp')
def D4pCoeffs():
    return op_coeff(4, 'dnp')

def laplDpCoeffs():
    return op_coeff(1, 'lapldnp')
def laplD2pCoeffs():
    return op_coeff(2, 'lapldnp')

def laplDmCoeffs():
    return op_coeff(-1, 'lapldnm') # CHANGING FIRST ARG DOESN'T SEEM TO CHANGE ANYTHING 🤔🤔🤔
def laplD2mCoeffs():
    return op_coeff(-2, 'lapldnm')


# MINUS (NOTE improve labelling of coefficients)
def DnMFactor(absmu):
    return 2**(2*absmu)
# D^{-2}
def D2mCoeffs():
    """Done to match my 'h' rather than 'g'"""
    return {0: DnMFactor(2)}
# D^{-3}
def D3mCoeffs():
    """Done to match my 'h' rather than 'g'"""
    return {0: DnMFactor(3)}
# D^{-4}
def D4mCoeffs():
    """Done to match my 'h' rather than 'g'"""
    return {0: DnMFactor(4)}



## RECURRENCE STUFF

# HELPERS

def parts_from_coeffs(coeffs):
    """
    Sort dictionary of coefficients
    into a suitable form
    """
    return [{0:coeffs[k]} for k in sorted(coeffs)]

def sum_parts(parts):
    # merge dicts, simplify+factor each poly
    r = parts[0]
    for part in parts[1:]:
        for k, v in part.items():
            r[k] = r[k] + v
            # r[k] = r.get(k, 0) + v
            r[k] = r[k].simplify().factor()
    return r

def integrate_same_l(expr_in, n_int, finish=True):
    terms = [
            {'q':n_int, 'p':0, 'd':0, 'c':1}
            ]
    expr_out = symbolic.build_recurrence(terms, expr_in, finish)

    return expr_out

def net_integral_same_l(exprs_in, n_int, 
                        finish=True,
                        sumup=True,
                        ):
    """
    Takes a list of expressions and integrates each
    member n_int times
    """
    rr = []
    for expr in exprs_in:
        rr.append(
            integrate_same_l(
                expr, 
                n_int, 
                finish=finish
            )
        )
    
    return sum_parts(rr) if sumup else rr 

def simplify_single(expr):
    r = expr 
    for k, v in expr.items():
        r[k] = v.simplify().factor()
    return r 

def qi_op_from_gen_d_op(d_op_generic, n_integ):

    if n_integ == 0:
        # Get the generic parts and sum
        # (No integral so get asrow)
        l_shift, thing = d_op_generic(asrow=True)
        r = sum_parts(thing)
    else:
        # Get the generic parts, do
        # n_integ net intgrals, and sum
        l_shift, thing = d_op_generic()
        r = net_integral_same_l(thing, n_integ, sumup=True)   

    return r, l_shift

def get_diags(d_op_generic, n_integ, *, name=None):
    """
    Compute the algebraic expressions, then show them
    and pass the right name
    (Essentially wraps PM's original showDiags)
    """

    # Compute the actual algebraic expressions
    r, l_shift = qi_op_from_gen_d_op(d_op_generic, n_integ)

    # Get the name of the function 
    # (to be fed to PM's original showDiags)
    if name is None:
        f = sys._getframe(1)           # caller's frame
        try:
            name = f.f_code.co_name     # e.g., "i3D3p"
        finally:
            del f                       # break ref cycle
    
    # Show the diags
    showDiags(r, w_layout, lshift=l_shift, name=name)

def get_offset_level(s): 
    return mod.SymbolicJacobi(
        a = w_alpha, 
        b = w_beta + s
    ) 


# Super generic (?) layer

def raise_l_plus_p_to_l_plus_q(parts, 
                               p=0, q=0, 
                               asrow=False):
    
    l_offsets = range(p+1, q+1)

    exprs_out = []
    for expr in parts:
        raised = expr
        for s in l_offsets:
            next_level = get_offset_level(s)

            finish = (s==q)
            ar = asrow if finish else False 

            raised = next_level.spectral_increase(raised, asrow=ar)

        exprs_out.append(raised)
        
    return exprs_out


 

def first_order_sum_terms(parts, asrow=False):
    """
    (a_1 Y + a_0 I) acting on J^{l+1}
    """

    # (1+x) P^{l+1}
    ## Get (1+x) P^{l+1} in terms of P^{l}
    partA = symbolic.spectral_decrease(parts[1], asrow)

    # I P^{l+1}
    ## Get I P^{l+1} in terms of P^{l}
    partB = symbolic.spectral_integral_decrease(parts[0], asrow)

    return [partA, partB]

def second_order_sum_terms(parts, asrow=False):
    """
    (a_2 Y^2 + a_1 I Y + a_0 I^2) acting on J^{l+2}
    """

    oneUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)

    # (1+x)^{2} P^{l+2}
    ## Get (1+x) P^{l+2} in terms of P^{l+1}
    partA = oneUp.spectral_decrease(parts[2], False)
    ## Get (1+x) P^{l+1} in terms of P^{l}
    partA = symbolic.spectral_decrease(partA, asrow)

    # I (1+x) P_{n}^{l+2}
    ## Get (1+x) P^{l+2} in terms of P^{l+1}
    partB = oneUp.spectral_decrease(parts[1], False)
    ## Get I P^{l+1} in terms of P^{l}
    partB = symbolic.spectral_integral_decrease(partB, asrow)

    # I I P_{n}^{l+2}
    ## Get I P^{l+2} in terms of P^{l+1}
    partC = oneUp.spectral_integral_decrease(parts[0], False)
    ## Get I P^{l+1} in terms of P^{l}
    partC = symbolic.spectral_integral_decrease(partC, asrow)

    return [partA, partB, partC]

def third_order_sum_terms(parts, asrow=False):
    """
    (a_3 Y^3 + a_2 I Y^2 + a_1 I^2 Y + a_0 I^3) acting on J^{l+3}
    """

    oneUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)
    twoUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 2)

    ## (1+x)^{3} P^{l+3} = (1+x) (1+x) (1+x) P^{l+3}
    # (1+x) P^{l+3} --> P^{l+2}
    partA = twoUp.spectral_decrease(parts[3], False)
    # (1+x) P^{l+2} --> P^{l+1}
    partA = oneUp.spectral_decrease(partA, False)
    # (1+x) P^{l+1} --> P^{l}
    partA = symbolic.spectral_decrease(partA, asrow)

    ## I (1+x)^{2} P^{l+3}
    # (1+x) P^{l+3} --> P^{l+2}
    partB = twoUp.spectral_decrease(parts[2], False)
    # (1+x) P^{l+2} --> P^{l+1}
    partB = oneUp.spectral_decrease(partB, False)
    # I P^{l+1} ------> P^{l}
    partB = symbolic.spectral_integral_decrease(partB, asrow)

    ## I^{2} (1+x) P_{n}^{l+3}
    # (1+x) P^{l+3} --> P^{l+2}
    partC = twoUp.spectral_decrease(parts[1], False)
    # I P^{l+2} ------> P^{l+1}
    partC = oneUp.spectral_integral_decrease(partC, False)
    # I P^{l+1} ------> P^{l}
    partC = symbolic.spectral_integral_decrease(partC, asrow)

    ## I3 P^{l+3}
    # I P^{l+3} ------> P^{l+2}
    partD = twoUp.spectral_integral_decrease(parts[0], False)
    # I P^{l+2} ------> P^{l+1}
    partD = oneUp.spectral_integral_decrease(partD, False)
    # I P^{l+1} ------> P^{l}
    partD = symbolic.spectral_integral_decrease(partD, asrow)

    return [partA, partB, partC, partD]

def fourth_order_sum_terms(parts, asrow=False):
    """
    (a_4 Y^4 + a_3 I Y^3 + a_2 I^2 Y^2 + a_1 I^3 Y + a_0 I^4) acting on J^{l+4}
    """

    oneUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 1)
    twoUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 2)
    threeUp = mod.SymbolicJacobi(a = w_alpha, b = w_beta + 3)

    ## (1+x)^{4} P^{l+4}
    # (1+x) P^{l+4} --> P^{l+3}
    partA = threeUp.spectral_decrease(parts[4], False)
    # (1+x) P^{l+3} --> P^{l+2}
    partA = twoUp.spectral_decrease(partA, False)
    # (1+x) P^{l+2} --> P^{l+1}
    partA = oneUp.spectral_decrease(partA, False)
    # (1+x) P^{l+1} --> P^{l}
    partA = symbolic.spectral_decrease(partA, asrow)

    ## I (1+x)^{3} P^{l+4}
    # (1+x) P^{l+4} --> P^{l+3}
    partB = threeUp.spectral_decrease(parts[3], False)
    # (1+x) P^{l+3} --> P^{l+2}
    partB = twoUp.spectral_decrease(partB, False)
    # (1+x) P^{l+2} --> P^{l+1}
    partB = oneUp.spectral_decrease(partB, False)
    # I P^{l+1} ------> P^{l}
    partB = symbolic.spectral_integral_decrease(partB, asrow)

    ## I^{2} (1+x)^{2} P^{l+4}
    # (1+x) P^{l+4} --> P^{l+3}
    partC = threeUp.spectral_decrease(parts[2], False)
    # (1+x) P^{l+3} --> P^{l+2}
    partC = twoUp.spectral_decrease(partC, False)
    # I P^{l+2} ------> P^{l+1}
    partC = oneUp.spectral_integral_decrease(partC, False)
    # I P^{l+1} ------> P^{l}
    partC = symbolic.spectral_integral_decrease(partC, asrow)

    ## I^{3} (1+x) P^{l+4}
    # (1+x) P^{l+4} --> P^{l+3}
    partD = threeUp.spectral_decrease(parts[1], False)
    # I P^{l+3} ------> P^{l+2}
    partD = twoUp.spectral_integral_decrease(partD, False)
    # I P^{l+2} ------> P^{l+1}
    partD = oneUp.spectral_integral_decrease(partD, False)
    # I P^{l+1} ------> P^{l}
    partD = symbolic.spectral_integral_decrease(partD, asrow)

    ## I^{4} P^{l+4}
    # I P^{l+4} ------> P^{l+3}
    partE = threeUp.spectral_integral_decrease(parts[0], False)
    # I P^{l+3} ------> P^{l+2}
    partE = twoUp.spectral_integral_decrease(partE, False)
    # I P^{l+2} ------> P^{l+1}
    partE = oneUp.spectral_integral_decrease(partE, False)
    # I P^{l+1} ------> P^{l}
    partE = symbolic.spectral_integral_decrease(partE, asrow)

    return [partA, partB, partC, partD, partE]





## GENERIC 'PLUS' FORMS

def D2p_generic(asrow=False):

    l_shift = 2

    parts = parts_from_coeffs(D2pCoeffs())
    exprs = second_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs

def D3p_generic(asrow=False):

    l_shift = 3

    parts = parts_from_coeffs(D3pCoeffs())
    exprs = third_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs

def D4p_generic(asrow=False):

    l_shift = 4

    parts = parts_from_coeffs(D4pCoeffs())
    exprs = fourth_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs

def laplDp_generic(asrow=False):

    l_shift = 1

    parts = parts_from_coeffs(laplDpCoeffs())
    parts = raise_l_plus_p_to_l_plus_q(parts,
                                       p = 1, 
                                       q = 2)
    exprs = second_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs

def laplD2p_generic(asrow=False):

    l_shift = 2

    parts = parts_from_coeffs(laplD2pCoeffs())
    parts = raise_l_plus_p_to_l_plus_q(parts,
                                       p = 2, 
                                       q = 3)
    exprs = third_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs


## GENERIC 'MINUS' FORMS

def D2m_generic(asrow=False):

    l_shift = -2

    parts = parts_from_coeffs(D2mCoeffs())
    exprs = raise_l_plus_p_to_l_plus_q(parts,
                                       p = -2, 
                                       q =  0,  
                                       asrow=asrow)

    return l_shift, exprs

def D3m_generic(asrow=False):

    l_shift = -3

    parts = parts_from_coeffs(D3mCoeffs())
    exprs = raise_l_plus_p_to_l_plus_q(parts,
                                       p = -3, 
                                       q =  0,  
                                       asrow=asrow)

    return l_shift, exprs

def D4m_generic(asrow=False):

    l_shift = -4

    parts = parts_from_coeffs(D4mCoeffs())
    exprs = raise_l_plus_p_to_l_plus_q(parts,
                                       p = -4, 
                                       q =  0,  
                                       asrow=asrow)

    return l_shift, exprs

def laplDm_generic(asrow=False):

    l_shift = -1

    parts = parts_from_coeffs(laplDmCoeffs())
    parts = raise_l_plus_p_to_l_plus_q(parts,
                                       p = -1, 
                                       q =  1)
    exprs = first_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs

def laplD2m_generic(asrow=False):

    l_shift = -2

    parts = parts_from_coeffs(laplD2mCoeffs())
    parts = raise_l_plus_p_to_l_plus_q(parts,
                                       p = -2, 
                                       q =  1)
    exprs = first_order_sum_terms(parts, asrow=asrow)

    return l_shift, exprs


## ACTUAL NEW OPERATORS

def i3D2p():
    """
    I3 of D^{2} operator:
    Single integral of generic d2p
    """
    get_diags(D2p_generic, 1)

def i3D2m():
    """
    I3 of D^{-2} operator:
    Single integral of generic d2m
    """
    get_diags(D2m_generic, 1)


def i4D2p():
    """
    I4 of D^{2} operator
    Double integral of generic d2p
    """
    get_diags(D2p_generic, 2)

def i4D2m():
    """
    I4 of D^{-2} operator
    Double integral of generic d2m
    """
    get_diags(D2m_generic, 2)


def i3D3p():
    """
    I3 of D^{3} operator
    Just generic d3p
    """
    get_diags(D3p_generic, 0)

def i3D3m():
    """
    I3 of D^{-3} operator
    Just generic d3m
    """
    get_diags(D3m_generic, 0)


def i4D3p():
    """
    I4 of D^{3} operator
    Single integral of generic d3p
    """
    get_diags(D3p_generic, 1)

def i4D3m():
    """
    I4 of D^{-3} operator
    Single integral of generic d3m
    """
    get_diags(D3m_generic, 1)


def i4D4p():
    """
    I4 of D^{4} operator
    Just generic d4p
    """
    get_diags(D4p_generic, 0)

def i4D4m():
    """
    I4 of D^{-4} operator
    Just generic d4m
    """
    get_diags(D4m_generic, 0)


def i3laplDp():
    """
    I3 of lapl D^{1} operator
    Just generic laplDp
    """
    get_diags(laplDp_generic, 0)

def i3laplDm():
    """
    I3 of lapl D^{-1} operator
    Just generic laplDm
    """
    get_diags(laplDm_generic, 0)
    return 


def i4laplDp():
    """
    I4 of lapl D^{1} operator
    Single integral of generic laplDp
    """
    get_diags(laplDp_generic, 1)

def i4laplDm():
    """
    I4 of lapl D^{-1} operator
    Single integral of generic laplDm
    """
    get_diags(laplDm_generic, 1)
    return 


def i4laplD2p():
    """
    I4 of lapl D^{2} operator
    Just generic laplD2p
    """
    get_diags(laplD2p_generic, 0)

def i4laplD2m():
    """
    I4 of lapl D^{-2} operator
    Just generic laplD2m
    """
    get_diags(laplD2m_generic, 0)
    return 





# def I_n__Y_q_m_n(expr_in, 
#                  n=0, q=0, 
#                  asrow_final=False):
#     """
#     I^n Y^(q-n) acting on J^(l+q) -----> J^l 
#     i.e. the n'th term of the q'th order sum
#     """

#     l_offsets = reversed(range(q))

#     expr = expr_in
#     for s in l_offsets:
#         l_level = get_offset_level(s)

#         finish = (s==0)
#         ar = asrow_final if finish else False 
        
#         if s > n-1:
#             expr = l_level.spectral_decrease(expr, asrow=ar)
#         else:
#             expr = l_level.spectral_integral_decrease(expr, asrow=ar)

#     return expr