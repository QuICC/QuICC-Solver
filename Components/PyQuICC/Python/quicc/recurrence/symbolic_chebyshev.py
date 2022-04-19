"""Module used for symbolic manipulations of sparse chebyshev expansions"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import copy

import quicc.recurrence.symbolic_base as base

class SymbolicChebyshev(base.SymbolicBase):
    """Class to compute the symbolic recurrence relation for a Chebyshev expansion"""

    def spectral_monomial(self, p, f, asrow = True):
        """Recurrence relation for the multiplication by a monomial of a Chebyshev expansion"""

        # Some assertion for safety
        assert (p >= 0)

        # End recurrence: monomial order zero is identity
        if p == 0:
            recurrence = dict(f)
        else:
            prev = self.spectral_monomial(p-1, f)
            recurrence = dict()

            for i in prev.keys():
                recurrence[i-1] = recurrence.get(i-1,0) + sympy.Rational(1,2)*prev[i]
                recurrence[i+1] = recurrence.get(i+1,0) + sympy.Rational(1,2)*prev[i]
            for i in recurrence.keys():
                recurrence[i] = recurrence[i].simplify().factor()

        return recurrence

    def spectral_integral(self, q, f, asrow = True):
        """Recurrence relation for the indefinite integral of a Chebyshev expansion"""

        # Some assertion for safety
        assert (q >= 0)

        # End recurrence: integration "power" zero is identity
        if q == 0:
            recurrence = dict(f)
        else:
            prev = self.spectral_integral(q-1, f, False)
            recurrence = dict()

            n = sympy.Symbol('n')
            for i in prev.keys():
                recurrence[i-1] = recurrence.get(i-1,0) - (1/(2*(n+i-1)))*prev[i]
                recurrence[i+1] = recurrence.get(i+1,0) + (1/(2*(n+i+1)))*prev[i]

            # Convert to row recurrence relation
            if asrow:
                old = dict(recurrence)
                recurrence = dict()
                for i in old.keys():
                    recurrence[-i] = old[i].subs(n, n-i)

            for i in recurrence.keys():
                recurrence[i] = recurrence[i].simplify().factor()

        return recurrence

    def change_variable(self, terms, type):
        """Compute change of variable from r = [0, 1] to x = [-1,1]"""

        new_terms = []
        if type in ['linear_r2x', 'linear']:
            a, b, x = sympy.symbols('a b x')
            for term in terms:
                poly = sympy.poly((a*x+b)**term['p'], x)
                pcoeffs = list(reversed(poly.coeffs()))
                for j,c in enumerate(pcoeffs):
                    t = term.copy()
                    t['c'] = c*t['c']*a**(t['q'] - t['d'])
                    t['p'] = j
                    new_terms.append(t)

        return new_terms
