"""Module used for symbolic manipulations for general integration by parts method"""

from __future__ import division
from __future__ import unicode_literals

import sympy
import copy

class SymbolicBase:
    """Base class to compute the symbolic recurrence relation"""

    # Use SymPy simplify before factor
    useSimplify = False

    def spectral_monomial(self, p, f, asrow = True):
        """This needs to be defined in polynomial specific implementation"""

        raise NotImplementedError("Monomial operation not defined")

    def spectral_integral(self, p, f, asrow = True):
        """This needs to be defined in polynomial specific implementation"""

        raise NotImplementedError("Integral operation not defined")

    def integrate(self, term, base = None):
        """Compute a symbolic integration by parts"""

        # Some assertion for safety
        assert (term['q'] >= 0)
        assert (term['p'] >= 0)
        assert (term['d'] >= 0)
        assert (term['q'] >= term['d'])

        if base is None: base = dict()

        # End recurrence: derivative order is zero
        if term['d'] == 0:
            base[term['q'],term['p']] = base.get((term['q'],term['p']),0) + term['c']
        else:
            # Compute first part of integration by parts
            t = {'q':term['q'] - 1, 'p':term['p'], 'd':term['d'] - 1, 'c':term['c']}
            base = self.integrate(t, base)

            # Compute second part of integration by parts
            if term['p'] - 1 >= 0:
                t = {'q':term['q'], 'p':term['p'] - 1, 'd':term['d'] - 1, 'c':-term['p']*term['c']}
                base = self.integrate(t, base)

        return base

    def spectral_intgxmult(self, q, p, f, finish = True):
        """Symbolic integration and multiplication by monomial"""

        # Some assertion for safety
        assert (q >= 0)
        assert (p >= 0)

        # Compute monomial multiplication
        asrow = (p == 0 or q == 0)
        recurrence = self.spectral_monomial(p, f, asrow)

        # Compute integration
        recurrence = self.spectral_integral(q, recurrence, asrow = finish)

        return recurrence

    def build_recurrence(self, terms, fs, finish = True):
        """Recursively build the full recurrence relation"""

        ops = dict();

        # Integrate and sum all terms
        for term in terms:
            op = self.integrate(term)
            for tup in op.keys():
                ops[tup] = ops.get(tup, 0) + op[tup]

        # Create recurrence relation
        recurrence = dict()
        for tup in ops.keys():
            for d,f in fs.items():
                rec = self.spectral_intgxmult(tup[0], tup[1], {d:f}, finish)
                for i in rec.keys():
                    recurrence[i] = recurrence.get(i,0) + ops[tup]*rec[i]

        for i in recurrence.keys():
            if self.useSimplify:
                recurrence[i] = recurrence[i].simplify().factor()
            else:
                recurrence[i] = recurrence[i].factor()

        return recurrence
