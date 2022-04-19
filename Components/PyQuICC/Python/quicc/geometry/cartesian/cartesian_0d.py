"""Module provides functions to generate sparse operators in a triply periodic cartesian box."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp


def zblk():
    """Create a block of zeros"""

    mat = spsp.lil_matrix((1,1))
    return mat

def qid():
    """Create an identity block"""

    return spsp.identity(1)

def lapl(k, l, m):
    """Create operator for triply periodic Laplacian"""

    return -(k**2 + l**2 + m**2)*spsp.identity(1)

def lapl2(k, l, m):
    """Create operator for triply periodic Laplacian^2"""

    return ((k**2 + l**2 + m**2)**2)*spsp.identity(1)
