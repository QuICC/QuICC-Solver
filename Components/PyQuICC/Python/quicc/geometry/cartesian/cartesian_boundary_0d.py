"""Module provides functions to generate the boundary conditions in a cartesian 0D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import quicc.base.utils as utils


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}
