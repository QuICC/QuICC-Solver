"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
from mpi4py import MPI

def write_header(f, name, nRes, nEigs, eq_params):
    """Write marginal curve file header"""

    if MPI.COMM_WORLD.Get_rank() == 0:

        # file is not empty add double newline (gnuplot reads it separately)
        if f.tell() != 0:
            f.write('\n\n')

        # First header
        f.write('#Results for: ' + name + '\n')

        header = []

        # Resolution
        for i in range(0, nRes):
            header.append('Res_' + str(i))

        # Wavenumber
        for i in range(0, nEigs):
            header.append('k_' + str(i))

        # Parameters
        for k,v in sorted(eq_params.items()):
            if k == 'rayleigh':
                continue
            header.append(k)

        # Critical Rayleigh number
        header.append('Rac')

        # Critical frequencies
        header.append('freq')

        # Mode index
        header.append('mode')

        # Convergence (growth rate)
        header.append('convergence')

        f.write('#'+'\t'.join(header) + '\n')

def write_results(f, res, eigs, eq_params, Rac, freq, mode, conv):
    """Write marginal curve point to file"""

    if MPI.COMM_WORLD.Get_rank() == 0:
        result = []

        # Resolution
        for r in res:
            result.append("{:d}".format(r))

        # Wavenumber
        for k in eigs:
            result.append("{:.14g}".format(k))

        # Parameters
        for k,v in sorted(eq_params.items()):
            if k == 'rayleigh':
                continue
            result.append("{:.14g}".format(v))

        # Critical rayleigh number
        result.append("{:.14g}".format(Rac))

        # Critical frequency
        result.append("{:.14g}".format(freq))

        # Mode index
        result.append("{:d}".format(mode))

        # Convergence (growth rate)
        result.append("{:.2g}".format(conv))

        f.write('\t'.join(result) + '\n')
