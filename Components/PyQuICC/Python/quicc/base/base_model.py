"""Module provides the functions required for any model"""

from __future__ import division
from __future__ import unicode_literals

verbose_write_mtx = False
if verbose_write_mtx:
    import scipy.io as io
    import os

    def make_single_name(fields, bcs):
        pid = str(os.getpid())
        return str(bcs["bcType"]) + "_" + fields[0][0] + '_' + fields[0][1] + '_' + pid

    def make_double_name(field_row, field_col, bcs):
        pid = str(os.getpid())
        return ''.join(field_row) + '_' + ''.join(field_col) + '_' + str(bcs["bcType"]) + "_" + pid

import scipy.sparse as spsp
import quicc.base.utils as utils


class BaseModel:
    """Base class for all the models"""

    linearize = False

    # Same as C++ ModelOperatorBoundary
    SOLVER_HAS_BC = "solver_has_bc"
    SOLVER_NO_TAU = "solver_no_tau"
    STENCIL = "stencil"
    FIELD_TO_RHS = "field_to_rhs"

    # Same as C++ CouplingIndexType
    SLOWEST_SINGLE_RHS = 0
    SLOWEST_MULTI_RHS = 1
    MODE = 2
    SINGLE = 3

    # Internal to Python
    EXPLICIT_LINEAR = 0
    EXPLICIT_NONLINEAR = 1
    EXPLICIT_NEXTSTEP = 2

    def time(self, res, eq_params, eigs, bcs, fields, restriction = None):
        """Create the time derivative operator"""

        mat = utils.build_diag_matrix(fields, self.time_block, (res,eq_params,eigs,bcs), restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_time_" + make_single_name(fields, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def implicit_linear(self, res, eq_params, eigs, bcs, fields, restriction = None):
        """Create the implicit linear operator"""

        mat = utils.build_block_matrix(fields, self.implicit_block, (res,eq_params,eigs,bcs), restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_linear_" + make_single_name(fields, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname  + ".mtx", mat)
        return mat

    def boundary(self, res, eq_params, eigs, bcs, fields, restriction = None):
        """Create the boundary operator"""

        mat = utils.build_block_matrix(fields, self.boundary_block, (res,eq_params,eigs,bcs), restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_boundary_" + make_single_name(fields, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname  + ".mtx", mat)
        return mat

    def explicit_linear(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.explicit_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_linear_" + make_double_name(field_row, field_col, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def explicit_nonlinear(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.nonlinear_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_nonlinear_" + make_double_name(field_row, field_col, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def explicit_nextstep(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.nextstep_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_nextstep_" + make_double_name(field_row, field_col, bcs)
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def compile_equation_info(self, res, field_row, is_complex, index_mode):
        """Collect all equation info together"""

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        lin_fields = self.explicit_fields(self.EXPLICIT_LINEAR, field_row)
        # Additional explicit nonlinear fields
        nl_fields = self.explicit_fields(self.EXPLICIT_NONLINEAR, field_row)
        # Additional explicit update for next step linear fields
        next_fields = self.explicit_fields(self.EXPLICIT_NEXTSTEP, field_row)

        return (is_complex, im_fields, lin_fields, nl_fields, next_fields, index_mode)

    def operator_info(self, res, eigs, bcs, field_row):
        """Combine block sizes into global operator info"""

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)

        # Compute block info
        info = self.block_size(res, eigs, bcs, field_row)

        # Compute system size
        sys_n = 0
        for f in im_fields:
            sys_n += self.block_size(res, eigs, bcs, f)[1]

        if sys_n == 0:
            sys_n = info[1]
        info = info + (sys_n,)

        return info

    def stability_sizes(self, res, eigs, bcs):
        """Get the block sizes in the stability calculation matrix"""

        # Block sizes
        blocks = []
        for f in self.stability_fields():
            blocks.append(self.block_size(res, eigs, bcs, f)[1])

        # Invariant size (local dimension in spectral space, no restriction)
        invariant = (res[0],)*len(self.stability_fields())

        # Index shift
        shift = 0

        return (blocks, invariant, shift)

    def automatic_parameters(self, eq_params):
        """Extend parameters with automatically computable values"""

        return dict()

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        raise NotImplementedError("Model cannot be used for linear stability calculations!")

    def block_size(self, res, eigs, bcs, field_row):
        """Create block size information"""

        raise NotImplementedError("Operator block sizes have not been defined!")

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        raise NotImplementedError("Model cannot be used for linear stability calculations!")

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        raise NotImplementedError("Model should implement this method!")

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        raise NotImplementedError("Model should implement this method!")

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for boundary operator"""

        raise NotImplementedError("Model should implement this method!")

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for implicit operator"""

        raise NotImplementedError("Model should implement this method!")

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        raise NotImplementedError("Model should implement this method!")

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nonlinear term"""

        raise NotImplementedError("Model should implement this method!")

    def nextstep_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nextstep update"""

        raise NotImplementedError("Model should implement this method!")

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""

        raise NotImplementedError("Stencil needs to be implemented in model!")
