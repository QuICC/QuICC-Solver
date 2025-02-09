"""Module provides the functions to generate the Boussinesq inviscid inertial wave in a sphere with Worland expansion (Toroidal/Poloidal formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.sphere_worland as geo
import quicc.base.base_model as base_model
from quicc.geometry.spherical.sphere_boundary_worland import no_bc


class BoussinesqIIWSphere(base_model.BaseModel):
    """Class to setup the Boussinesq inviscid inertial wave in a sphere with Worland expansion (Toroidal/Poloidal formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return []

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","tor"), ("velocity","pol")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocity","tor"), ("velocity","pol")]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit linear dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, eigs, bcs, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row == ("velocity","tor"):
                shift_r = 0
            elif field_row == ("velocity","pol"):
                shift_r = 1
            else:
                shift_r = 0

            gal_n = (res[0] - shift_r)

        else:
            gal_n = tau_n
            shift_r = 0

        block_info = (tau_n, gal_n, (shift_r,0,0), 1)
        return block_info

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return geo.stencil(res[0], res[1], m, bc, make_square)

    def stability_sizes(self, res, eigs, bcs):
        """Get the block sizes in the stability calculation matrix"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        # Block sizes
        blocks = []
        for f in self.stability_fields():
            blocks.append(self.block_size(res, eigs, bcs, f)[1]*(res[1]-m))

        # Invariant size (local dimension in spectral space, no restriction)
        invariant = (res[0],)*len(self.stability_fields())

        # Index shift
        shift = m

        return (blocks, invariant, shift)

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_SINGLE_RHS

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:0, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-10, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:0}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:10}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 0
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 1

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {0:0, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-10, 'rt':1}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 0
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 1

        else:
            bc = no_bc()

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = geo.i1(res[0], res[1], m, bc, 1j*m, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i1coriolis(res[0], res[1], m, bc, -1.0, l_zero_fix = 'zero', restriction = restriction)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = geo.i2coriolis(res[0], res[1], m, bc, 1.0, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i2lapl(res[0], res[1], m, bc, 1j*m, l_zero_fix = 'zero', restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = geo.i1(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("velocity","pol"):
            mat = geo.i2lapl(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'set', restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block of boundary operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row in [("velocity","tor"), ("velocity","pol")] and field_row == field_col:
            mat = geo.zblk(res[0], res[1], m, bc, l_zero_fix = 'zero', restriction = restriction)
        else:
            mat = geo.zblk(res[0], res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat
