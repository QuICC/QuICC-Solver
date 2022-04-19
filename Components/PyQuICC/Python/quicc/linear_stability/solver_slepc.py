"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import random
import numpy as np
import numpy.linalg as nplin
import scipy.linalg as splin
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import sys, slepc4py
slepc4py.init(sys.argv)

from mpi4py import MPI

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

class GEVPSolver:
    """GEVP Solver using on SLEPc"""

    def __init__(self, shift_range = None, tol = 1e-8, ellipse_radius = None, fixed_shift = False, target = None, euler = None, conv_idx = 1, spectrum_conv = 1, geometry = None, impose_symmetry = False, use_spherical_evp = False, whichTarget = None):
        """Initialize the SLEPc solver"""

        self.tol = tol
        self.ellipse_radius = ellipse_radius
        self.fixed_shift = fixed_shift
        self.target = target
        self.euler = euler
        self.conv_idx = conv_idx
        self.spectrum_conv = spectrum_conv
        self.geometry = geometry
        self.impose_symmetry = impose_symmetry
        self.use_spherical_evp = use_spherical_evp
        self.whichTarget = whichTarget

        if shift_range is None:
            #self.shift_range = (1e-2, 0.2)
            #self.shift_range = (-1e-1, 1e-1)
            self.shift_range = (-1e-2, 1e-2)
        else:
            self.shift_range = shift_range

        if self.fixed_shift:
            self.setRandomShift()

        self.create_eps()

    def create_eps(self):
        """Create SLEPc's eigensolver"""

        opts = PETSc.Options()
        opts["mat_mumps_icntl_14"] = 200
        opts["mat_mumps_icntl_29"] = 2
        opts["mat_mumps_cntl_1"] = 0.1

        if self.ellipse_radius is not None:
            opts['rg_type'] = 'ellipse'
            opts['rg_ellipse_center'] = 0
            opts['rg_ellipse_radius'] = self.ellipse_radius
            opts['rg_ellipse_vscale'] = 1.0

        self.E = SLEPc.EPS()
        self.E.create()

        self.E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
        if self.target is not None:
            self.E.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
            #self.E.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_IMAGINARY)
        else:
            self.E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
            #self.E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)

        # Override target type
        if self.whichTarget is not None:
            self.E.setWhichEigenpairs(self.whichTarget)

        self.E.setBalance(SLEPc.EPS.Balance.TWOSIDE)
        self.E.setTolerances(tol = self.tol)

        ST = self.E.getST()
        ST.setType('sinvert')
        if not self.fixed_shift:
            self.setRandomShift()
        if self.target is not None:
            self.E.setTarget(self.target)
            self.shift = self.target
        ST.setShift(self.shift)

        KSP = ST.getKSP() # Krylov subspace method
        KSP.setType('preonly')

        PC = KSP.getPC() # Preconditioning
        PC.setType('lu')
        PC.setFactorSolverType('mumps')

        PC.setFromOptions()
        KSP.setFromOptions()
        ST.setFromOptions()

        self.E.setFromOptions()

    def setRandomShift(self):
        """Compute random spectra transform shift and store it"""

        if self.shift_range[0] == self.shift_range[1]:
            self.shift = self.shift_range[0]
        else:
            rnd = PETSc.Random()
            rnd.create(comm = MPI.COMM_SELF)
            rnd.setType(PETSc.Random.Type.RAND)
            rnd.setInterval(self.shift_range)
            self.shift = rnd.getValueReal()

    def set_initial_vector(self, v, sizes):
        """Create an initial vector with converged spectrum"""

        v.set(1.0)
        data = v.getArray()
        rstart, rend = v.getOwnershipRange()
        start = 0
        if self.geometry in ['sphere_worland', 'shell'] and self.impose_symmetry:
            par = [1, 0, 0]
            gal = [2, 4, 2]
            for i, sze in enumerate(sizes[0]):
                for j in range(0, sze):
                    if start + j >= rstart and start + j < rend:
                        if (j//(sizes[1][0]-gal[i]))%2 == par[i]:
                            data[start+j-rstart] = ((2*np.random.ranf() - 1.0)+(2*np.random.ranf() - 1.0)*1j)*10**(-10.0*(j)/float(sze))*10**(-10*(j%sizes[1][0])/float(sizes[1][0]))
                        else:
                            data[start+j-rstart] = 0.0
                    elif start + j >= rend:
                        break
                start = start + sze
        else:
            for i, sze in enumerate(sizes[0]):
                for j in range(0, sze):
                    if start + j >= rstart and start + j < rend:
                        data[start+j-rstart] = ((2*np.random.ranf() - 1.0)+(2*np.random.ranf() - 1.0)*1j)*10**(-10.0*(j)/float(sze))*10**(-10*(j%sizes[1][0])/float(sizes[1][0]))
                    elif start + j >= rend:
                        break
                start = start + sze
        Print("    Generated spectrum")
        v.createWithArray(data)

    def update_eps(self, A, B, nev, sizes, initial_vector = None):
        """Create SLEPc eigensolver"""

        self.create_eps()

        self.E.setOperators(A,B)
        if initial_vector is not None:
            self.E.setInitialSpace(initial_vector)
        elif self.euler is not None:
            Print('Initialise vector with Euler steps')
            # Create initial vector through timestepping
            euler_nstep = self.euler[0]
            euler_step = self.euler[1]
            vrnd, v = A.getVecs()
            self.set_initial_vector(vrnd, sizes)
            vrnd = -(1.0/euler_step)*B*vrnd
            ksp = PETSc.KSP().create()
            ksp.setType('preonly')
            pc = ksp.getPC()
            pc.setType('lu')
            pc.setFactorSolverPackage('mumps')
            AA = A - (1.0/euler_step)*B
            ksp.setOperators(AA)
            ksp.setFromOptions()
            for i in range(0,euler_nstep):
                ksp.solve(vrnd, v)
                vrnd = (-1.0/euler_step)*B*v
                #vrnd.normalize()
            self.E.setInitialSpace(v)
            Print('Done')

        self.E.setDimensions(nev = nev)

    def eigenvalues(self, system, nev, initial_vector = None):
        """Compute eigenvalues using SLEPc"""

        if self.geometry in ['shell', 'sphere_worland'] and self.use_spherical_evp:
            return self.eigenvalues_spherical(system, nev, initial_vector = initial_vector)
        else:
            return self.eigenvalues_simple(system, nev, initial_vector = initial_vector)

    def eigenpairs(self, system, nev, initial_vector = None):
        """Compute eigenpairs using SLEPc"""

        if self.geometry in ['shell', 'sphere_worland'] and self.use_spherical_evp:
            return self.eigenpairs_spherical(system, nev, initial_vector = initial_vector)
        else:
            return self.eigenpairs_simple(system, nev, initial_vector = initial_vector)

    def eigenvalues_simple(self, system, nev, initial_vector = None):
        """Compute eigenvalues using SLEPc"""

        pA, pB = self.petsc_operators(*system)

        self.update_eps(pA, pB, nev, system[-1], initial_vector = initial_vector)

        self.E.solve()
        nconv = self.E.getConverged()

        if nconv >= nev:
            eigs = np.array(np.zeros(nev), dtype=complex)
            for i in range(nev):
                eigs[i] = self.E.getEigenvalue(i)

            return eigs
        else:
            return None

    def eigenpairs_simple(self, system, nev, initial_vector = None):
        """Compute eigenpairs using SLEPc"""

        pA, pB = self.petsc_operators(*system)

        self.update_eps(pA, pB, nev, system[-1], initial_vector = initial_vector)

        self.E.solve()
        nconv = self.E.getConverged()

        if nconv >= nev:
            # Create the results vectors
            vr, wr = pA.getVecs()
            vi, wi = pA.getVecs()
            eigs = np.array(np.zeros(nev), dtype=complex)
            vects = []
            for i in range(nev):
                eigs[i] = self.E.getEigenpair(i, vr, vi)
                vects.append(vr + 1j*vi)

            return (eigs, vects)
        else:
            return (None, None)

    def eigenvalues_spherical(self, system, nev, initial_vector = None):
        """Compute eigenvalues using SLEPc"""

        pair = self.eigenpairs(system, nev, initial_vector = initial_vector)

        return pair[0]

    def eigenpairs_spherical(self, system, nev, initial_vector = None):
        """Compute eigenpairs using SLEPc"""

        pA, pB = self.petsc_operators(*system)

        self.update_eps(pA, pB, nev, system[-1], initial_vector = initial_vector)

        comm = MPI.COMM_WORLD
        c_nev = 1
        pair = (None, None)
        tmp_eigs = []
        def_space = []
        bad_vects = 0
        vects = []
        for attempts in range(0,10):
            for i in range(0, len(def_space)):
                self.E.setDeflationSpace(def_space[i])
            if initial_vector is not None:
                self.E.setInitialSpace(initial_vector)
            elif attempts > 0:
                vtmpL, vtmpR = pA.getVecs()
                self.set_initial_vector(vtmpL, system[-1])
                self.E.setInitialSpace(vtmpL)
            self.E.setDimensions(nev = c_nev)

            self.E.solve()
            nconv = self.E.getConverged()

            if nconv >= c_nev:
                # Create the results vectors
                vr, wr = pA.getVecs()
                vi, wi = pA.getVecs()
                rstart, rend = vr.getOwnershipRange()
                for i in range(nconv):
                    lmb = self.E.getEigenpair(i, vr, vi)
                    if lmb.imag < 0:
                        convratio = 1e99
                        start = 0
                        for sze in system[-1][0]:
                            spec_info = np.zeros(2)
                            tmp_info = np.zeros(2)
                            rank_start = max(rstart, start)
                            rank_end = min(rend,start+sze)
                            # Get max of spectrum
                            if rank_start < rank_end:
                                tmp_info[0] = np.max(np.abs(vr.getArray()[rank_start-rstart:rank_end-rstart]))
                            # Get min of spectrum
                            rank_start = max(rstart, start+sze-self.conv_idx)
                            if rank_start < rank_end:
                                tmp_info[1] = np.max(np.abs(vr.getArray()[rank_start-rstart:rank_end-rstart]))
                            comm.Allreduce(tmp_info, spec_info, MPI.MAX)
                            convratio = min(convratio, spec_info[0]/spec_info[1])
                            start = start + sze
                        Print('Spectral convergence test for ' + str(lmb) + ': ' + str(convratio))
                        if convratio > self.spectrum_conv:
                            tmp_eigs.append(lmb)
                            vects.append(vr + 1j*vi)
                        else:
                            def_space.append(vr + 1j*vi)
                            bad_vects += 1
                    else:
                        def_space.append(vr + 1j*vi)
                        bad_vects += 1
                Print('Converged eigenvalues: ' + str(tmp_eigs))

                if len(tmp_eigs) >= nev:
                    pair = (np.array(tmp_eigs[0:nev], dtype=complex), vects)
                    break
                else:
                    c_nev = max(nev, 1 + len(vects))
                    bad_vects = 0
                    if len(vects) > 0:
                        initial_vector = vects[0]
                        tmp_eigs = []
                        vects = []
                        #def_space.append(vects)

        return pair

    def restrict_operators(self, sizes):
        """Compute restriction for operators"""

        if MPI.COMM_WORLD.Get_size() > 1:
            pTmp = PETSc.Vec().create()
            pTmp.setSizes(np.sum(sizes[0]))
            pTmp.setUp()
            rstart, rend = pTmp.getOwnershipRange()

            restrict = []
            bstart = 0
            tot = 0
            for s, l in zip(sizes[0], sizes[1]):
                if rstart < tot + s:
                    bstart = max(np.floor((rstart-tot)/l),0) + sizes[2]
                    bend = min(np.ceil((rend-tot)/l),np.ceil(s/l)) + sizes[2]
                    restrict.append(np.arange(bstart,bend))
                else:
                    restrict.append(np.array([]))
                tot = tot + s

        else:
            restrict = None

        return restrict

    def petsc_operators(self, opA, opB, opC, sizes):
        """Convert SciPy operators to PETSc operators"""

        # Build operator restriction
        #restrict = self.restrict_operators(sizes)
        restrict = None

        # Setup A matrix
        A = (opA(restriction = restrict) + opC(restriction = restrict)).transpose().tocsr()
        pA = PETSc.Mat().create()
        pA.setSizes(A.shape)
        pA.setUp()

        # Fill A matrix
        rstart, rend = pA.getOwnershipRange()
        pA.createAIJ(size=A.shape, nnz=A.getnnz(1)[rstart:rend+1], csr=(A.indptr[rstart:rend+1] - A.indptr[rstart], A.indices[A.indptr[rstart]:A.indptr[rend]], A.data[A.indptr[rstart]:A.indptr[rend]]))
        A = None
        pA.assemblyBegin()

        # Setup B matrix
        B = opB(restriction = restrict).transpose().tocsr()
        pA.assemblyEnd()
        pA.transpose()

        pB = PETSc.Mat().create()
        pB.setSizes(B.shape)
        pB.setUp()

        # Fill B matrix
        rstart, rend = pB.getOwnershipRange()
        pB.createAIJ(size=B.shape, nnz=B.getnnz(1)[rstart:rend+1], csr=(B.indptr[rstart:rend+1] - B.indptr[rstart], B.indices[B.indptr[rstart]:B.indptr[rend]], B.data[B.indptr[rstart]:B.indptr[rend]]))
        pB.assemble()
        B = None
        pB.transpose()

        return (pA, pB)

class NewtonNoneError(Exception):
    pass

class NewtonDoneError(Exception):
    pass

def newton(func, x0, args=(), tol=1.48e-8, maxiter=50, step = 1e-4):
    """
    Find a zero using secant method.
    Find a zero of the function `func` given a nearby starting point `x0`.
    Parameters
    ----------
    func : function
        The function whose zero is wanted. It must be a function of a
        single variable of the form f(x,a,b,c...), where a,b,c... are extra
        arguments that can be passed in the `args` parameter.
    x0 : float
        An initial estimate of the zero that should be somewhere near the
        actual zero.
    args : tuple, optional
        Extra arguments to be used in the function call.
    tol : float, optional
        The allowable error of the zero value.
    maxiter : int, optional
        Maximum number of iterations.
    Returns
    -------
    zero : float
        Estimated location where function is zero.
    """
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    # Secant method
    p0 = x0
    q0 = func(*((p0,) + args))
    if q0 >= 0:
        p1 = x0*(1.0 - step) + step
    else:
        p1 = x0*(1.0 + step) + step
    q1 = func(*((p1,) + args))
    for iter in range(maxiter):
        if q1 == q0:
            if p1 != p0:
                msg = "Tolerance of %s reached" % (p1 - p0)
                warnings.warn(msg, RuntimeWarning)
            return (p1 + p0)/2.0
        else:
            p = p1 - q1*(p1 - p0)/(q1 - q0)
            if p < 0:
                factor = 15.0/10.0
                if q1 < 0:
                    p = p1*factor
                else:
                    p = p1/factor
        if abs(p - p1) < tol:
            return p
        p0 = p1
        q0 = q1
        p1 = p

        try:
            q1 = func(*((p1,) + args))
        except NewtonNoneError:
            if abs(p1/p0) > 10:
                p1 = 10*p0
                q1 = func(*((p1,) + args))
            elif abs(p1/p0) > 2:
                p1 = (p0 + p1)/2.0
                q1 = func(*((p1,) + args))
            else:
                raise
        except:
            raise

    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)
