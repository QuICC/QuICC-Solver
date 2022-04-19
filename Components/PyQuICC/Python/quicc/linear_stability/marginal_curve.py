"""Module provides the functor required to trace a marginal curve"""

from __future__ import division
from __future__ import unicode_literals

import copy
import numpy as np
import scipy.optimize as optimize
import functools

import h5py
import tables

import quicc.linear_stability.solver_slepc as solver_mod
import quicc.linear_stability.viewer as viewer
import quicc.linear_stability.io as io
from quicc.linear_stability.solver_slepc import Print
from mpi4py import MPI
from petsc4py import PETSc

PETSc.ScalarType = np.complex64

class MarginalCurve:
    """The marginal curve"""

    def __init__(self, gevp_opts, mode, rtol = 1e-7, evp_tol = 1e-8):
        """Initialize the marginal curve"""

        # Open file for IO and write header
        if MPI.COMM_WORLD.Get_rank() == 0:
            self.out = open('marginal_curve.dat', 'a', 1)
        else:
            self.out = None
        io.write_header(self.out, gevp_opts['model'].__class__.__name__, len(gevp_opts['res']), len(gevp_opts['eigs']), gevp_opts['eq_params'])

        self.point = MarginalPoint(gevp_opts, mode, rtol = rtol, out_file = self.out, evp_tol = evp_tol)
        self.rtol = rtol
        self.guess = None

    def trace(self, ks, guess = None, initial_guess = 1.0):
        """Trace the marginal curve"""

        Print("Tracing marginal curve")
        Print("----------------------")

        data_k = np.zeros(len(ks))
        data_Ra = np.zeros(len(ks))
        data_freq = np.zeros(len(ks))
        Rac = guess
        for i, k in enumerate(ks):
            if i > 3:
                Rac = p(k)
            Rac, freq = self.point(k, guess = Rac, initial_guess = initial_guess)
            data_k[i] = k
            data_Ra[i] = Rac
            data_freq[i] = freq
            if i > 2:
                lb = max(0, i-3)
                z = np.polyfit(data_k[lb:i+1], data_Ra[lb:i+1], min(i,3))
                p = np.poly1d(z)

        return (data_k, data_Ra, data_freq)

    def minimum(self, ks, Ras, xtol = 1e-8, only_int = False):
        """Compute the critical value (minimum of curve)"""

        # Open file for IO and write header
        if MPI.COMM_WORLD.Get_rank() == 0:
            min_out = open('marginal_minimum.dat', 'a', 1)
        else:
            min_out = None
        io.write_header(min_out, self.point._gevp.model.__class__.__name__, len(self.point._gevp.res), len(self.point._gevp.eigs), self.point._gevp.eq_params)

        Print("Finding minimum")
        Print("---------------")

        # Get bounds for minimum
        imin = np.argmin(Ras)
        a = ks[imin-1]
        b = ks[imin]
        c = ks[imin+1]
        self.guess = Ras[imin]
        if imin == 0 or imin == len(Ras)-1:
            Print("Minimum is not in provided range")
            return (None, None, None)
        else:
            if only_int:
                self.point(b, guess = Ras[imin])
                min_kc = b
                min_Rac = Ras[imin]
                min_fc = self.point.frequency()

            else:
                # Find minimum
                opt = optimize.minimize_scalar(self._tracker, [a, b, c], options = {'xtol':xtol})
                min_kc = opt.x
                min_Rac = opt.fun
                min_fc = self.point.frequency()

            Print("Minimumt at kc: {:g}, Rac: {:g}, fc: {:g}".format(min_kc, min_Rac, min_fc))
            io.write_results(min_out, self.point._gevp.res, self.point._gevp.wave(min_kc), self.point._gevp.eq_params, min_Rac, min_fc, self.point.mode, self.point.growthRate())
            return (min_kc, min_Rac, min_fc)

    def view(self, ks, Ras, fs, minimum = None, plot = True, save_pdf = False):
        """Plot the marginal curve"""

        if plot or save_pdf:
            import matplotlib.pylab as pl

            pl.subplot(1,2,1)
            pl.plot(ks, Ras, 'b')
            if minimum is not None:
                pl.plot(minimum[0], minimum[1], 'r+', markersize=14)
            pl.title('Critical Rayleigh number')
            pl.xlabel('Wavenumber')
            pl.ylabel('Ra')
            pl.tight_layout()

            pl.subplot(1,2,2)
            pl.plot(ks, fs, 'g')
            pl.title('Critical frequency')
            pl.xlabel('Wavenumber')
            pl.ylabel('frequency')
            if save_pdf:
                pl.savefig('marginal_curve.pdf', bbox_inches='tight', dpi = 200)
            if plot:
                pl.tight_layout()
                pl.show()
                pl.clf()

    def _tracker(self, k):

        return self.point(k, guess = self.guess)[0]

class MarginalPoint:
    """A point on the marginal curve"""

    def __init__(self, gevp_opts, mode, rtol = 1e-7, out_file = None, evp_tol = 1e-8):
        """Initialize the marginal point"""

        self.out = out_file
        gevp_opts['tol'] = evp_tol
        self._gevp = GEVP(**gevp_opts)
        self.mode = mode
        self.Rac = None
        self.fc = None
        self.a = None
        self.b = None
        self.rtol = rtol
        self.nev = max(mode+1, 1)
        self.best = None

    def __call__(self, k, guess = None, initial_guess = 1.0):
        """Compute marginal point in given interval [a, b]"""

        if guess is None:
            Ra = initial_guess
            factor = 2.0
        else:
            Ra = guess
            factor = 1.1

        Print("Computing marginal curve at k = {:g}".format(k))
        Print("\t- initial guess Ra = {:g}".format(Ra))

        # Set wavenumter
        self._gevp.setEigs(k)

        # Search interval (newton)
        self._gevp.best = None
        Print("\t- finding interval")
        for i in range(0,16):
            self.findInterval(Ra)
            if self.a != None and self.b != None:
                break
            elif self.best is not None and self.best[2]:
                Ra = self.best[1]
                self.best = (self.best[0], self.best[1], False)
            elif self.best is not None:
                if self.best[0] > 0:
                    Ra = self.best[1]/(1.1*factor)
                else:
                    Ra = self.best[1]*(1.1*factor)
                self.best = None
            else:
                Ra = factor*Ra

        if self.a is None or self.b is None:
            self.Rac = float('NaN')
            self.rc = float('NaN')
            self.a = None
            self.b = None
        else:
            # Refine zero with Brent's method
            Print("\t- Brent's method: [{:g}, {:g}]".format(self.a,self.b))
            Rac = optimize.brentq(self._tracker, self.a, self.b, rtol = self.rtol)
            self.a = None
            self.b = None
            self.Rac = Rac
            self.fc = self.frequency()

        # Write results to file and show on console output
        io.write_results(self.out, self._gevp.res, self._gevp.eigs, self._gevp.eq_params, self.Rac, self.fc, self.mode, self.growthRate())
        Print("\t- Ra = {:g}, freq = {:g}, conv = {:g}".format(self.Rac, self.fc, self.growthRate()))

        return (self.Rac, self.fc)

    def findInterval(self, guess):
        """Find the sign changing interval"""

        try:
            #Rac = optimize.newton(self._signTracker, guess)
            Rac = solver_mod.newton(self._signTracker, guess, step = 1e-4)
        except (solver_mod.NewtonDoneError, solver_mod.NewtonNoneError):
            pass
        except:
            raise

    def eigenvector(self):
        """Compute the eigenvectors"""

        self._gevp.solve(self.Rac, self.nev, with_vectors = True)

        return self._gevp.evp_vec[:,self.mode]

    def eigenvalue(self):
        """Compute the eigenvalues"""

        self._gevp.solve(self.Rac, self.nev)

        return self._gevp.evp_lmb[self.mode]

    def growthRate(self):
        """Get the growth rate"""

        self._gevp.solve(self.Rac, self.nev)

        return self._gevp.evp_lmb[self.mode].real

    def frequency(self):
        """Compute the frequency rate"""

        self._gevp.solve(self.Rac, self.nev, use_vector = self.mode)

        return self._gevp.evp_lmb[self.mode].imag

    def _tracker(self, Ra):
        """Track the growth rate to find critical value"""

        for i in range(0, 10):
            self._gevp.solve(Ra, self.nev, with_vectors = True, use_vector = self.mode)
            if self._gevp.evp_lmb is not None:
                break

        growth = self._gevp.evp_lmb[self.mode].real
        if self.best is None:
            self.best = (abs(growth), Ra, False)
        elif abs(growth) < self.best[0]:
            self.best = (abs(growth), Ra, True)

        return growth

    def _signTracker(self, Ra):
        """Track the growth rate to find critical value"""

        # Get change in rayleigh number
        if self._gevp.eq_params['rayleigh'] is None:
            dRa = Ra
        else:
            dRa = Ra - self._gevp.eq_params['rayleigh']

        # Solve
        self._gevp.solve(Ra, self.nev)

        if self._gevp.evp_lmb is None:
            raise solver_mod.NewtonNoneError

        growth = self._gevp.evp_lmb[self.mode].real
        if growth > 0 and growth < 100:
            self.b = Ra
        elif growth < 0 and growth > -100:
            self.a = Ra

        if self.a is not None and self.b is not None:
            raise solver_mod.NewtonDoneError

        # Check for change of sign (require small change in rayleigh and small growth rate)
        if abs(dRa/Ra) < 0.01 and abs(self._gevp.evp_lmb[self.mode].real) < 0.1:
            self._gevp.solve(Ra+dRa, self.nev)
            if self._gevp.evp_lmb is not None and growth*self._gevp.evp_lmb[self.mode].real < 0:
                self.a = min(Ra, Ra+dRa)
                self.b = max(Ra, Ra+dRa)
                raise solver_mod.NewtonDoneError

        if self.best is None:
            self.best = (abs(growth), Ra, False)
        elif abs(growth) < self.best[0]:
            self.best = (abs(growth), Ra, True)

        return growth

class GEVP:
    """Class to represent a marginal point on a marginal curve"""

    def __init__(self, model, res, eq_params, eigs, bcs, wave = None, tol = 1e-8, ellipse_radius = None, fixed_shift = False, target = None, euler = None, conv_idx = 1, spectrum_conv = 1, geometry = None, impose_symmetry = False, use_spherical_evp = False):
        """Initialize the marginal point variables"""

        self.model = copy.copy(model)
        self.res = res
        self.eq_params = eq_params.copy()
        self.eq_params['rayleigh'] = None
        self.eigs = eigs
        self.bcs = bcs.copy()
        self.no_bcs = bcs.copy()
        self.no_bcs['bcType'] = model.SOLVER_NO_TAU
        self.fields = self.model.stability_fields()
        self.evp_lmb = None
        self.evp_vec = None
        self.changed = True
        self.solver = solver_mod.GEVPSolver(tol = tol, ellipse_radius = ellipse_radius, fixed_shift = fixed_shift, target = target, euler = euler, conv_idx = conv_idx, spectrum_conv = spectrum_conv, geometry = geometry, impose_symmetry = impose_symmetry, use_spherical_evp = use_spherical_evp, whichTarget = None)
        if wave is None:
            self.wave = self.defaultWave
        else:
            self.wave = wave

    def setupProblem(self, Ra):
        """Build A and B matrices for the given rayleigh number"""
        # generate operators opA, opB, opC that args are already set by functools.partial(), but no calculation has been
        # at this stage. functions are given in the first arg in partial().

        self.eq_params['rayleigh'] = Ra

        opA = functools.partial(self.model.implicit_linear, self.res, self.eq_params, self.eigs, self.no_bcs, self.fields)
        opB = functools.partial(self.model.time, self.res, self.eq_params, self.eigs, self.no_bcs, self.fields)
        opC = functools.partial(self.model.boundary, self.res, self.eq_params, self.eigs, self.bcs, self.fields)
        sizes = self.model.stability_sizes(self.res, self.eigs, self.bcs)

        return (opA, opB, opC, sizes)

    def setEigs(self, k):
       """Set the wave number"""

       self.eigs = self.wave(k)
       self.changed = True

    def solve(self, Ra, nev, with_vectors = False, use_vector = None):
        """Solve the GEVP and store eigenvalues"""

        if self.changed or Ra != self.eq_params['rayleigh'] or (with_vectors and self.evp_vec is None):
            self.changed = False
            problem = self.setupProblem(Ra)

            if use_vector is not None and self.evp_vec is not None:
                initial_vector = self.evp_vec[use_vector]
            else:
                initial_vector = None

            if with_vectors:
                self.evp_lmb, self.evp_vec = self.solver.eigenpairs(problem, nev, initial_vector = initial_vector)
            else:
                self.evp_lmb = self.solver.eigenvalues(problem, nev, initial_vector = initial_vector)
                self.evp_vec = None
            Print("\t\t (Ra = {:.14g}, ev = ".format(Ra) + str(self.evp_lmb) + ")")

    def viewOperators(self, Ra, spy = True, write_mtx = True):
        """Spy and/or write the operators to MatrixMarket file"""

        # Create operators
        if spy or write_mtx:
            opA, opB, opC, sizes = self.setupProblem(Ra)
            restrict = self.solver.restrict_operators(sizes)
            A = opA(restriction = restrict)
            B = opB(restriction = restrict)
            C = opC(restriction = restrict)

            viewer.viewOperators(A, B, C, show = spy, save = write_mtx)

    def saveHdf5_h5py(self, spectra, geometry, postfix = ''):
        """Save spectra to HDF5 file"""

        h5_file = h5py.File('eigensolution' + postfix + '.hdf5', 'w')

        eig = int(self.eigs[0])

        # Create header
        h5_file.attrs['header'] = np.string_('StateFile'.encode('ascii'))
        if geometry == 'shell':
            h5_file.attrs['type'] = np.string_('SLFm'.encode('ascii'))
        elif geometry == 'sphere_worland':
            h5_file.attrs['type'] = np.string_('WLFm'.encode('ascii'))
        h5_file.attrs['version'] = np.string_('1.0'.encode('ascii'))

        # Create physical group
        group = h5_file.create_group('physical')
        for k,p in self.eq_params.items():
            group.create_dataset(k, (), 'f8', data = p)
        # ... add boundary conditions
        for k,b in self.bcs.items():
            if k != 'bcType':
                group.create_dataset('bc_' + k, (), 'f8', data = b)

        # Create run group
        group = h5_file.create_group('run')
        group.create_dataset('time', (), 'f8', data = -1)
        group.create_dataset('timestep', (), 'f8', data = -1)

        # Create truncation group
        base_res = [self.res[0]-1, self.res[1]-1, eig]
        group = h5_file.create_group('truncation')
        for sg in ['physical', 'spectral', 'transform']:
            subgroup = group.create_group(sg)
            file_res = base_res
            for i, d in enumerate(['dim1D', 'dim2D', 'dim3D']):
                subgroup.create_dataset(d, (), 'i4', data = file_res[i])

        # Compute dataset shape
        ls = 0
        for m in range(eig+1):
            ls = ls + self.res[1] - m
        data_shape = (ls, self.res[0])

        # Create data groups
        for k,f in spectra.items():
            name = k[0]
            group = h5_file.require_group(name)
            if k[1] != '':
                name = name + '_' + k[1]
            dset = group.create_dataset(name, data_shape, f.dtype)
            dset[-self.res[1]+eig:,:] = np.reshape(f, (self.res[1]-eig, self.res[0]), order='C')

        h5_file.close()

    def saveHdf5_tables(self, spectra, geometry, postfix = ''):
        """Save spectra to HDF5 file"""

        eig = int(self.eigs[0])

        h5_file = tables.open_file('eigensolution' + postfix + '.hdf5', mode = 'w')

        # Create header
        h5_file.set_node_attr('/', 'header', 'StateFile'.encode('ascii'))
        if geometry == 'shell':
            h5_file.set_node_attr('/', 'type', 'SLFm'.encode('ascii'))
        elif geometry == 'sphere_worland':
            h5_file.set_node_attr('/', 'type', 'WLFm'.encode('ascii'))
        h5_file.set_node_attr('/', 'version', '1.0'.encode('ascii'))

        # Create physical group
        group = h5_file.create_group('/', 'physical')
        for k,p in self.eq_params.items():
            h5_file.create_array(group, k, p)
        # ... add boundary conditions
        for k,b in self.bcs.items():
            if k != 'bcType':
                h5_file.create_array(group, 'bc_' + k, b)

        # Create run group
        group = h5_file.create_group('/', 'run')
        h5_file.create_array(group, 'time', -1)
        h5_file.create_array(group, 'timestep', -1)

        # Create truncation group
        group = h5_file.create_group('/', 'truncation')
        for sg in ['physical', 'spectral', 'transform']:
            subgroup = h5_file.create_group(group, sg)
            file_res = [self.res[0]-1, self.res[1]-1, eig]
            for i,d in enumerate(['dim1D', 'dim2D', 'dim3D']):
                h5_file.create_array(subgroup, d, file_res[i])

        # Compute dataset shape
        ls = 0
        for m in range(eig+1):
            ls = ls + self.res[1] - m
        data_shape = (ls, self.res[0])

        # Create data groups
        for g in set([k[0] for k,f in spectra.items()]):
            group = h5_file.create_group('/', g)

        # Create datasets
        for k,f in spectra.items():
            name = k[0]
            if k[1] != '':
                name = name + '_' + k[1]
            zz = np.zeros(data_shape, dtype=f.dtype)
            dset = h5_file.create_array(h5_file.get_node('/'+k[0]), name, zz)
            if eig == 0:
                f.imag = 0
            dset[-self.res[1]+eig:,:] = np.reshape(f, (self.res[1]-eig, self.res[0]), order='C')

        h5_file.close()

    def viewSpectra(self, viz_mode, plot = True, save_pdf = False):
        """Plot the spectra of the eigenvectors"""

        if self.evp_lmb is not None:
            # Gather full field across ranks
            rank = MPI.COMM_WORLD.Get_rank()
            scatter, v = PETSc.Scatter.toZero(self.evp_vec[viz_mode])
            scatter.scatter(self.evp_vec[viz_mode], v, False, PETSc.Scatter.Mode.FORWARD)
            full_vec = v.getArray()

            # Output is only done on rank 0
            if rank == 0:
                # Extra different fields
                start = 0
                stop = 0
                sol_spec = dict()
                for i,f in enumerate(self.fields):
                    stop = stop + self.model.stability_sizes(self.res, self.eigs, self.bcs)[0][i]
                    sol_spec[f] = full_vec[start:stop]
                    start = stop

                if plot or save_pdf:
                    Print("\nVisualizing spectra of mode: " + str(self.evp_lmb[viz_mode]))

                    fid = ""
                    if 'prandtl' in self.eq_params:
                        fid = fid + "Pr{:g}".format(self.eq_params['prandtl'])
                    if 'taylor' in self.eq_params:
                        fid = fid + "Ta{:g}".format(self.eq_params['taylor'])
                    viewer.viewSpectra(sol_spec, show = plot, save = save_pdf, fid = fid)

                return sol_spec

            else:
                return None

        else:
            return None

    def viewPhysical(self, viz_mode, geometry, plot = True, save_pdf = False, save_profile = True, save_hdf5 = True):
        """Plot eigenvectors in physical space"""

        if self.evp_lmb is not None:
            if viz_mode >= 0:
                # Extra different fields
                sol_spec = self.viewSpectra(viz_mode, plot = False, save_pdf = False)

                rank = MPI.COMM_WORLD.Get_rank()
                if rank == 0:
                    if self.model.use_galerkin:
                        sol_spec = self.apply_stencil(sol_spec)

                    # Save spectra to HDF5
                    if save_hdf5:
                        self.saveHdf5_tables(sol_spec, geometry)

                    Print("\nVisualizing physical data of mode: " + str(self.evp_lmb[viz_mode]))
                    fid = ""
                    if 'prandtl' in self.eq_params:
                        fid = fid + "Pr{:g}".format(self.eq_params['prandtl'])
                    if 'taylor' in self.eq_params:
                        fid = fid + "Ta{:g}".format(self.eq_params['taylor'])
                    if 'ekman' in self.eq_params:
                        fid = fid + "E{:g}".format(self.eq_params['ekman'])
                    viewer.viewPhysical(sol_spec, geometry, self.res, self.eigs, self.eq_params, show = plot, save = save_pdf, fid = fid)
            else:
                for i in range(0, len(self.evp_vec)):
                    # Extra different fields
                    sol_spec = self.viewSpectra(i, plot = False, save_pdf = False)

                    rank = MPI.COMM_WORLD.Get_rank()
                    if rank == 0:
                        if self.model.use_galerkin:
                            sol_spec = self.apply_stencil(sol_spec)

                        # Save spectra to HDF5
                        if save_hdf5:
                            self.saveHdf5_tables(sol_spec, geometry, postfix = '_mode_' + str(i))

                        Print("\nVisualizing physical data of mode: " + str(self.evp_lmb[i]))
                        fid = ""
                        if 'prandtl' in self.eq_params:
                            fid = fid + "Pr{:g}".format(self.eq_params['prandtl'])
                        if 'taylor' in self.eq_params:
                            fid = fid + "Ta{:g}".format(self.eq_params['taylor'])
                        if 'ekman' in self.eq_params:
                            fid = fid + "E{:g}".format(self.eq_params['ekman'])
                        fid += "_" + str(i)
                        viewer.viewPhysical(sol_spec, geometry, self.res, self.eigs, self.eq_params, show = plot, save = save_pdf, fid = fid)

    def apply_stencil(self, sol_vec):
        """Apply Galerkin stencil to recover physical fields"""

        for k,v in sol_vec.items():
            S = self.model.stencil(self.res, self.eq_params, self.eigs, self.bcs, k, False)
            sol_vec[k] = S*v

        return sol_vec

    def defaultWave(self, k):
        """Default conversion for wavenumber index"""

        return list(k)

def default_options():
    """Create dictionary of default options"""

    opts = dict()

    opts['fixed_shift'] = False     # Compute random shift only once
    opts['target'] = None           # Compute eigenvalues around target
    opts['euler'] = None            # Compute implicit Euler steps before starting calculation
    opts['conv_idx'] = 1            # Number of index at end of spectrum to be converged
    opts['spectrum_conv'] = 1       # Ratio between min and max of spectrum
    opts['ellipse_radius'] = None   # Restrict eigenvalue search to be within radius
    opts['mode'] = 0                # Mode to track
    opts['root_tol'] = 1e-8         # Tolerance used in root finding algorighm
    opts['evp_tol'] = 1e-10         # Tolerance used in eigensolver

    opts['geometry'] = 'c1d'        # Geometry of the system

    opts['point'] = False           # Compute marginal curve point
    opts['point_k'] = 1           # Compute marginal curve point
    opts['curve'] = False           # Trace marginal curve
    opts['curve_points'] = [1]      # Wavenumbers for marginal curve
    opts['minimum'] = False         # Compute marginal curve minimum
    opts['minimum_int'] = False     # Only used integer wave numbers
    opts['solve'] = False           # Solve fixed point
    opts['solve_nev'] = 1           # Solve fixed point

    opts['impose_symmetry'] = False # Solve fixed point
    opts['use_spherical_evp'] = False # Solve fixed point

    opts['plot_spy'] = False        # Plot matrix spy
    opts['plot_point'] = False      # Plot solution for marginal curve point
    opts['plot_curve'] = False      # Plot marginal curve
    opts['viz_mode'] = 0        # Mode to visualize/save
    opts['show_spectra'] = False    # Plot solution spectra
    opts['show_physical'] = False   # Plot solution in physical space

    opts['write_mtx'] = False       # Write matrix to MatrixMarket files
    opts['save_curve'] = False      # Save marginal curve plot to pdf
    opts['save_spectra'] = False    # Save spectra plot to pdf
    opts['save_physical'] = False   # Save physical space plots to pdf
    opts['data_spectra'] = False    # Save spectral data to file(s)
    opts['data_physical'] = False   # Save physical data to file(s)
    opts['save_hdf5'] = False   # Save physical data to file(s)

    return opts

def compute(gevp_opts, marginal_opts):
    """Compute eigen solutions requested by options"""

    # Make options consistent
    marginal_opts['minimum'] = marginal_opts['minimum'] and marginal_opts['curve']
    marginal_opts['plot_point'] = marginal_opts['plot_point'] and (marginal_opts['point'] or marginal_opts['minimum'])
    marginal_opts['plot_curve'] = marginal_opts['plot_curve'] and marginal_opts['curve']
    marginal_opts['show_spectra'] = marginal_opts['show_spectra'] and marginal_opts['solve']
    marginal_opts['show_physical'] = marginal_opts['show_physical'] and marginal_opts['solve']
    marginal_opts['save_spectra'] = marginal_opts['save_spectra'] and marginal_opts['solve']
    marginal_opts['save_physical'] = marginal_opts['save_physical'] and marginal_opts['solve']

    # Move options to GEVP
    gevp_opts['ellipse_radius'] = marginal_opts['ellipse_radius']
    gevp_opts['fixed_shift'] = marginal_opts['fixed_shift']
    gevp_opts['target'] = marginal_opts['target']
    gevp_opts['euler'] = marginal_opts['euler']
    gevp_opts['conv_idx'] = marginal_opts['conv_idx']
    gevp_opts['spectrum_conv'] = marginal_opts['spectrum_conv']
    gevp_opts['geometry'] = marginal_opts['geometry']
    gevp_opts['impose_symmetry'] = marginal_opts['impose_symmetry']
    gevp_opts['use_spherical_evp'] = marginal_opts['use_spherical_evp']

    if marginal_opts['point'] or marginal_opts['curve']:
        # Create marginal curve object
        curve = MarginalCurve(gevp_opts, mode = marginal_opts['mode'], rtol = marginal_opts['root_tol'], evp_tol = marginal_opts['evp_tol'])

    if marginal_opts['point'] and not marginal_opts['curve']:
        # Compute marginal curve at a single point
        kp = marginal_opts['point_k']
        Rac, evp_freq = curve.point(kp, guess = gevp_opts['eq_params']['rayleigh'])
        kc = kp

    if marginal_opts['curve']:
        # Trace marginal curve for a set of wave indexes
        ks = marginal_opts['curve_points']
        (data_k, data_Ra, data_freq) = curve.trace(ks, initial_guess = gevp_opts['eq_params']['rayleigh'])

        if marginal_opts['minimum']:
            # Compute minimum of marginal curve
            kc, Rac, fc = curve.minimum(data_k, data_Ra, only_int = marginal_opts['minimum_int'])

        if marginal_opts['plot_curve'] or marginal_opts['save_curve']:
            # Plot marginal curve and minimum
            if marginal_opts['minimum']:
                min_point = None
            else:
                min_point = (kc, Rac)
            curve.view(data_k, data_Ra, data_freq, minimum = min_point, plot = marginal_opts['plot_curve'], save_pdf = marginal_opts['save_curve'])

    if marginal_opts['plot_spy'] or marginal_opts['write_mtx'] or marginal_opts['solve'] or marginal_opts['minimum']:
        if marginal_opts['plot_point']:
            Ra = Rac
            kp = kc
        else:
            kp = marginal_opts['point_k']
            Ra = gevp_opts['eq_params']['rayleigh']
        Print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(kp))
        gevp_opts['eigs'] = gevp_opts['wave'](kp)
        gevp_opts['tol'] = marginal_opts['evp_tol']
        gevp = GEVP(**gevp_opts) # create solver GEVP object and its solver with gevp.solver

    if marginal_opts['plot_spy'] or marginal_opts['write_mtx']:
        gevp.viewOperators(Ra, spy = marginal_opts['plot_spy'], write_mtx = marginal_opts['write_mtx'])

    if marginal_opts['solve']:
        gevp.solve(Ra, max(marginal_opts['solve_nev'], marginal_opts['viz_mode']+3), with_vectors = True)
        Print("Found eigenvalues:")
        Print(gevp.evp_lmb)

    if marginal_opts['show_spectra'] or marginal_opts['save_spectra'] or marginal_opts['data_spectra']:
        gevp.viewSpectra(marginal_opts['viz_mode'], plot = marginal_opts['show_spectra'], save_pdf = marginal_opts['save_spectra'])

    if marginal_opts['show_physical'] or marginal_opts['save_physical'] or marginal_opts['data_physical'] or marginal_opts['save_hdf5']:
        gevp.viewPhysical(marginal_opts['viz_mode'], marginal_opts['geometry'], plot = marginal_opts['show_physical'], save_pdf = marginal_opts['save_physical'], save_hdf5 = marginal_opts['save_hdf5'])
