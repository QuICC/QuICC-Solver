Linear stability solver
=======================

The Stability component provides as simple linear stability solver to be used together with the C++ backend implementations of the Models.
The underlying generalized eigenvalue problem is solved using PETSc/SLEPc.

General setup
=============

PETSc and SLEPc needs to be available in order to use the linear stability solver. The location of PETSc/SLEPc needs to be provided to cmake in order to enable the stability solver:

```bash
#$>cmake . -DPETSC_DIR=/path/to/petsc -DPETSC_ARCH="MyPETScArch" -DSLEPC_DIR=/path/to/slepc -DSLEPC_ARCH="MySLEPcArch"
```

In the case of a prefix install of PETSc/SLEPc the two `_ARCH` options can be omited. If the configuration is successful the configured models should have a new "Stability" executable. For example if the BoussinesqSphereRTC model is configured, the stability solver can be build using:

```bash
#$>make -j 6 BoussinesqSphereRTCImplicitStability
```

Model setup
===========

As part of the linear stability solver, the "useLinearized()" flag is enabled. In order to provide the linearized system the model can be modified by using this flag as a guard. See for example the BoussinesqSphereRTC model and search for `useLinearized()`.

MPI vs Serial PETSc/SLEPc
=========================

It is generally easier to setup an MPI version of PETSc/SLEPc. For a serial execution, the linear stability solver can be used with a single rank.

Problem setup
=============

The parameters for the physical model are read from the `parameters.cfg` file. The Rayleigh number from the file is used as initial guess. The interpretation of the resolution depends on the problem. For a spherical setup, `dim1D` and `dim2D` are the radial and harmonic degree truncation as usual. `dim3D` is the harmonic order of the problem to solve.

In addition to the parameters for the nonlinear simulations, the following parameters are added:
- omega: the target frequency for the GEVP solver
- growth_rate: the target growth rate for the GEVP solver
- tolerance: tolerance for converged eigenvalues
- max_iteration: max iterations for eigensolver
- nev: the number of eigenvalues to compute
- sort: sorting to apply to the eigenpairs. This is applied after the eigenvalues have been computed.
    - 0: no additional sorting
    - 1: decreasing real part
- stability_mode: computation mode for the solver:
    - 0: solve GEVP with provided parameters
    - 1: solve GEVP with provided parameters and save eigenfunction(s)
    - 2: compute critical Rayleigh number starting from provided initial guess
- write_mtx: write matrices to MatrixMarket format

# Legacy python base margincal curve tracing tool

QuICC provides a set of Python scripts to do linear stability calculations. The tools require PETSc and SLEPc, as well as their Python bindings petsc4py and slepc4py. As those library are quite complex and provide a lot of different options during installation, obtaining the right setup can be challenging. In order to simplify the setup, a docker image which contains all dependencies is provided. The image is available on the docker hub as "quicc/stability". In order to use it, follow the instructions above to setup docker but use "quicc/stability" in place of "quicc/buildbase". 

## Run the marginal curve tracing tool

Let's assume the configured build directory is called `build_dir` and the model used is `BoussinesqSphereRTC`. To trace a marginal curve do the following:
 1. Change to build directory
    ```bash
    cd /path/to/build_dir
    ```
 2. Copy linear stability files to build directory
    ```bash
    make quicc_pyquicc quicc_boussinesqspherertc_updatepy
    ```
 3. The scripts are installed under `lib64/python/linear_stability`
    ```bash
    cd lib64/python/linear_stability/boussinesq/sphere/rtc/
    ```
 4. Define the parameters and options in `trace_marginal.py`. The script contains a few examples for different parameters, resolutions, etc. At the bottom of the file are the main options for the tracing tool. The options are described in `lib64/python/quicc/linear_stability/marginal_curve.py` in `default_options` line 597. The main options are
    - `marginal_options['curve']`: True/False trace marginal curve
    - `marginal_options['minimum']`: True/False compute minimum
    - `marginal_options['solve']`: True/False compute single point
    - `marginal_options['solve_nev']`: Number of eigenvalues to compute
    - `marginal_options['curve_points']`: list of modes where to compute the marginal curve
 5. Tell Python where to find the modules
    ```bash
    export PYTHONPATH=/path/to/build_dir/lib64/python
    ```
 6. Run the marginal curve tracing script:
    ```bash
    python trace_marginal.py -st_type sinvert -eps_target 0 -eps_target_real
    ```
