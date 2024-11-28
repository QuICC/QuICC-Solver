Linear stability solver
=======================

This component provides as simple linear stability solver to be used together with the C++ backend implementations of the Models.
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
- nev: the number of eigenvalues to compute
- sort: sorting to apply to the eigenpairs. This is applied after the eigenvalues have been computed.
    - 0: no additional sorting
    - 1: decreasing real part
- stability_mode: computation mode for the solver:
    - 0: solve GEVP with provided parameters and save eigenfunction(s)
    - 1: compute critical Rayleigh number starting from provided initial guess
