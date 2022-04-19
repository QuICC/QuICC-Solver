Linear stability     {#mhdPTutLinStab}
================

One of the advantages in having the model matrices setup through python is that it allows to do simple linear stability calculation separately. The tools are finished and it is currently not possible to trace a marginal curve.

The linear stability tools are avaible in the Python directory under quicc/linear\_stability. In currently has a not so robust port of MATLAB's sptarn routine to compute the eigenvalues of the generalized eigenvalue.

For each implemented simulation, there is a small Python script under Python/scripts that allows to compute the eigenvalues. All the scripts follow the naming convention: trace_SimName.py.

There are differences between the different scripts, but the general structure is the same. They all allow to compute the matrices, write to a MatrixMarket file (to be read for example by MATLAB), compute eigenvalues and do some simple visualization of the eigenmodes.
