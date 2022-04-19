Platform {#mhdPTutPlatform}
========

The platform files located under *cmaked.d/platforms/* are simple scripts defining platform specific parameters.

More precisely it defines:

   - the list of available compilers

   - the list of available implementation of the smart pointers

   - the list of available FFT implementations and the corresponding library name to link with, as well as list of possible plans

   - the list of available linear algebra libraries

   - the list of available linear sparse solver as well as the corresponding library name to link width and include paths

   - the list of available state file formats. Only HDF5 is currently supported

   - the list of availabe multiple precision libraries. Currently not used in the code.

   - compiler settings and optimization flags as well as corresponding MPI compiler settings

   - Python library and include paths as well as version to use
