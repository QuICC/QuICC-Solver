Configure {#mhdPTutConfig}
========

CMake is used to generate the Makefile that will allow the compilation of the simulations. The build process uses an out-of-source approach to keep 
the sources clean. The following steps are required to create and configure the build. To fix the notations, we assume:

   - The sources are in: /home/myself/QuICC

   - The binaries should be somewhere under: /scratch/myself

   1. Create a new directory where you'll like all the binaries to be (Warning: directory might grow large!):

      #$>mkdir /scratch/myself/QuICC

   2. Binaries with different options/parallelisations need to have separate directories

      #$>mkdir /scratch/myself/QuICC/Release-Serial #(serial version)

      #$>mkdir /scratch/myself/QuICC/Release-Tubular #(tubular parallelisation)

   3. Change into one of the newly create directories,as an example Release-Tubular/:

      #$>cd /scratch/myself/QuICC/Release-Tubular

   4. Initialise CMake and start the configuration process:

      #$>ccmake /home/myself/QuICC

   5. Set the build type and the platform. Press 'c' to reconfigure. More options should now be avaiable (depending on platform)

   6. Configure the sources and press 'c' (possibly multiple times) and generate generate the Makefile by pressing 'g'

   7. Once back to the shell, the build directory should be ready to compile the sources


Available options
=================

CMAKE\_BUILD\_TYPE
------------------

   - Debug: Resulting executable will include debug symbols during compilation as well as additional debugging code 
   - Release: Generates code for production. It also desactivates a few safeguards.
   - RelWithDebInfo: Mix between Debug and Release. Useful for debugging with optimizations

CMAKE\_INSTALL\_PREFIX
----------------------

   Leave at default values unless you know what you're doing.

QUICC\_PLATFORM
--------------------
   
   The platforms files are located in code directory under cmake.d/platforms. The define platform specific information like compilers, path to libraries, etc. Below are a few examples.
   
   - Desktop
   - Fortytwo
   - Garulfo
   - Janus
   - Lothlorien
   - Theodosia

QUICC\_BOUNDARYMETHOD
-------------------------

   - Tau: The simulation will use a spectral Tau method for the boundary conditions 
   - Galerkin: The simulation will use a spectral Galerkin method where possible. This is still in a test phase.

QUICC\_COMPILER
--------------------

   Compiler to use. This option is platform specific. Below are few examples.

   - GCC
   - Intel
   - Clang

QUICC\_FFT
---------------

   FFT library to use for the transforms. This option is platform specific. There are currently two implemented interfaces. The main one is FFTW. There is also a GPU version based on cuFFT but does not (yet) provide any performance boost.

   - FFTW
   - cuFFT

QUICC\_FFTPLAN
------------------

   FFTW plan generation flag. 

   - Fast: Plan is very fast to be generated but might not be optimal.
   - Medium: Plan is slower to generate, but might be faster than the one from the Fast option. Can be useful, but might become extremely slow in MPI code.
   - Slow: Plan can be very slow to generate. It involves testing a plethora of possible algorithms for the FFT plan. Rarely useful for our FFT sizes.

QUICC\_LARGEIO
-------------------
   
   Format used for the large state files.

   - HDF5

QUICC\_LINALG
------------------

   Linear algebra library to use. This option will eventually be removed. 

   - Eigen: Use the Eigen template library.
   - (LAPACK): Obsolete option to LAPACK and will be eventually removed

QUICC\_MEMORYUSAGE
-----------------------

   Placeholder for restricted memory algorithms.

   - High:  Will always produce the fastest code but requires the most memory.
   - (Limited): Will eventually provide a slower code path to reduce memory usage.

QUICC\_MPIALGO
-------------------

   Parallelisation algorith to use. Note that not all options will work for all simulations. Compilation will fail abort if this option is not compatible with the simulation.

   - Serial: Generate a serial code with no dependence on MPI.
   - Auto:  Used for parallel mode, ought to select the best scheme depending on setup. Not yet recommended
   - Single1D: Parallelise between the first and second transform
   - Single2D: Parallelise between the second and the third transform
   - Fixed: Parallelise between the first and the second transform but with additional constrains required for problems with 2 finite directions.
   - Tubular: Pencil like parallelisation. When avaiable this will nearly always be the fastest

QUICC\_MULTPRECISION
-------------------------

   Placeholder to allow initialisation calculations to be done in arbitrary precission.

   - (ON) 
   - OFF 

QUICC\_PROFILE
-------------------
   
   Activate the performance profiler.

   - ON: Will open an additional option requesing the type of profiler: Coarse or Detailled. Coarse profiler should not have a significant overhead.
   - OFF

QUICC\_SPEIGSOLVER
-----------------------

   Obsolete option that will be removed

   - ARPACK

QUICC\_SPLINALG
-------------------

   Sparse solver. The available options are platform dependent.

   - UmfPack: UmfPack is the most likely to be available on a standard install. Resonably fast for small problems.
   - SuperLU: Usually slower than UmfPack
   - SparseLU: Eigen's own implementation of SuperLU, similar speed but not as mature yet.
   - MUMPS: Fastest of the freely avaiable solvers. Currently in test phase. Should only be used for testing.
   - [Pardiso]: Usually fast but requires a license is not well tested with this code
   - [MKLPardiso]: Part of the intel MKL library (older version of Pardiso). Usually fast, but not well tested with this code.
   - [SPQR]: Experimental QR solver. Do not use except for testing
   - [SparseQR]: Experimental QR solver. Do not use except for testing
   - [BiCGStab]: Experimental iterative solver. Do not use except for testing

QUICC\_STORAGEPROFILE
--------------------------

   Memory profiler. This is currently not providing very useful information as it hasn't been updated in a while.

   - OFF
   - (ON)

QUICC\_TEST
----------------

   Activate the TestSuite avaialbe as a separate repository QuICC\_TestSuite.

   - OFF
   - (ON)

QUICC\_TIMESTEPPER
-----------------------

   Selection of implemented timesteppers:

   - ImExRK3: Spalart's RK3 timestep scheme
   - ImSBDF2: Standard SBDF2 scheme, does currently not support a clean restart.

QUICC\_TRANSGROUPER
------------------------

   MPI parallelisation communication grouping strategy.

   - Auto: Automatic selection
   - Equation: Group communication by equation. Usually fastest for small number of CPUs.
   - Single1D: Group communication between first and second transform
   - Single2D: Group communication between second and third transform
   - Transform: Group communication after each transform. Usually fastest for large number of CPUs.
