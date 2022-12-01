# Quasi-Inverse Convection Code

### QuICC is a simulation framework for solving the Navier-Stokes equations in spherical, cylindrical and Cartesian geometry.

## Features:

   - Fully spectral using Fourier, Chebyshev, Worland and Associated Legendre basis
   - Multiple geometries: spherical shell, full sphere, cylinder, plane layer, cartesian box
   - Sparse formulation of the linear operators
   - High-level Python interface for the linear operators
   - Nonlinear simulations
   - Linear stability calculations using PETSc/SLEPc
   - 1D, 2D and 2.5D parallelisation

## Quick Start

0. Upload your SSH key to your github profile
1. Clone repository
   ```bash
   git clone -b main git@github.com:QuICC/QuICC.git
   ```

2. Create build directory
   ```bash
   mkdir </path/to/Builds>
   cd </path/to/Builds>
   ```
3. Configure build with CMake specifying the model(s) to build
   ```bash
   cmake </path/to/QuICC> -DQUICC_MODEL=<GreatSimulation1>\;<GreatSimulation2>
   ```
   The physical model descriptions are stored in separate git repositories. A list of existing models is obtained [here](https://github.com/QuICC?q=Model-+in%3Aname&type=&language=). All official model repositories follow the naming "Model-<GreatSimulation>".

   The registered models are cloned into Models/.

   Custom models can be added into Models/.

4. Compile the model executables
   ```bash
   make <GreatSimulation>
   ```
5. Create configuration XML
   ```bash
   ./<GreatSimulation>Config
   mv parameters_GEOMETRY.cfg parameters.cfg
   edit parameters.cfg
   ```
6. Create initial state
   ```bash
   ./<GreatSimulation>State
   mv state0000.hdf5 state_initial.hdf5
   ```
7. Run simulation
   ```bash
   ./<GreatSimulation>Model
   ```
8. Create physical space data for visualization
   ```bash
   ln -s state0042.hdf5 state4Visu.hdf5
   ./<GreatSimulation>Visu
   ```
9. visualize *visState0000.hdf5*

## Benchmarks

A set of benchmarks is provided in the [Benchmarks](https://github.com/QuICC/Benchmarks) repository.


## Tests

Tests are provided for some of the components.
First, one need to enable the compilation, for instance
```bash
cmake </path/to/QuICC> -DQUICC_TESTSUITE_<component>=ON
```
then, after compiling, one can run all the available tests

```bash
ctest --output-on-failure
```
or a specific test
```bash
ctest -R <only-this-test> --output-on-failure
```

## Code completion in Vim

Copy or link the compilation database "compilation_commands.json" into the root of QuICC in order to use [YouCompleteMe](https://github.com/ycm-core/YouCompleteMe)

The database can be generated in a build directory using Clang++ as compiler with
#$> cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .

## Custom models
The default repository and branch can be changed, for instance
```bash
cmake .. -DQUICC_MODEL=<myModel> -DQUICC_MODEL_GIT_REPO_BASE=<myRepo> -DQUICC_MODEL_GIT_BRANCH=<myBranch>
```
Remember that the repository full name should be in the form `<myRepo><myModel>`!

Additionally, they can be specified on a per model basis, for instance
```bash
cmake .. -DQUICC_MODEL=<myModel1>\;<myModel2> -DQUICC_MODEL_GIT_BRANCH_<myModel1>=<myBranch>
```

## Platform specific builds
To check the options selected, use the verbose mode:
```bash
cmake </path/to/QuICC> --log-level=VERBOSE
```

For instructions about known clusters refer to `README_clusters.md`.


### General optimization
It is always better to compile for the target architecture (this is not necessary on Piz-Daint as CC will do it for you)
```bash
cmake </path/to/QuICC> -DCMAKE_CXX_FLAGS="-march=native"
```
and to enable vectorized code generation in the `Eigen` library
```bash
cmake </path/to/QuICC> -DQUICC_EIGEN_ENABLE_VECTORIZATION=ON
```

## Forcing specific libraries

### Boost
```bash
cmake </path/to/QuICC> -DBOOST_ROOT=</path/to/boost>
```

### Python
```bash
cmake </path/to/QuICC> -DPython_ROOT_DIR=</path/to/python>
```

## Known issues

### MKL
MKL has a FFTW interface that, if loaded in front, can create issues by intercepting calls for FFTW. This is not a problem in itself except that MKL does not provide all the plans that are used in QuICC.
One could still use MKL for the BLAS back-end and ensure that the FFTW shared object is loaded first as follows:
```bash
LD_PRELOAD=</path/to/FFTW>/libfftw3.so <executable>
```


## Debugging

### gdb hook
By setting the enviroment variable `QUICC_GDB_HOOK` to a non-negative integer value, that process will print its PID and wait
```bash
QUICC_GDB_HOOK=0 mpirun -n 4 <executable>
```
## Profiling

Enabled by selecting the desired back-end and granularity level, for instance
```bash
cmake </path/to/QuICC> -DQUICC_PROFILE_BACKEND=native -DQUICC_PROFILE_LEVEL=0
```
### barrier
By setting the enviroment variable `QUICC_PROFILE_MPI_BARRIER` to either `after` and/or `before` a mpi barrier will be set before and/or after each region
