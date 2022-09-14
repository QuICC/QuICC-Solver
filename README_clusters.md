# Piz-Daint on GPU nodes

```bash
module load daint-gpu
module switch PrgEnv-cray PrgEnv-gnu
module load cray-fftw cray-hdf5-parallel cray-python cray-tpsl Boost
module use /apps/daint/UES/eurohack/modules/all
module load CMake/3.21.2

cmake </path/to/QuICC> -DCMAKE_CXX_COMPILER=CC \
-DQUICC_USE_MPI=ON \
-DQUICC_MULTPRECISION=ON \
-DQUICC_EIGEN_ENABLE_VECTORIZATION=ON \
-DQUICC_MODEL=<GreatSimulation>

make <GreatSimulation><Implementation>
```

# Piz-Daint on MC nodes

```bash
module load daint-mc
module switch PrgEnv-cray PrgEnv-gnu
module load cray-fftw cray-hdf5-parallel cray-python cray-tpsl Boost
module use /apps/daint/UES/eurohack/modules/all
module load CMake/3.21.2

cmake </path/to/QuICC> -DCMAKE_CXX_COMPILER=CC \
-DQUICC_USE_MPI=ON \
-DQUICC_MULTPRECISION=ON \
-DQUICC_EIGEN_ENABLE_VECTORIZATION=ON \
-DQUICC_MODEL=<GreatSimulation>

make <GreatSimulation><Implementation>
```

# Euler

```bash
env2lmod
module load cmake/3.20.3 gcc/8.2.0 openmpi openblas fftw hdf5 boost python

cmake </path/to/QuICC> -DQUICC_USE_MPI=ON \
-DQUICC_MULTPRECISION=ON \
-DQUICC_MODEL=<GreatSimulation>

make <GreatSimulation>
```
the default `h5py` is built with the incorrect hdf5 library, we need to install it ourselves

```bash
python -m venv quicc_env
. quicc_env/bin/activate
python -m pip install h5py
```
afterwards you'll only need to activate the python env
```bash
. quicc_env/bin/activate
```
