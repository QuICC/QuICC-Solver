###################################################
#--------------------- COMPILER ------------------#
###################################################

message(STATUS "Global setup")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

#
# MPI
#
find_package(MPI)


#
# BLAS
#

# look for BLAS, MKL specifically for the 32 bit version
set(BLA_VENDOR "Intel10_64lp")
find_package(BLAS)
if(BLAS_FOUND)
  message(VERBOSE "MKL BLAS")
  set(_BLAS_MKL "True")
else()
  message(VERBOSE "Generic BLAS")
  unset(BLA_VENDOR)
  find_package(BLAS REQUIRED)
endif()

# add target if needed
set(_BLAS_TARGET "BLAS::BLAS")
if(NOT TARGET ${_BLAS_TARGET})
  message(VERBOSE "BLAS target is missing")
  add_library(${_BLAS_TARGET} INTERFACE IMPORTED)
else()
message(VERBOSE "BLAS target exists")
# check tgt properties
get_target_property(BLAS_INTERFACE_LINK_LIBS ${_BLAS_TARGET} INTERFACE_LINK_LIBRARIES)
message(DEBUG "BLAS_INTERFACE_LINK_LIBS: ${BLAS_INTERFACE_LINK_LIBS}")
get_target_property(BLAS_INTERFACE_INCLUDE_DIRECTORIES ${_BLAS_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
  message(DEBUG "BLAS_INTERFACE_INCLUDE_DIRECTORIES: ${BLAS_INTERFACE_INCLUDE_DIRECTORIES}")
  # set tgt properties
  set_property(TARGET ${_BLAS_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
  set_property(TARGET ${_BLAS_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS})
endif()

if(_BLAS_MKL)
    target_compile_definitions(${_BLAS_TARGET} INTERFACE "QUICC_USING_MKL_BLAS")
endif()

#
# Eigen
#
include(BundleEigen)

#
# Kokkos
#
include(setup/Kokkos)

#
# Cuda
#
find_package(CUDAToolkit 11.3)

###################################################
#------------ THREADS PARALLELISATION ------------#
###################################################

quicc_create_option(NAME QUICC_THREADS
                    OPTS none pthread openmp
                    LABEL "Threads paralellization")
quicc_add_definition(QUICC_THREADS)


###################################################
#------------- FFT PLAN COMPUTATION --------------#
###################################################

quicc_create_option(
  NAME QUICC_FFTW_THREADS
  OPTS none pthread omp
  LABEL "FFTW threading")

if(QUICC_FFTW_THREADS STREQUAL "none")
  find_package(FFTW)
else()
  find_package(FFTW COMPONENTS ${QUICC_FFTW_THREADS} REQUIRED)
endif()
quicc_create_option(NAME QUICC_FFTPLAN
                    OPTS Fast Medium Slow
                    LABEL "FFT plan")

###################################################
#------------ VKFFT FFT IMPLEMENTATION -----------#
###################################################

if(QUICC_USE_VKFFT)
  include(BundleVkFFT)
endif()

###################################################
#------------ PFSOLVE JWT IMPLEMENTATION -----------#
###################################################
option(QUICC_USE_PFSOLVE "Enable PfSolve backend for JWT" OFF)

if(QUICC_USE_PFSOLVE)
  include(BundlePfSolve)
endif()
###################################################
#----- SPARSE LINEAR ALGEBRA IMPLEMENTATION ------#
###################################################

quicc_create_option(NAME QUICC_SPLINALG
                    OPTS SparseLU MUMPS UmfPack
                    LABEL "Sparse linear algebra")
quicc_add_definition(QUICC_SPLINALG)

if(QUICC_SPLINALG STREQUAL "MUMPS")
  option(QUICC_MPISPSOLVE "Use MPI sparse solver?" OFF)
  if(QUICC_MPISPSOLVE)
    add_definitions("-DQUICC_MPISPSOLVE")
  endif(QUICC_MPISPSOLVE)
endif(QUICC_SPLINALG STREQUAL "MUMPS")


###################################################
#---- SPARSE SPD LINEAR ALGEBRA IMPLEMENTATION ---#
###################################################

quicc_create_option(NAME QUICC_SPSPDLINALGS
                    OPTS SimplicialLDLT SimplicialLLT SparseLU UmfPack MUMPS
                    LABEL "Sparse SPD linear algebra"
                    ADVANCED)
quicc_add_definition(QUICC_SPSPDLINALGS)

###################################################
#---- SPARSE TRI LINEAR ALGEBRA IMPLEMENTATION ---#
###################################################

quicc_create_option(NAME QUICC_SPTRILINALG
                    OPTS SparseLU UmfPack MUMPS
                    LABEL "Sparse triangular linear algebra"
                    ADVANCED)
quicc_add_definition(QUICC_SPTRILINALG)


# Look for linear algebra libraries that are in use
include(MatchAny)

match_any(NAME "QUICC_" STRING "UmfPack")
if(UmfPack_IS_USED)
  find_package(UMFPACK REQUIRED)
endif()

match_any(NAME "QUICC_" STRING "MUMPS")
if(MUMPS_IS_USED)
  message(SEND_ERROR "MUMPS not supported yet.")
endif()

###################################################
#------------------ LARGE IO FORMAT --------------#
###################################################

if(MPI_FOUND)
  set(HDF5_PREFER_PARALLEL "True")
endif()

find_package(HDF5)
include(MarkAsAdvancedAll)
mark_as_advanced_all(HDF5)

if(HDF5_FOUND)
  #TODO, check parallel vs serial in accordance to QUICC_USE_MPI
  set(QUICC_LARGEIO "HDF5")
  quicc_add_definition(QUICC_LARGEIO)
  # add modern target if it does not exist (pre CMake 3.19)
  if(NOT TARGET hdf5::hdf5)
    add_library(hdf5::hdf5 INTERFACE IMPORTED)
    set_property(TARGET hdf5::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
    set_property(TARGET hdf5::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
  endif()
endif()


###################################################
#-------------- HDF5 COMPLEX FORMAT --------------#
###################################################

quicc_create_option(NAME QUICC_HDF5_CMPLX
                    OPTS "Array" "Struct"
                    LABEL "HDF5 complex format"
                    ADVANCED)
quicc_add_definition(QUICC_HDF5_CMPLX)


###################################################
#----------------- SH NORMALIZATION --------------#
###################################################

quicc_create_option(NAME QUICC_SH_NORM
                    OPTS "Unity" "Schmidt"
                    LABEL "Spherical harmonics normalization"
                    ADVANCED)
quicc_add_definition(QUICC_SH_NORM)


###################################################
#---------- TRANSFORM TREE OPTIMIZATION ----------#
###################################################

# Disable by default as it doesn't work for all cases
option(QUICC_OPTIMIZE_TREE "Optimize transform tree?" ON)
mark_as_advanced(FORCE QUICC_OPTIMIZE_TREE)
if(QUICC_OPTIMIZE_TREE)
  add_definitions("-DQUICC_OPTIMIZE_TREE")
endif(QUICC_OPTIMIZE_TREE)


###################################################
#------------------ EMBEDDED PYTHON --------------#
###################################################

find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)

###################################################
#-------------------- SLEPc/PETSc ----------------#
###################################################
include(cmake.d/PkgConfigSLEPc.cmake)

# Enable stability solver if SLEPc/PETSc are available
if(QUICC_HAVE_SLEPC)
  message(STATUS "SLEPc is present. Enabling stability solver.")

  # Enable stability solver
  set(QUICC_HAVE_STABILITY_SOLVER ON CACHE INTERNAL "Enable stability solver")

# Disable stability solver if SLEPc/PETSc is not available
else()
  set(QUICC_HAVE_STABILITY_SOLVER OFF CACHE INTERNAL "Enable stability solver")

  message(STATUS "SLEPc was not found. Stability solver will not be available.")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
