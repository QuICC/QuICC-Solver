# Setup SLEPc/PETSc targets via pkg-config

# In order to find PETSc and SLEPc the variables 
# PETCS_DIR, PETSC_ARCH and SLEPC_DIR, SLEPC_ARCH
# need to be set

# set root of location to find PETSc's pkg-config
set(PETSC ${PETSC_DIR}/${PETSC_ARCH})
set(SLEPC ${SLEPC_DIR}/${PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig:${SLEPC}/lib/pkgconfig)

find_package(PkgConfig REQUIRED)
pkg_search_module(SLEPC IMPORTED_TARGET slepc)

if(SLEPC_FOUND)
  set(QUICC_HAVE_SLEPC ON)
else()
  set(QUICC_HAVE_SLEPC OFF)
endif()
