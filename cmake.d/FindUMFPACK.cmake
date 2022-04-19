# FindUMFPACK
# -----------
#
# This module looks for the UMFPACK library.
#
# The following variables are set
#
# ::
#
#   UMFPACK_FOUND           - True if double precision UMFPACK library is found
#   UMFPACK_LIBRARIES       - The required libraries
#   UMFPACK_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
#   UMFPACK::UMFPACK


message(VERBOSE "Looking for UMFPACK")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# set paths to look for library
set(_UMFPACK_PATHS ${UMFPACK_ROOT} $ENV{UMFPACK_ROOT} $ENV{UMFPACKDIR})

if(_UMFPACK_PATHS)
    # disable default paths if ROOT is set
    set(_UMFPACK_DEFAULT_PATH_SWITCH NO_DEFAULT_PATH)
else()
    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PKG_UMFPACK QUIET "UMFPACK")
    endif()
    set(_UMFPACK_PATHS ${PKG_UMFPACK_LIBRARY_DIRS})
    set(_UMFPACK_INCLUDE_PATHS ${PKG_UMFPACK_INCLUDE_DIRS})
endif()


find_library(
    UMFPACK_LIBRARIES
    NAMES "umfpack"
    HINTS ${_UMFPACK_PATHS} ENV LIBRARY_PATH
    PATH_SUFFIXES "lib" "lib64"
    ${_UMFPACK_DEFAULT_PATH_SWITCH}
)
find_path(UMFPACK_INCLUDE_DIRS
    NAMES "umfpack.h"
    HINTS ${_UMFPACK_PATHS} ${_UMFPACK_INCLUDE_PATHS} ENV C_INCLUDE_PATH
    PATH_SUFFIXES "include"
    ${_UMFPACK_DEFAULT_PATH_SWITCH}
)

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK REQUIRED_VARS UMFPACK_INCLUDE_DIRS UMFPACK_LIBRARIES )

# add target to link against
if(UMFPACK_FOUND)
  message(VERBOSE "UMFPACK FOUND")
  message(VERBOSE "UMFPACK INCLUDE: ${UMFPACK_INCLUDE_DIRS}")
  message(VERBOSE "UMFPACK LIBS: ${UMFPACK_LIBRARIES}")
  set(_UMFPACK_TARGET "UMFPACK::UMFPACK")
  if(NOT TARGET ${_UMFPACK_TARGET})
      add_library(${_UMFPACK_TARGET} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_UMFPACK_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${UMFPACK_LIBRARIES})
  set_property(TARGET ${_UMFPACK_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${UMFPACK_INCLUDE_DIRS})
else()
  message(VERBOSE "UMFPACK NOT FOUND")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(UMFPACK_FOUND UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS pkgcfg_lib_PKG_UMFPACK_UMFPACK)
