# FindMETIS
# -----------
#
# This module looks for the METIS library.
#
# The following variables are set
#
# ::
#
#   METIS_FOUND           - True METIS library is found
#   METIS_LIBRARIES       - The required libraries
#   METIS_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
#   METIS::METIS


message(VERBOSE "Looking for METIS")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# set paths to look for library
set(_METIS_PATHS ${METIS_ROOT} $ENV{METIS_ROOT} $ENV{METISDIR})

if(_METIS_PATHS)
    # disable default paths if ROOT is set
    set(_METIS_DEFAULT_PATH_SWITCH NO_DEFAULT_PATH)
else()
    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PKG_METIS QUIET "METIS")
    endif()
    set(_METIS_PATHS ${PKG_METIS_LIBRARY_DIRS})
    set(_METIS_INCLUDE_PATHS ${PKG_METIS_INCLUDE_DIRS})
endif()


find_library(
    METIS_LIBRARIES
    NAMES "metis"
    HINTS ${_METIS_PATHS} ENV LIBRARY_PATH
    PATH_SUFFIXES "lib" "lib64"
    ${_METIS_DEFAULT_PATH_SWITCH}
)
find_path(METIS_INCLUDE_DIRS
    NAMES "metis.h"
    HINTS ${_METIS_PATHS} ${_METIS_INCLUDE_PATHS} ENV C_INCLUDE_PATH
    PATH_SUFFIXES "include"
    ${_METIS_DEFAULT_PATH_SWITCH}
)

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS REQUIRED_VARS METIS_INCLUDE_DIRS METIS_LIBRARIES )

# add target to link against
if(METIS_FOUND)
  message(VERBOSE "METIS FOUND")
  message(VERBOSE "METIS VERSION: ${METIS_VERSION}")
  message(VERBOSE "METIS INCLUDE: ${METIS_INCLUDE_DIRS}")
  message(VERBOSE "METIS LIBS: ${METIS_LIBRARIES}")
  set(_METIS_TARGET "METIS::METIS")
  if(NOT TARGET ${_METIS_TARGET})
      add_library(${_METIS_TARGET} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_METIS_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${METIS_LIBRARIES})
  set_property(TARGET ${_METIS_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS})
else()
  message(VERBOSE "METIS NOT FOUND")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(METIS_FOUND METIS_LIBRARIES METIS_INCLUDE_DIRS pkgcfg_lib_PKG_METIS_METIS)
