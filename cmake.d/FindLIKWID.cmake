# FindLIKWID
# -----------
#
# This module looks for the LIKWID library.
#
# The following variables are set
#
# ::
#
#   LIKWID_FOUND           - True if double precision LIKWID library is found
#   LIKWID_LIBRARIES       - The required libraries
#   LIKWID_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
#   LIKWID::LIKWID


message(VERBOSE "Looking for LIKWID")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# set paths to look for library
set(_LIKWID_PATHS ${LIKWID_ROOT} $ENV{LIKWID_ROOT} $ENV{LIKWIDDIR})

if(_LIKWID_PATHS)
    # disable default paths if ROOT is set
    set(_LIKWID_DEFAULT_PATH_SWITCH NO_DEFAULT_PATH)
else()
    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PKG_LIKWID QUIET "LIKWID")
    endif()
    set(_LIKWID_PATHS ${PKG_LIKWID_LIBRARY_DIRS})
    set(_LIKWID_INCLUDE_PATHS ${PKG_LIKWID_INCLUDE_DIRS})
endif()


find_library(
    LIKWID_LIBRARIES
    NAMES "likwid"
    HINTS ${_LIKWID_PATHS} ENV LIBRARY_PATH
    PATH_SUFFIXES "lib" "lib64"
    ${_LIKWID_DEFAULT_PATH_SWITCH}
)
find_path(LIKWID_INCLUDE_DIRS
    NAMES "likwid.h"
    HINTS ${_LIKWID_PATHS} ${_LIKWID_INCLUDE_PATHS} ENV C_INCLUDE_PATH
    PATH_SUFFIXES "include"
    ${_LIKWID_DEFAULT_PATH_SWITCH}
)

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIKWID REQUIRED_VARS LIKWID_INCLUDE_DIRS LIKWID_LIBRARIES )

# add target to link against
if(LIKWID_FOUND)
  message(VERBOSE "LIKWID FOUND")
  message(VERBOSE "LIKWID VERSION: ${LIKWID_VERSION}")
  message(VERBOSE "LIKWID INCLUDE: ${LIKWID_INCLUDE_DIRS}")
  message(VERBOSE "LIKWID LIBS: ${LIKWID_LIBRARIES}")
  set(_LIKWID_TARGET "LIKWID::LIKWID")
  if(NOT TARGET ${_LIKWID_TARGET})
      add_library(${_LIKWID_TARGET} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_LIKWID_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${LIKWID_LIBRARIES})
  set_property(TARGET ${_LIKWID_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${LIKWID_INCLUDE_DIRS})
else()
  message(VERBOSE "LIKWID NOT FOUND")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(LIKWID_FOUND LIKWID_LIBRARIES LIKWID_INCLUDE_DIRS pkgcfg_lib_PKG_LIKWID_LIKWID)
