# FindKaHIP
# -----------
#
# This module looks for the KaHIP library.
#
# The following variables are set
#
# ::
#
#   KAHIP_FOUND           - True KAHIP library is found
#   KAHIP_LIBRARIES       - The required libraries
#   KAHIP_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
#   KAHIP::KAHIP


message(VERBOSE "Looking for KAHIP")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# set paths to look for library
set(_KAHIP_PATHS ${KAHIP_ROOT} $ENV{KAHIP_ROOT} $ENV{KAHIPDIR})

if(_KAHIP_PATHS)
    # disable default paths if ROOT is set
    set(_KAHIP_DEFAULT_PATH_SWITCH NO_DEFAULT_PATH)
else()
    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PKG_KAHIP QUIET "KAHIP")
    endif()
    set(_KAHIP_PATHS ${PKG_KAHIP_LIBRARY_DIRS})
    set(_KAHIP_INCLUDE_PATHS ${PKG_KAHIP_INCLUDE_DIRS})
endif()


find_library(
    KAHIP_LIBRARIES
    NAMES "kahip"
    HINTS ${_KAHIP_PATHS} ENV LIBRARY_PATH
    PATH_SUFFIXES "lib" "lib64"
    ${_KAHIP_DEFAULT_PATH_SWITCH}
)
find_path(KAHIP_INCLUDE_DIRS
    NAMES "kaHIP_interface.h"
    HINTS ${_KAHIP_PATHS} ${_KAHIP_INCLUDE_PATHS} ENV C_INCLUDE_PATH
    PATH_SUFFIXES "include"
    ${_KAHIP_DEFAULT_PATH_SWITCH}
)

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KAHIP REQUIRED_VARS KAHIP_INCLUDE_DIRS KAHIP_LIBRARIES )

# add target to link against
if(KAHIP_FOUND)
  message(VERBOSE "KAHIP FOUND")
  message(VERBOSE "KAHIP VERSION: ${KAHIP_VERSION}")
  message(VERBOSE "KAHIP INCLUDE: ${KAHIP_INCLUDE_DIRS}")
  message(VERBOSE "KAHIP LIBS: ${KAHIP_LIBRARIES}")
  set(_KAHIP_TARGET "KAHIP::KAHIP")
  if(NOT TARGET ${_KAHIP_TARGET})
      add_library(${_KAHIP_TARGET} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_KAHIP_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${KAHIP_LIBRARIES})
  set_property(TARGET ${_KAHIP_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${KAHIP_INCLUDE_DIRS})
else()
  message(VERBOSE "KAHIP NOT FOUND")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(KAHIP_FOUND KAHIP_LIBRARIES KAHIP_INCLUDE_DIRS pkgcfg_lib_PKG_KAHIP_KAHIP)
