#  Copyright (c) 2019 ETH Zurich, Simon Frasch
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.


#.rst:
# FindFFTW
# -----------
#
# This module looks for the fftw3 library.
#
# The following variables are set
#
# ::
#
#   FFTW_FOUND           - True if double precision fftw library is found
#   FFTW_LIBRARIES       - The required libraries
#   FFTW_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
#   FFTW::FFTW


message(VERBOSE "Looking for FTTW")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# set paths to look for library
set(_FFTW_PATHS ${FFTW_ROOT} $ENV{FFTW_ROOT})

if(_FFTW_PATHS)
    # disable default paths if ROOT is set
    set(_FFTW_DEFAULT_PATH_SWITCH NO_DEFAULT_PATH)
else()
    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(PKG_FFTW QUIET "fftw3")
    endif()
    set(_FFTW_PATHS ${PKG_FFTW_LIBRARY_DIRS})
    set(_FFTW_INCLUDE_PATHS ${PKG_FFTW_INCLUDE_DIRS})
endif()


find_library(
    FFTW_LIBRARIES
    NAMES "fftw3"
    HINTS ${_FFTW_PATHS} ENV LIBRARY_PATH 
    PATH_SUFFIXES "lib" "lib64"
    ${_FFTW_DEFAULT_PATH_SWITCH} 
)
find_path(FFTW_INCLUDE_DIRS
    NAMES "fftw3.h"
    HINTS ${_FFTW_PATHS} ${_FFTW_INCLUDE_PATHS} ENV C_INCLUDE_PATH
    PATH_SUFFIXES "include" "include/fftw"
    ${_FFTW_DEFAULT_PATH_SWITCH}
)

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW REQUIRED_VARS FFTW_INCLUDE_DIRS FFTW_LIBRARIES )

# add target to link against
if(FFTW_FOUND)
  message(VERBOSE "FFTW FOUND")
  message(VERBOSE "FFTW INCLUDE: ${FFTW_INCLUDE_DIRS}")
  message(VERBOSE "FFTW LIBS: ${FFTW_LIBRARIES}")
  set(_FFTW_TARGET "FFTW::FFTW")
  if(NOT TARGET ${_FFTW_TARGET})
      add_library(${_FFTW_TARGET} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_FFTW_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES ${FFTW_LIBRARIES})
  set_property(TARGET ${_FFTW_TARGET} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIRS})
else()
  message(VERBOSE "FFTW NOT FOUND")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(FFTW_FOUND FFTW_LIBRARIES FFTW_INCLUDE_DIRS pkgcfg_lib_PKG_FFTW_fftw3)
