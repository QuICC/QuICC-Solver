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
#
# The following import targets are created (if requested as COMPONENTS)
#
# ::
#
#   FFTW::omp
#   FFTW::pthread
#   FFTW::mpi
#

message(VERBOSE "Looking for FTTW")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# We haven't found FFTW yet. Clear its state in case it is set in the parent
# scope somewhere else. We can't rely on it because different components may
# have been requested for this call.
set(FFTW_FOUND OFF)
set(FFTW_LIBRARIES)


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

# serial library
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
set(FFTW_REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

# Components
message(DEBUG "FFTW_FIND_COMPONENTS: ${FFTW_FIND_COMPONENTS}")
foreach(_comp IN LISTS FFTW_FIND_COMPONENTS)
    set(FFTW_${_comp}_LIBRARY)
    find_library(
        FFTW_${_comp}_LIBRARY
        NAMES "fftw3_${_comp}"
        HINTS ${_FFTW_PATHS} ENV LIBRARY_PATH
        PATH_SUFFIXES "lib" "lib64"
        ${_FFTW_DEFAULT_PATH_SWITCH}
    )
    list(APPEND FFTW_REQUIRED_VARS FFTW_${_comp}_LIBRARY)

    if (FFTW_${_comp}_LIBRARY)
        set(FFTW_${_comp}_FOUND YES)
    endif()
endforeach()

# Check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW
    REQUIRED_VARS ${FFTW_REQUIRED_VARS}
    HANDLE_COMPONENTS)


# Add target to link against
if(FFTW_FOUND)
  message(VERBOSE "FFTW INCLUDE: ${FFTW_INCLUDE_DIRS}")
  message(VERBOSE "FFTW LIBS: ${FFTW_LIBRARIES}")
  set(_main_tgt "FFTW::FFTW")
  if(NOT TARGET ${_main_tgt})
      add_library(${_main_tgt} INTERFACE IMPORTED)
  endif()
  set_property(TARGET ${_main_tgt} PROPERTY INTERFACE_LINK_LIBRARIES ${FFTW_LIBRARIES})
  set_property(TARGET ${_main_tgt} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIRS})
endif()

# Add components
foreach(_comp IN LISTS FFTW_FIND_COMPONENTS)
    if(FFTW_${_comp}_FOUND)
        set(_tgt "FFTW::${_comp}")
        if(NOT TARGET ${_tgt})
            add_library(${_tgt} INTERFACE IMPORTED)
        endif()
        set_property(TARGET ${_tgt} PROPERTY INTERFACE_LINK_LIBRARIES ${FFTW_${_comp}_LIBRARY})
        set_property(TARGET ${_tgt} PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIRS})
    endif()
    target_link_libraries(${_main_tgt} INTERFACE ${_tgt})
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# prevent clutter in cache
mark_as_advanced(FFTW_FOUND FFTW_LIBRARIES FFTW_INCLUDE_DIRS pkgcfg_lib_PKG_FFTW_fftw3)
