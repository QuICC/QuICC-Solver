#
# Utilities to add libraries that depend on Kokkos
#
# Usage:
#
# # Somewhere upstream
# include(cmake.d/setup/Kokkos.cmake)
# ...
# # At the library level
# include(cmake.d/setup/KokkosAddLibrary.cmake)
# quicc_add_library_kokkos(my_lib_kk src1 src2 ...)
# ...
# # Alternatively
# include(cmake.d/setup/KokkosAddLibrary.cmake)
# quicc_add_library_kokkos(my_lib_kk "")
# add_subdirectory(src)
#     ... inside src dir ...
#     quicc_target_sources_kokkos(my_lib_kk src1 src2 ...)
#     ... out of src dir ...
# quicc_fix_target_sources_kokkos(my_lib_kk)
#

#
# Add library that depends on Kokkos.
# If Kokkos uses CUDA it requires special handling
# in order to compile with nvcc
#
function(quicc_add_library_kokkos libname)
  # check precondition
  # QuICC::Kokkos target must be defined
  if(NOT TARGET QuICC::Kokkos)
    message(FATAL_ERROR "Kokkos was not setup")
  endif()


  # parse inputs
  # cmake_parse_arguments(QAKL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_library_kokkos")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  # message(DEBUG "QAKL_SRC: ${QAKL_SRC}")
  message(DEBUG "libname: ${libname}")
  message(DEBUG "ARGN: ${ARGN}")

  message(DEBUG "Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")
  if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(${ARGN} PROPERTIES LANGUAGE CUDA)
  endif()
  add_library(${libname} STATIC ${ARGN})
  target_link_libraries(${libname} PRIVATE QuICC::Kokkos)


  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_add_library_kokkos)

#
# Add sources to existing library that depends on Kokkos.
# If this is called in a subdirectory with respect to
# quicc_add_library_kokkos, the quicc_fix_target_sources_kokkos
# needs to be called as well.
#
function(quicc_target_sources_kokkos libname)

  message(DEBUG "quicc_target_sources_kokkos")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "libname: ${libname}")
  message(DEBUG "ARGN: ${ARGN}")

  message(DEBUG "Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")
  if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(${ARGN} PROPERTIES LANGUAGE CUDA)
  endif()

  target_sources(${libname} PRIVATE ${ARGN})

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_target_sources_kokkos)

#
# Unfortunately, set_source_files_properties is directory scoped.
# This function is a patch that needs to be called from the
# directory where quicc_add_library_kokkos is called
# and AFTER adding the target sources
#
# https://gitlab.kitware.com/cmake/cmake/-/issues/20128
#
function(quicc_fix_target_sources_kokkos libname)

  message(DEBUG "quicc_fix_target_sources_kokkos")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "libname: ${libname}")
  message(DEBUG "ARGN: ${ARGN}")

  # get sources
  get_target_property(_src ${libname} SOURCES)
  message(DEBUG "_src: ${_src}")

  message(DEBUG "Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")
  if(Kokkos_ENABLE_CUDA)
    set_source_files_properties(${_src} PROPERTIES LANGUAGE CUDA)
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction(quicc_fix_target_sources_kokkos)
