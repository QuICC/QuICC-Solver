#
# This is a somewhat complex setup to import the Kokkos library.
# The goal is to compartimentalize compilation units with Kokkos.
# The latters will be compiled with CUDA compiler while the rest
# of the compilation units are compiled with the CXX compiler.
#
# We define a target QuICC::Kokkos with the appropriate properties.
# All libraries that depends on Kokkos need to be added as follows
#
#      include(KokkosAddLibrary)
#      quicc_add_kokkos_library(libname src1 src2 ...)
#

option(QUICC_USE_KOKKOS "Enable Kokkos" OFF)
include(CMakeDependentOption)
cmake_dependent_option(QUICC_USE_KOKKOS_CUDA "Enable CUDA backend for Kokkos" OFF
    "QUICC_USE_KOKKOS" OFF)
cmake_dependent_option(QUICC_USE_KOKKOS_HIP "Enable HIP backend for Kokkos" OFF
    "QUICC_USE_KOKKOS" OFF)
cmake_dependent_option(QUICC_USE_KOKKOS_CUDAUVM "Kokkos use as default memory space CUDAUVM" OFF
    "QUICC_USE_KOKKOS;QUICC_USE_KOKKOS_CUDA" OFF)
cmake_dependent_option(QUICC_USE_KOKKOS_KERNELS "Enable Kokkos Kernels" ON
    "QUICC_USE_KOKKOS" OFF)

if(QUICC_USE_KOKKOS)
    message(STATUS "Setup Kokkos")
    list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

    if(NOT TRILINOS_DIR)
        message(DEBUG "Setting TRILINOS_DIR to $ENV{TRILINOS_DIR}")
        set(TRILINOS_DIR
            $ENV{TRILINOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    if(NOT KOKKOS_DIR)
        message(DEBUG "Setting KOKKOS_DIR to $ENV{KOKKOS_DIR}")
        set(KOKKOS_DIR
            $ENV{KOKKOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    find_package(
        Kokkos
        HINTS
        ${KOKKOS_DIR}
        ${KOKKOS_DIR}/lib/CMake/Kokkos
        ${KOKKOS_DIR}/lib64/CMake/Kokkos
        ${TRILINOS_DIR}
        ${TRILINOS_DIR}/lib/cmake/Kokkos
        ${TRILINOS_DIR}/lib64/cmake/Kokkos
        REQUIRED)

    message(VERBOSE "Found Kokkos")

    # check what was found
    message(VERBOSE "Kokkos_CXX_FLAGS: ${Kokkos_CXX_FLAGS}")
    message(VERBOSE "Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
    message(VERBOSE "Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
    message(VERBOSE "Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
    message(VERBOSE "Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
    message(VERBOSE "Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")

    # _KK_TARGET is set as a local variable
    # do not use outside this file
    set(_KK_TARGET "Kokkos::kokkos")

    if(Kokkos_ENABLE_OPENMP)
      set(_openmp "-fopenmp")
      # we need to be sure that all targets link against opemp
      add_link_options(${_openmp})
    endif()

    if(NOT TARGET ${_KK_TARGET})
        message(DEBUG "Kokkos target is not defined")
        add_library(${_KK_TARGET} INTERFACE IMPORTED)
        set_target_properties(${_KK_TARGET}
            PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES ${Kokkos_INCLUDE_DIRS} ${Kokkos_TPL_INCLUDE_DIRS}
                INTERFACE_LINK_LIBRARIES ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES}
                INTERFACE_LINK_DIRECTORIES ${Kokkos_LIBRARY_DIRS}
                INTERFACE_COMPILE_OPTIONS ${_openmp}
        )
    else()
        message(DEBUG "Kokkos target is defined")
    endif()

    # Check what the (imported) target does
    get_target_property(Kokkos_INTERFACE_COMPILE_OPTIONS ${_KK_TARGET} INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Kokkos_INTERFACE_COMPILE_OPTIONS: ${Kokkos_INTERFACE_COMPILE_OPTIONS}")
    get_target_property(Kokkos_INTERFACE_LINK_LIBRARIES ${_KK_TARGET} INTERFACE_LINK_LIBRARIES)
    message(DEBUG "Kokkos_INTERFACE_LINK_LIBRARIES: ${Kokkos_INTERFACE_LINK_LIBRARIES}")
    get_target_property(Kokkos_INTERFACE_INCLUDE_DIRECTORIES ${_KK_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    message(DEBUG "Kokkos_INTERFACE_INCLUDE_DIRECTORIES: ${Kokkos_INTERFACE_INCLUDE_DIRECTORIES}")

    if(QUICC_USE_KOKKOS_CUDA)
        if(NOT DEFINED Kokkos_ENABLE_CUDA OR NOT ${Kokkos_ENABLE_CUDA})
            message(
                FATAL_ERROR
                    "Enable Kokkos Cuda or unset QUICC_HAS_CUDA_BACKEND to continue with OpenMP!"
            )
        endif()
        message(VERBOSE "Kokkos CUDA Enabled = ${Kokkos_ENABLE_CUDA}")
        kokkos_check(OPTIONS CUDA_LAMBDA)

        # get cuda flags from the wrapper
        # alternatively we can strip Kokkos_INTERFACE_COMPILE_OPTIONS
        # when defined
        execute_process(
            COMMAND ${Kokkos_CXX_COMPILER} --show
            OUTPUT_VARIABLE _wrapper_command
            ERROR_QUIET)
        string(REGEX REPLACE [[\n\v\c\c]] "" _wrapper_flags ${_wrapper_command})
        string(STRIP "${_wrapper_flags}" _wrapper_flags)
        message(DEBUG "_wrapper_flags ${_wrapper_flags}")

        # this could be done per target if we need to compile other parts of QuICC
        # with different CUDA settings
        set(CMAKE_CUDA_FLAGS "${_wrapper_flags} ${_openmp}")

        # strip target
        include(StripTarget)
        strip_target(${_KK_TARGET} LANGUAGE CUDA)


    else()
        string(FIND "${CMAKE_CXX_FLAGS}" "${_openmp}" _pos)
        if(_pos EQUAL -1)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_openmp}")
        endif()

        # strip target
        include(StripTarget)
        if(QUICC_USE_KOKKOS_HIP)
            strip_target(${_KK_TARGET} LANGUAGE HIP)
        else()
            strip_target(${_KK_TARGET} LANGUAGE CXX)
        endif()
    endif()

    get_target_property(Kokkos_INTERFACE_LINK_LIBRARIES ${_KK_TARGET} INTERFACE_LINK_LIBRARIES)
    message(DEBUG "Kokkos_INTERFACE_LINK_LIBRARIES: ${Kokkos_INTERFACE_LINK_LIBRARIES}")

    add_library(QuICC::Kokkos INTERFACE IMPORTED)
    set_target_properties(QuICC::Kokkos
        PROPERTIES
	    INTERFACE_LINK_LIBRARIES "${Kokkos_INTERFACE_LINK_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${Kokkos_INTERFACE_INCLUDE_DIRECTORIES}"
	    INTERFACE_COMPILE_OPTIONS "${_openmp}"
    )

    # done with setting up Kokkos target
    unset(_KK_TARGET)

    # defs to disable code sections
    target_compile_definitions(QuICC::Kokkos INTERFACE QUICC_USE_KOKKOS)
    if(QUICC_USE_KOKKOS_CUDA)
        target_compile_definitions(QuICC::Kokkos INTERFACE QUICC_USE_KOKKOS_CUDA)
    endif()
    if(QUICC_USE_KOKKOS_HIP)
        target_compile_definitions(QuICC::Kokkos INTERFACE QUICC_USE_KOKKOS_HIP)
    endif()

    list(POP_BACK CMAKE_MESSAGE_INDENT)

    #
    # to be separated
    #

    message(STATUS "Setup Kokkos Kernels")
    list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

    # Kokkos Kernels
    if(QUICC_USE_KOKKOS_KERNELS)
        find_package(
            KokkosKernels
            HINTS
            ${KOKKOS_DIR}
            ${KOKKOS_DIR}/lib/CMake/KokkosKernels
            ${KOKKOS_DIR}/lib64/cmake/KokkosKernels
            ${TRILINOS_DIR}
            ${TRILINOS_DIR}/lib/cmake/KokkosKernels
            ${TRILINOS_DIR}/lib64/cmake/KokkosKernels
            REQUIRED)

        message(VERBOSE "Found Kokkos Kernels")
        message(VERBOSE "KokkosKernels_LIBRARIES = ${KokkosKernels_LIBRARIES}")
        message(VERBOSE "KokkosKernels_LIBRARY_DIRS = ${KokkosKernels_LIBRARY_DIRS}")

        # _KKK_TARGET is set as a local variable
        # do not use outside this file
        set(_KKK_TARGET "Kokkos::KokkosKernels")
        if(NOT TARGET ${_KKK_TARGET})
            message(DEBUG "Kokkos kernel target is not defined")
            add_library(${_KKK_TARGET} INTERFACE IMPORTED)
            set_property(
                TARGET ${_KKK_TARGET}
                PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${KokkosKernels_INCLUDE_DIRS} ${KokkosKernels_TPL_INCLUDE_DIRS})
            set_property(
                TARGET ${_KKK_TARGET}
                PROPERTY INTERFACE_LINK_LIBRARIES ${KokkosKernels_LIBRARIES} ${KokkosKernels_TPL_LIBRARIES})
            set_property(
                TARGET ${_KKK_TARGET}
                PROPERTY INTERFACE_LINK_DIRECTORIES ${KokkosKernels_LIBRARY_DIRS})
        endif()

        # done with setting up Kokkos Kernels target
        unset(_KKK_TARGET)
    endif()

    list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()
