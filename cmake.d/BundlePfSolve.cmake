#
# Get system PfSolve if available otherwise fetch it
#
message(DEBUG "PfSolve")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_PFSOLVE_VERSION "f93789ca8e47faea5cf263122244c9b2b247eac0")

option(QUICC_USE_SYSTEM_PFSOLVE "Use system installed PfSolve." OFF)
option(QUICC_USE_QUADDOUBLEDOUBLE_PFSOLVE "Use quad double-double version of PfSolve." OFF)

add_library(PfSolve INTERFACE IMPORTED)
set(QUICC_MESSAGE_QUIET ON)
set(VKFFT_BACKEND 1 CACHE STRING "0 - Vulkan, 1 - CUDA, 2 - HIP, 3 - OpenCL")
add_definitions(-DVKFFT_BACKEND=${VKFFT_BACKEND})
add_definitions(-DQUICC_USE_PFSOLVE)
if(QUICC_USE_QUADDOUBLEDOUBLE_PFSOLVE)
    add_definitions(-DPFSOLVE_FP128)
endif() 

if(QUICC_USE_SYSTEM_PFSOLVE)
    find_path(PFSOLVE_DIR
    NAMES pfSolve.h
    HINTS "${PFSOLVE_SOURCE_DIR}"
    DOC "pfSolve directory"
    )
    if(NOT PFSOLVE_DIR)
    message(FATAL_ERROR "Not found PfSolve")
    endif()
    if(NOT PFSOLVE_SOURCE_DIR)
        set(PFSOLVE_SOURCE_DIR "${PFSOLVE_DIR}")
    endif()
    if(NOT PFSOLVE_KERNELS_DIR)
        set(PFSOLVE_KERNELS_DIR "${PFSOLVE_SOURCE_DIR}")
    endif()
    add_definitions(-DPFSOLVE_KERNELS_DIR="${PFSOLVE_KERNELS_DIR}")
else()
    include(FetchContent)
    FetchContent_Declare(
        PfSolve
        GIT_REPOSITORY "https://github.com/QuICC/PfSolve"
        GIT_TAG ${QUICC_PFSOLVE_VERSION}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    FetchContent_Populate(PfSolve)
    set(PFSOLVE_SOURCE_DIR "${FETCHCONTENT_BASE_DIR}/pfsolve-src/PfSolve")
    set(PFSOLVE_KERNELS_DIR "${FETCHCONTENT_BASE_DIR}/pfsolve-build")
    add_definitions(-DPFSOLVE_KERNELS_DIR="${PFSOLVE_KERNELS_DIR}")
endif()
    
if(${VKFFT_BACKEND} EQUAL 0)
elseif(${VKFFT_BACKEND} EQUAL 1)
    enable_language(CUDA)
    add_definitions(-DCUDA_TOOLKIT_ROOT_DIR="${CUDA_TOOLKIT_ROOT_DIR}")
elseif(${VKFFT_BACKEND} EQUAL 2)
    find_package(hip)
endif()
if(${VKFFT_BACKEND} EQUAL 0)
elseif(${VKFFT_BACKEND} EQUAL 1)
    target_link_libraries(PfSolve INTERFACE
    cuda
    nvrtc
    )
    target_include_directories(PfSolve INTERFACE
    "${CUDAToolkit_INCLUDE_DIRS}"
    )
elseif(${VKFFT_BACKEND} EQUAL 2)
    find_library(HIP_HIPRTC_LIB libhiprtc hiprtc HINTS "${HIP_LIB_INSTALL_DIR}/lib64" "${LIBHIPRTC_LIBRARY_DIR}" "${HIP_LIB_INSTALL_DIR}/lib/x64" /usr/lib64 /usr/local/cuda/lib64)
    target_link_libraries(PfSolve INTERFACE
    hip::host
    ${HIP_HIPRTC_LIB}
    )
    target_include_directories(PfSolve INTERFACE
    "${hip_INCLUDE_DIR}"
    )
endif()
if(QUICC_USE_QUADDOUBLEDOUBLE_PFSOLVE)
    target_link_libraries(PfSolve INTERFACE
    quadmath
    )
endif()

if((${VKFFT_BACKEND} EQUAL 1) OR (${VKFFT_BACKEND} EQUAL 2))
    target_include_directories(PfSolve INTERFACE "${PFSOLVE_SOURCE_DIR}")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
