#
# Get system VkFFT if available otherwise fetch it
#
message(DEBUG "VkFFT")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_VKFFT_VERSION "v1.3.2")

option(QUICC_USE_SYSTEM_VKFFT "Use system installed VkFFT." OFF)
option(QUICC_USE_DOUBLEDOUBLE_VKFFT "Use double-double version of VkFFT." OFF)

add_library(VkFFT INTERFACE IMPORTED)
set(QUICC_MESSAGE_QUIET ON)
set(VKFFT_BACKEND 1 CACHE STRING "0 - Vulkan, 1 - CUDA, 2 - HIP, 3 - OpenCL")
add_definitions(-DVKFFT_BACKEND=${VKFFT_BACKEND})
if(QUICC_USE_DOUBLEDOUBLE_VKFFT)
    add_definitions(-DVKFFT_USE_DOUBLEDOUBLE_FP128)
endif() 

if(QUICC_USE_SYSTEM_VKFFT)
    find_path(VKFFT_DIR
    NAMES vkFFT.h
    HINTS "${VKFFT_SOURCE_DIR}"
    DOC "vkFFT directory"
    )

    if(NOT VKFFT_DIR)
    message(FATAL_ERROR "Not found VkFFT")
    endif()
    if(NOT VKFFT_SOURCE_DIR)
        set(VKFFT_SOURCE_DIR "${VKFFT_DIR}")
    endif()
    if(NOT VKFFT_KERNELS_DIR)
        set(VKFFT_KERNELS_DIR="${VKFFT_SOURCE_DIR}")
    endif()
    add_definitions(-DVKFFT_KERNELS_DIR="${VKFFT_KERNELS_DIR}")
else()
    include(FetchContent)
    FetchContent_Declare(
        VkFFT
        GIT_REPOSITORY https://github.com/DTolm/VkFFT
        GIT_TAG ${QUICC_VKFFT_VERSION}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    FetchContent_Populate(VkFFT)
    set(VKFFT_SOURCE_DIR "${FETCHCONTENT_BASE_DIR}/vkfft-src/vkFFT")
    set(VKFFT_KERNELS_DIR="${FETCHCONTENT_BASE_DIR}/vkfft-build")
    add_definitions(-DVKFFT_KERNELS_DIR="${VKFFT_KERNELS_DIR}")
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
    target_link_libraries(VkFFT INTERFACE
    cuda
    nvrtc
    )
    target_include_directories(VkFFT INTERFACE
    "${CUDAToolkit_INCLUDE_DIRS}"
    )
elseif(${VKFFT_BACKEND} EQUAL 2)
    find_library(HIP_HIPRTC_LIB libhiprtc hiprtc HINTS "${HIP_LIB_INSTALL_DIR}/lib64" "${LIBHIPRTC_LIBRARY_DIR}" "${HIP_LIB_INSTALL_DIR}/lib/x64" /usr/lib64 /usr/local/cuda/lib64)
    target_link_libraries(VkFFT INTERFACE
    hip::host
    ${HIP_HIPRTC_LIB}
    )
    target_include_directories(VkFFT INTERFACE
    "${hip_INCLUDE_DIR}"
    )
endif()
if(QUICC_USE_DOUBLEDOUBLE_VKFFT)
    target_link_libraries(VkFFT INTERFACE
    quadmath
    )
endif()
if((${VKFFT_BACKEND} EQUAL 1) OR (${VKFFT_BACKEND} EQUAL 2))
    target_include_directories(VkFFT INTERFACE "${VKFFT_SOURCE_DIR}")
    #target_include_directories(VkFFT INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/QuICC/Math/PfSolve>")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
