#
# Get system Eigen if available otherwise fetch it
#
message(DEBUG "Eigen")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_EIGEN_VERSION "3.4.0")

option(QUICC_USE_SYSTEM_EIGEN "Use system installed Eigen." OFF)

if(QUICC_USE_SYSTEM_EIGEN)
    find_package(Eigen3 ${QUICC_EIGEN_VERSION} NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(SEND_ERROR "To use bundled Eigen add: -DQUICC_USE_SYSTEM_EIGEN=OFF or set Eigen3_DIR.")
    endif()
else()
    if(NOT TARGET Eigen3::Eigen)
        include(MessageQuICC)
        set(QUICC_MESSAGE_QUIET ON)

        include(FetchContent)
            FetchContent_Declare(
                eigen
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen
                GIT_TAG ${QUICC_EIGEN_VERSION}
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
            )

        # save CMAKE_POLICY_DEFAULT_CMP0077
        set(_CMAKE_POLICY_DEFAULT_CMP0077 ${CMAKE_POLICY_DEFAULT_CMP0077})
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        # save BUILD_TESTING value
        set(_QUICC_BUILD_TESTING ${BUILD_TESTING})
        set(BUILD_TESTING OFF)
        set(EIGEN_BUILD_DOC OFF)
        set(EIGEN_BUILD_PKGCONFIG OFF)

        FetchContent_MakeAvailable(eigen)

        # hide Eigen vars
        include(MarkAsAdvancedAll)
        mark_as_advanced_all(EIGEN)
        mark_as_advanced_all(PASTIX)
        mark_as_advanced_all(QT)
        mark_as_advanced_all(FETCHCONTENT)

        # restore BUILD_TESTING
        set(BUILD_TESTING ${_QUICC_BUILD_TESTING})

        # restore CMAKE_POLICY_DEFAULT_CMP0077
        set(CMAKE_POLICY_DEFAULT_CMP0077 ${_CMAKE_POLICY_DEFAULT_CMP0077})

        unset(QUICC_MESSAGE_QUIET)
    endif()
endif()

# Enable/disable internal Threading
option(QUICC_EIGEN_ENABLE_OPENMP "Enable/disable Eigen internal threading" OFF)
mark_as_advanced(QUICC_EIGEN_ENABLE_OPENMP)

if(NOT ${QUICC_EIGEN_ENABLE_OPENMP})
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
    # retrieve original target, we cannot set properties on aliased targets
    get_target_property(_aliased Eigen3::Eigen ALIASED_TARGET)
    target_compile_options(${_aliased} INTERFACE "-DEIGEN_DONT_PARALLELIZE")
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
else()
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
    # retrieve original target, we cannot set properties on aliased targets
    get_target_property(_aliased Eigen3::Eigen ALIASED_TARGET)
    target_compile_options(${_aliased} INTERFACE "-fopenmp")
    target_link_options(${_aliased} INTERFACE "-fopenmp")
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
endif()

# Enable/disable vectorization
option(QUICC_EIGEN_ENABLE_VECTORIZATION "Enable/disable Eigen internal SIMD vectorization" OFF)
mark_as_advanced(QUICC_EIGEN_ENABLE_VECTORIZATION)

if(NOT ${QUICC_EIGEN_ENABLE_VECTORIZATION})
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
    # retrieve original target, we cannot set properties on aliased targets
    get_target_property(_aliased Eigen3::Eigen ALIASED_TARGET)
    target_compile_options(${_aliased} INTERFACE "-DEIGEN_DONT_VECTORIZE")
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
endif()

# Enable/disable cuda
option(QUICC_EIGEN_ENABLE_CUDA "Enable/disable Eigen internal CUDA" OFF)
mark_as_advanced(QUICC_EIGEN_ENABLE_CUDA)

if(NOT ${QUICC_EIGEN_ENABLE_CUDA})
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
    # retrieve original target, we cannot set properties on aliased targets
    get_target_property(_aliased Eigen3::Eigen ALIASED_TARGET)
    target_compile_options(${_aliased} INTERFACE "-DEIGEN_NO_CUDA")
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
endif()

# Enable/disable MKL
option(QUICC_EIGEN_ENABLE_MKL "Enable/disable Eigen MKL backend" OFF)
mark_as_advanced(QUICC_EIGEN_ENABLE_MKL)

if(${QUICC_EIGEN_ENABLE_MKL})
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")
    # retrieve original target, we cannot set properties on aliased targets
    get_target_property(_aliased Eigen3::Eigen ALIASED_TARGET)
    target_compile_options(${_aliased} INTERFACE "-DEIGEN_USE_MKL_ALL")
    get_target_property(Eigen_INTERFACE_COMPILE_OPTIONS Eigen3::Eigen INTERFACE_COMPILE_OPTIONS)
    message(DEBUG "Eigen_INTERFACE_COMPILE_OPTIONS: ${Eigen_INTERFACE_COMPILE_OPTIONS}")

    # need MKL
    set(BLA_VENDOR Intel10_64lp)
    find_package(BLAS REQUIRED)
    target_link_libraries(${_aliased} INTERFACE BLAS::BLAS)
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
