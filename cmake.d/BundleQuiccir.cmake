#
# Get system Quiccir if available otherwise fetch it
#
message(DEBUG "Quiccir")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_QUICCIR_VERSION "integrate")

option(QUICC_USE_SYSTEM_QUICCIR "Use system installed Quiccir." OFF)

if(QUICC_USE_SYSTEM_QUICCIR)
    find_package(Quiccir ${QUICC_QUICCIR_VERSION} REQUIRED)
    # find_package(Quiccir ${QUICC_QUICCIR_VERSION} REQUIRED)
    if(NOT Quiccir_FOUND)
        message(SEND_ERROR "To use bundled Quiccir add: -DQUICC_USE_SYSTEM_QUICCIR=OFF or set Quiccir_DIR.")
    endif()
else()
    if(NOT TARGET Quiccir)
        set(QUICC_MESSAGE_QUIET ON)

        include(FetchContent)
            FetchContent_Declare(
                quiccir
                GIT_REPOSITORY https://github.com/QuICC/quiccir
                GIT_TAG "${QUICC_QUICCIR_VERSION}"
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
            )

        # set some vars
        # if(QUICC_MPI)
        #     set(Quiccir_PARALLEL_HDF5 ON CACHE BOOL "")
        # endif()
        # # requires system and serialization which are not header only
        # set(Quiccir_USE_BOOST OFF CACHE BOOL "")

        # set(Quiccir_UNIT_TESTS OFF CACHE BOOL "")
        # set(Quiccir_EXAMPLES OFF CACHE BOOL "")
        # set(Quiccir_BUILD_DOCS OFF CACHE BOOL "")

        FetchContent_MakeAvailable(quiccir)

        # hide LLVM vars
        include(MarkAsAdvancedAll)
        mark_as_advanced_all(LLVM)

        unset(QUICC_MESSAGE_QUIET)
    endif()
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
