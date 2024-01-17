#
# Get system HighFive if available otherwise fetch it
#
message(DEBUG "HighFive")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_HIGHFIVE_VERSION "2.8.0")

option(QUICC_USE_SYSTEM_HIGHFIVE "Use system installed HighFive." OFF)

if(QUICC_USE_SYSTEM_HIGHFIVE)
    find_package(HighFive ${QUICC_HIGHFIVE_VERSION} REQUIRED)
    if(NOT HighFive_FOUND)
        message(SEND_ERROR "To use bundled HighFive add: -DQUICC_USE_SYSTEM_HIGHFIVE=OFF or set HighFive_DIR.")
    endif()
else()
    if(NOT TARGET HighFive)
        set(QUICC_MESSAGE_QUIET ON)

        include(FetchContent)
            FetchContent_Declare(
                h5
                GIT_REPOSITORY https://github.com/BlueBrain/HighFive
                GIT_TAG "v${QUICC_HIGHFIVE_VERSION}"
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
            )

        # set some vars
        if(QUICC_MPI)
            set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "")
        endif()
        # requires system and serialization which are not header only
        set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "")

        set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "")
        set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "")
        set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "")

        FetchContent_MakeAvailable(h5)

        # hide HighFive vars
        include(MarkAsAdvancedAll)
        mark_as_advanced_all(HIGHFIVE)

        unset(QUICC_MESSAGE_QUIET)
    endif()
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
