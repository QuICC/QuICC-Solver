#
# Get system Quiccir if available otherwise fetch it
#
message(DEBUG "BundleQuiccir")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_QUICCIR_VERSION "nl-ops")

option(QUICC_USE_SYSTEM_QUICCIR "Use system installed quiccir." OFF)

if(QUICC_USE_SYSTEM_QUICCIR)
    # find_package(quiccir ${QUICC_QUICCIR_VERSION} REQUIRED)
    find_package(quiccir REQUIRED)
    if(NOT quiccir_FOUND)
        message(SEND_ERROR "To use bundled quiccir add: -DQUICC_USE_SYSTEM_QUICCIR=OFF or set quiccir_DIR.")
    endif()
else()
    if(NOT TARGET quiccir)
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
        set(QUICCIR_BUILD_TEST OFF CACHE BOOL "")
        set(QUICCIR_BUILD_OPT OFF CACHE BOOL "")
        set(QUICCIR_BUILD_MINIAPP OFF CACHE BOOL "")

        FetchContent_MakeAvailable(quiccir)

        # Hide LLVM vars
        include(MarkAsAdvancedAll)
        mark_as_advanced_all(LLVM)
        mark_as_advanced_all(MLIR)

        unset(QUICC_MESSAGE_QUIET)
    endif()
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
