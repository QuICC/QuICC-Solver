#
# Get system Catch2 if available otherwise fetch it
#
message(DEBUG "Catch2")
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

set(QUICC_CATCH2_VERSION "2.13.7")

option(QUICC_USE_SYSTEM_CATCH2 "Use system installed Catch2." OFF)

if(QUICC_USE_SYSTEM_CATCH2)
    find_package(Catch2 ${QUICC_CATCH2_VERSION})
    if(NOT Catch2_FOUND)
        message(SEND_ERROR "To use bundled Catch2 add: -DQUICC_USE_SYSTEM_Catch=OFF or set Catch2_DIR.")
    endif()
else()
    if(NOT TARGET Catch2::Catch2)
        include(FetchContent)
            FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG "v${QUICC_CATCH2_VERSION}"
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
            )

        # save CMAKE_POLICY_DEFAULT_CMP0077
        set(_CMAKE_POLICY_DEFAULT_CMP0077 ${CMAKE_POLICY_DEFAULT_CMP0077})
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

        set(CATCH_INSTALL_DOCS OFF)
        set(CATCH_INSTALL_HELPERS OFF)
        set(CATCH_BUILD_TESTING OFF)

        FetchContent_MakeAvailable(Catch2)

        # hide Catch vars
        include(MarkAsAdvancedAll)
        mark_as_advanced_all(CATCH)

        # restore CMAKE_POLICY_DEFAULT_CMP0077
        set(CMAKE_POLICY_DEFAULT_CMP0077 ${_CMAKE_POLICY_DEFAULT_CMP0077})

    endif()
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
