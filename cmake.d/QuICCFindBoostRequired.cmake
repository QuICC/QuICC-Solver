#
# Convenience wrapper to findBoost
#

macro(quicc_find_boost_required)
    if(NOT TARGET Boost::headers)
        if(NOT BOOST_ROOT)
            message(VERBOSE "setting BOOST_ROOT")
            list(APPEND _ALL_PATHS $ENV{CPATH} $ENV{C_INCLUDE_PATH} $ENV{CPLUS_INCLUDE_PATH})
            if(NOT "${_ALL_PATHS}" STREQUAL "")
                string(REPLACE ":" ";" _ALL_PATHS "${_ALL_PATHS}")
                include(ListFindRegex)
                quicc_list(FIND_REGEX _ALL_PATHS "boost" BOOST_ROOT)
                message(VERBOSE "BOOST_ROOT: ${BOOST_ROOT}")
            endif()
        endif()
        set(_version "1.78")
        find_package(Boost ${_version})
        if(NOT Boost_FOUND)
            message(FATAL_ERROR "Could not find Boost, required >= ${_version}, try to specify path: -DBOOST_ROOT=</path/to/boost>")
        endif()
        # Imported target does not necessarily have global scope
        set_target_properties(Boost::headers PROPERTIES IMPORTED_GLOBAL TRUE)
    endif()
endmacro()
