#
# Utility to add library with modern target and export
#
# target
#     target library
# PARENT
#     parent library, needed only if current library is a component
# INCLUDE
#     list of include paths, requires PUBLIC or PRIVATE attribute
# SRC
#     list of sources
# LINK
#     list of libraries to link paths, requires PUBLIC or PRIVATE attribute
# CLASSICHEADER
#     Add QuICC to install header path


# support function to split public and private lists
function(_split_list _list _list_public _list_private)
    message(DEBUG "_list: ${_list}")
    list(FIND _list "PRIVATE" _pripos)
    list(FIND _list "PUBLIC" _pubpos)

    if(${_pubpos} GREATER_EQUAL 0)
      if(${_pripos} GREATER ${_pubpos})
        list(SUBLIST _list ${_pubpos} ${_pripos} _list)
      else()
        list(SUBLIST _list ${_pubpos} -1 _list)
      endif()
    endif()
    list(POP_FRONT _list _list)
    message(DEBUG "PUBLIC list: ${_list}")
    set(${_list_public} "${_list}" PARENT_SCOPE)

    if(${_pripos} GREATER_EQUAL 0)
      if(${_pripos} LESS ${_pubpos})
        list(SUBLIST _list ${_pripos} ${_pubpos} _list)
      else()
        list(SUBLIST _list ${_pripos} -1 _list)
      endif()
    endif()
    list(POP_FRONT _list _list)
    message(DEBUG "PRIVATE list: ${_list}")
    set(${_list_private} "${_list}" PARENT_SCOPE)
endfunction()


function(quicc_add_library target)
  # parse inputs
  set(options CLASSICHEADER KOKKOSLIB)
  set(oneValueArgs PARENT)
  set(multiValueArgs INCLUDE SRC LINK)
  cmake_parse_arguments(QAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_library")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QAL_PARENT: ${QAL_PARENT}")
  message(DEBUG "QAL_INCLUDE: ${QAL_INCLUDE}")
  message(DEBUG "QAL_SRC: ${QAL_SRC}")
  message(DEBUG "QAL_CLASSICHEADER: ${QAL_CLASSICHEADER}")

  # Create component library if there is a parent
  if(QAL_PARENT)
    set(QUICC_CURRENT_SUBCOMPONENT_LIB ${target})
    set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "${QAL_PARENT}::" )
    set(_lib ${QAL_PARENT}_${QUICC_CURRENT_SUBCOMPONENT_LIB})
  else()
    set(QUICC_CURRENT_COMPONENT_LIB_NAMESPACE "" )
    set(_lib ${target})
  endif()

  if(QAL_KOKKOSLIB AND QUICC_USE_KOKKOS)
    include(setup/KokkosAddLibrary)
    quicc_add_library_kokkos(${_lib} ${QAL_SRC})
  else()
    add_library(${_lib} ${QAL_SRC})
  endif()
  set_target_properties(${_lib} PROPERTIES LINKER_LANGUAGE CXX)

  if(QAL_INCLUDE)
    # parse private/public
    _split_list("${QAL_INCLUDE}" _include_pub _include_pri)

    if(_include_pub)
      target_include_directories(${_lib} PUBLIC
          "$<BUILD_INTERFACE:${_include_pub}>"
          "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
    endif()

    if(_include_pri)
      target_include_directories(${_lib} PRIVATE
          "$<BUILD_INTERFACE:${_include_pri}>"
          "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
    endif()
  endif()

  if(QAL_LINK)
     _split_list("${QAL_LINK}" _link_pub _link_pri)
    target_link_libraries(${_lib} PUBLIC ${_link_pub})
    target_link_libraries(${_lib} PRIVATE ${_link_pri})
  endif()

  # Alias
  set(_alias ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${target})
  message(DEBUG "_alias: ${_alias}")
  add_library(${_alias} ALIAS ${_lib})

  # Export info
  if(QAL_CLASSICHEADER)
      list(TRANSFORM _include_pub APPEND "/QuICC")
  endif()

  quicc_export_target(${_lib}
    COMPONENT ${QUICC_CURRENT_COMPONENT_LIB_NAMESPACE}${target}
    DIRECTORIES
      ${_include_pub}
    FILES_MATCHING_PATTERN "*.h*"
      )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
