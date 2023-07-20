#
# Utility to generate required classes for tags
#
# NAMESPACE
#     namespace of the tags relative to QuICC
# BASECLASS
#     name of base class of tags
# COMMON_DIR
#     relative path to Common component root
# TAGS
#     tag class names, tag is created from class name: CamelCase -> camel_case
# EXCLUDED
#     files excluded for glob
# PREFIX
#     add prefix to tag
# VALUE
#     tags also carries a numerical value
# REGISTRATOR
#     name of convenience header to register all tags. registerAll by default
#
function(quicc_register_tags)
  # parse inputs
  set(oneValueArgs NAMESPACE BASECLASS COMMON_DIR PREFIX VALUE REGISTRATOR)
  set(multiValueArgs EXCLUDED TAGS)
  cmake_parse_arguments(QRT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_register_tags")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QRT_NAMESPACE: ${QRT_NAMESPACE}")
  message(DEBUG "QRT_BASECLASS: ${QRT_BASECLASS}")
  message(DEBUG "QRT_COMMON_DIR: ${QRT_COMMON_DIR}")
  if(QRT_VALUE)
    message(DEBUG "${QRT_VALUE}")
  endif()
  if(QRT_PREFIX)
    message(DEBUG "${QRT_PREFIX}")
  endif()

  if(NOT QRT_REGISTRATOR)
    set(QRT_REGISTRATOR "registerAll")
  endif()

  if(NOT QRT_COMMON_DIR)
    set(QRT_COMMON_DIR ".")
  endif()

  if(NOT QRT_COMMON_DIR)
    set(QRT_EXCLUDED "")
  endif()

  # File configuration
  # @idNS@ Namespace path
  # @IDNS@ Namespace path for header guard (/ -> _)
  # @cxxNS@ C++ namespace
  # @cxxNS_@ C++ namespace closing brackets
  # @idID@ Base class name
  # @IDID@ Base class name in upper case
  # @idARG@ Argument name for value
  # @_idCOMMA@ Introduce comma in function arguments
  # @_idNAN@   Value for bad value
  # @_idSIG@   Argument signature
  # @_idARGDefault@ Argument default value

  # Configure namespace path and guard namespace
  set(idNS "${QRT_NAMESPACE}")
  string(TOUPPER "${idNS}" _tmp)
  string(REPLACE "/" "_" IDNS "${_tmp}")
  message(DEBUG "idNS, IDNS: ${idNS}, ${IDNS}")

  # set baseclase name and upper case version
  set(idID "${QRT_BASECLASS}")
  string(TOUPPER ${QRT_BASECLASS} IDID)
  message(DEBUG "idID, IDID: ${idID}, ${IDID}")

  # set baseclase name and upper case version
  if(QRT_VALUE)
    set(idARG "${QRT_VALUE}")
    set(_idSIG "const MHDFloat ${idARG}")
    set(_idARGDefault " = std::numeric_limits<MHDFloat>::signaling_NaN()")
    set(_idCOMMA ", ")
    set(_idNAN "-424242.424242")
  endif()

  # Configure C++ namespace and closing brackets
  set(_tmp "namespace ${idNS} {")
  string(REPLACE "/" " { \n\nnamespace " cxxNS "${_tmp}")
  # ... extract / from path
  string(REGEX REPLACE "[^ /]" "" _tmp "${idNS}")
  set(_tmp "}${_tmp}")
  string(REPLACE "/" "\n}" cxxNS_ "${_tmp}")

  # glob all ID files
  if(QRT_TAGS)
    foreach(idCLASS ${QRT_TAGS})
      string(TOUPPER ${idCLASS} IDCLASS)

      # Create tag from class name
      string(REGEX REPLACE "([A-Z][a-z0-9]?)" "_\\1" _idTAG "${idCLASS}")
      string(TOLOWER ${_idTAG} _idTAG)
      string(REGEX REPLACE "__" "_" _idTAG "${_idTAG}")
      string(REGEX REPLACE "^_" "" _idTAG "${_idTAG}")
      set(idTAG "${QRT_PREFIX}${_idTAG}")

      # Configure tag file
      configure_file(
        "${QRT_COMMON_DIR}/include/QuICC/IdTools/IdClass.hpp.in"
        "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/${idCLASS}.hpp"
        )
      list(APPEND AllFiles "${idCLASS}.hpp")

      unset(idCLASS)
      unset(idTAG)
      unset(IDCLASS)
    endforeach()
  else()
    file(GLOB AllFiles include/QuICC/${QRT_NAMESPACE}/*.hpp)
  endif()

  # generate header and id list
  set(_idList "")
  set(_idHeader "")
  foreach(fullname ${AllFiles})
    get_filename_component(fname ${fullname} NAME)
    if(NOT fname IN_LIST QRT_EXCLUDED)
      set(_idHeader "${_idHeader}\n#include \"QuICC/${QRT_NAMESPACE}/${fname}\"")
      get_filename_component(cname ${fullname} NAME_WE)
      set(_idList "${_idList}\n      ${cname}::id();")
    endif()
  endforeach(fullname)

  set(_registrator ${QRT_REGISTRATOR})
  string(TOUPPER "${_registrator}" _REGISTRATOR)
  configure_file(
    "${QRT_COMMON_DIR}/include/QuICC/IdTools/registerAll.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/${QRT_REGISTRATOR}.hpp"
    )

  # Configure ID interface registration files
  if(QRT_VALUE)
    configure_file(
      "${QRT_COMMON_DIR}/include/QuICC/IdTools/IValuedId.hpp.in"
      "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/${QRT_BASECLASS}.hpp"
      )
  else()
    configure_file(
      "${QRT_COMMON_DIR}/include/QuICC/IdTools/IId.hpp.in"
      "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/${QRT_BASECLASS}.hpp"
      )
  endif()

  configure_file(
    "${QRT_COMMON_DIR}/include/QuICC/IdTools/IRegisterId.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/IRegisterId.hpp"
    )
  configure_file(
    "${QRT_COMMON_DIR}/include/QuICC/IdTools/ICreator.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/ICreator.hpp"
    )
  configure_file(
    "${QRT_COMMON_DIR}/include/QuICC/IdTools/CreatorImpl.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/CreatorImpl.hpp"
    )
  configure_file(
    "${QRT_COMMON_DIR}/include/QuICC/IdTools/Coordinator.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QRT_NAMESPACE}/Coordinator.hpp"
    )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
