#
# Utility to generate required classes for tags
#
# NAMESPACE
#     namespace of the tags relative to QuICC
# BASECLASS
#     name of base class of tags
# COMMON_DIR
#     relative path to Common component root
# EXCLUDED
#     files excluded for glob
# VALUE
#     tags also carries a numerical value
#
function(quicc_register_tags)
  # parse inputs
  set(oneValueArgs NAMESPACE BASECLASS COMMON_DIR VALUE)
  set(multiValueArgs EXCLUDED)
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_register_tags")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QAT_NAMESPACE: ${QAT_NAMESPACE}")
  message(DEBUG "QAT_BASECLASS: ${QAT_BASECLASS}")
  message(DEBUG "QAT_COMMON_DIR: ${QAT_COMMON_DIR}")
  if(QAT_VALUE)
    message(DEBUG "${QAT_VALUE}")
  endif()

  if(NOT QAT_COMMON_DIR)
    set(QAT_COMMON_DIR ".")
  endif()

  if(NOT QAT_COMMON_DIR)
    set(QAT_EXCLUDED "")
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

  # Configure namespace path and guard namespace
  set(idNS "${QAT_NAMESPACE}")
  string(TOUPPER "${idNS}" _tmp)
  string(REPLACE "/" "_" IDNS "${_tmp}")
  message(DEBUG "idNS, IDNS: ${idNS}, ${IDNS}")

  # set baseclase name and upper case version
  set(idID "${QAT_BASECLASS}")
  string(TOUPPER ${QAT_BASECLASS} IDID)
  message(DEBUG "idID, IDID: ${idID}, ${IDID}")

  # set baseclase name and upper case version
  if(QAT_VALUE)
    set(idARG "${QAT_VALUE}")
    set(_idSIG "const MHDFloat ${idARG}")
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
  file(GLOB AllFiles include/QuICC/${QAT_NAMESPACE}/*.hpp)

  # generate header and id list
  set(_idList "")
  set(_idHeader "")
  foreach(fullname ${AllFiles})
    get_filename_component(fname ${fullname} NAME)
    if(NOT fname IN_LIST QAT_EXCLUDED)
      set(_idHeader "${_idHeader}\n#include \"QuICC/${QAT_NAMESPACE}/${fname}\"")
      get_filename_component(cname ${fullname} NAME_WE)
      set(_idList "${_idList}\n      ${cname}::id();")
    endif()
  endforeach(fullname)

  # Configure registration file
  configure_file(
    "${QAT_COMMON_DIR}/include/QuICC/IdTools/registerAll.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/registerAll.hpp"
    )

  # Configure ID interface registration files
  if(QAT_VALUE)
    configure_file(
      "${QAT_COMMON_DIR}/include/QuICC/IdTools/IValuedId.hpp.in"
      "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/${QAT_BASECLASS}.hpp"
      )
  else()
    configure_file(
      "${QAT_COMMON_DIR}/include/QuICC/IdTools/IId.hpp.in"
      "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/${QAT_BASECLASS}.hpp"
      )
  endif()

  configure_file(
    "${QAT_COMMON_DIR}/include/QuICC/IdTools/IRegisterId.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/IRegisterId.hpp"
    )
  configure_file(
    "${QAT_COMMON_DIR}/include/QuICC/IdTools/ICreator.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/ICreator.hpp"
    )
  configure_file(
    "${QAT_COMMON_DIR}/include/QuICC/IdTools/CreatorImpl.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/CreatorImpl.hpp"
    )
  configure_file(
    "${QAT_COMMON_DIR}/include/QuICC/IdTools/Coordinator.hpp.in"
    "${PROJECT_BINARY_DIR}/include/QuICC/${QAT_NAMESPACE}/Coordinator.hpp"
    )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
