#
# Utility to add tests
#
# target
#     name of the test and [target].cpp is the name of the source
# COMMAND
#     command passed to add_test plus optionally [--types=<TYPES>]
# KEYWORD
#     variable for file configuration
# TYPES
#     list of optional types for test executables (see COMMAND)
# IDS
#     list of id option
# ULP
#     set single ulp option
# ULPS
#     list of ulps, should match IDS size
# DISABLED
#     if set to true, set test as disabled
#
function(quicc_add_test target)
  # parse inputs
  set(oneValueArgs COMMAND KEYWORD ULP DISABLED)
  set(multiValueArgs TYPES IDS ULPS)
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_test")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAT_COMMAND: ${QAT_COMMAND}")

  string(REGEX REPLACE "Test" "" ${QAT_KEYWORD} ${target})
  message(DEBUG "QAT_KEYWORD: ${QAT_KEYWORD} : ${${QAT_KEYWORD}}")
  message(DEBUG "QAT_ULP: ${QAT_ULP}")
  message(DEBUG "QAT_ULPS: ${QAT_ULPS}")
  message(DEBUG "QAT_DISABLED: ${QAT_DISABLED}")

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp")
    message(VERBOSE "Adding ${target}")
    set(_keyword ${${QAT_KEYWORD}})
    set(testfile ${target}.cpp)
  else()
    message(VERBOSE "Generating ${target} from template")
    string(REGEX REPLACE ":+" ";" _tlist ${${QAT_KEYWORD}})
    string(REGEX REPLACE ":+" ";" _klist ${QAT_KEYWORD})
    if(${CMAKE_VERSION} VERSION_LESS 3.17)
      # patch branch
      foreach(_var_0 IN LISTS _tlist)
        message(DEBUG "_var_0: ${_var_0}")
        # get position
        list(FIND _tlist ${_var_0} _pos)
        message(DEBUG "_pos: ${_pos}")
        if(${_pos} GREATER "-1")
          # get corresponding item in second list
          list(GET _klist ${_pos} _var_1)
          set(${_var_1} ${_var_0})
          message(DEBUG "_var_1: ${_var_1}: ${${_var_1}}")
        endif()
      endforeach()
    else()
      foreach(_var IN ZIP_LISTS _tlist _klist)
        set(${_var_1} ${_var_0})
        message(DEBUG "_var_0: ${_var_0}")
        message(DEBUG "_var_1: ${_var_1}: ${${_var_1}}")
      endforeach()
    endif()

    set(CatchTestName ${${QAT_KEYWORD}})
    message(DEBUG "CatchTestName: ${CatchTestName}")
    string(REGEX REPLACE "::" "/" _cppname ${CatchTestName})
    string(REGEX REPLACE ":" "_" _cppname ${_cppname})
    message(DEBUG "_cppname: ${_cppname}")

    string(REGEX REPLACE ":+" "_" _keyword ${CatchTestName})
    message(DEBUG "_keyword: ${_keyword}")
    set(testfile
      "${CMAKE_CURRENT_BINARY_DIR}/${_cppname}Test.cpp")
    message(DEBUG "testfile: ${testfile}")
    configure_file(
      "TemplateTest.cpp.in"
      "${testfile}"
    )
  endif()

  target_sources(${QAT_COMMAND} PRIVATE
    ${testfile}
  )

  # check for variants
  if(NOT ${QAT_ULP} STREQUAL "")
    set(_ulp_cmd "--ulp ${QAT_ULP}")
    set(_ulp_name "_ulp${QAT_ULP}")
  endif()

  list(LENGTH QAT_TYPES _types_len)
  message(DEBUG "_types_len: ${_types_len}")
  list(LENGTH QAT_IDS _ids_len)
  message(DEBUG "_ids_len: ${_ids_len}")
  if(${_types_len} LESS 1 AND ${_ids_len} LESS 1)
    set(_testname "${QAT_COMMAND}_${_keyword}${_ulp_name}" )
    message(DEBUG "_testname: ${_testname}")
    add_test(
      NAME ${_testname}
      COMMAND ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests ${_ulp_cmd}
      WORKING_DIRECTORY
      "${QUICC_WORK_DIR}"
    )
    if(QAT_DISABLED)
      set_tests_properties(${_testname} PROPERTIES DISABLED ${QAT_DISABLED})
    endif()
  elseif(${_ids_len} LESS 1)
    foreach(_type ${QAT_TYPES})
      set(_testname "${QAT_COMMAND}_${_keyword}_type${_type}${_ulp_name}" )
      message(DEBUG "_testname: ${_testname}")
      add_test(
      NAME ${_testname}
      COMMAND ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests --type=${_type} ${_ulp_cmd}
      WORKING_DIRECTORY
      "${QUICC_WORK_DIR}"
      )
      if(QAT_DISABLED)
        set_tests_properties(${_testname} PROPERTIES DISABLED ${QAT_DISABLED})
      endif()
    endforeach()
  elseif(${_types_len} LESS 1)
    list(LENGTH QAT_ULPS _ulps_len)
    foreach(_id ${QAT_IDS})
      # check for custom ulp
      if(_ulps_len GREATER 1)
        if(NOT _ulps_len EQUAL _ids_len)
          message(SEND_ERROR "ULPS and IDS must have the same")
        endif()
        # get current index
        list(FIND QAT_IDS ${_id} _index)
        list(GET QAT_ULPS ${_index} _ulp)
        set(_ulp_cmd "--ulp ${_ulp}")
        set(_ulp_name "_ulp${_ulp}")
      endif()
      set(_testname "${QAT_COMMAND}_${_keyword}_id${_id}${_ulp_name}" )
      message(DEBUG "_testname: ${_testname}")
      add_test(
      NAME ${_testname}
      COMMAND ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests --id=${_id} ${_ulp_cmd}
      WORKING_DIRECTORY
      "${QUICC_WORK_DIR}"
      )
      if(QAT_DISABLED)
        set_tests_properties(${_testname} PROPERTIES DISABLED ${QAT_DISABLED})
      endif()
    endforeach()
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
