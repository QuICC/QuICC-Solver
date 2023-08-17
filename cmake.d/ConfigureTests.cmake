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
# STEPS
#     number of steps for the profiling test, -1 means no profiling
# SPLITS
#     list of cpus:rank pairs, should match IDS size
# PERFONLY
#     if true add only perf test
# OPTIONS
#     add generic options given as cmd:value converted to --cmd value
#

# support function
function(__add_test _testname)
  # parse inputs
  set(oneValueArgs DIS PRF)
  set(multiValueArgs COMM STP)
  cmake_parse_arguments(_QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "__add_test")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "_QAT_COMM: ${_QAT_COMM}")
  message(DEBUG "_QAT_DIS: ${_QAT_DIS}")
  message(DEBUG "_QAT_STP: ${_QAT_STP}")
  message(DEBUG "_QAT_PRF: ${_QAT_PRF}")

  # performance only test
  if(NOT _QAT_PRF)
    add_test(
      NAME ${_testname}
      COMMAND ${_QAT_COMM}
      WORKING_DIRECTORY
      "${QUICC_WORK_DIR}"
    )
  endif()
  if(_QAT_DIS)
    set_tests_properties(${_testname} PROPERTIES DISABLED ${_QAT_DIS})
  endif()
  # perf test only if number of steps is defined
  if(_QAT_STP)
    if(${_QAT_STP} GREATER "0")
      add_test(
        NAME prof_${_testname}
        COMMAND ${_QAT_COMM} --timeOnly --iter ${_QAT_STP}
        WORKING_DIRECTORY
        "${QUICC_WORK_DIR}"
      )
    endif()
  endif()
endfunction()


# main funtion
function(quicc_add_test target)
  # parse inputs
  set(options PERFONLY)
  set(oneValueArgs COMMAND KEYWORD ULP DISABLED)
  set(multiValueArgs TYPES IDS ULPS STEPS SPLITS OPTIONS)
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_test")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAT_COMMAND: ${QAT_COMMAND}")

  string(REGEX REPLACE "Test" "" ${QAT_KEYWORD} ${target})
  message(DEBUG "QAT_KEYWORD: ${QAT_KEYWORD} : ${${QAT_KEYWORD}}")
  message(DEBUG "QAT_ULP: ${QAT_ULP}")
  message(DEBUG "QAT_ULPS: ${QAT_ULPS}")
  message(DEBUG "QAT_STEPS: ${QAT_STEPS}")
  message(DEBUG "QAT_SPLITS: ${QAT_SPLITS}")
  message(DEBUG "QAT_DISABLED: ${QAT_DISABLED}")
  message(DEBUG "QAT_PERFONLY: ${QAT_PERFONLY}")
  message(DEBUG "QAT_OPTIONS: ${QAT_OPTIONS}")

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
    string(REGEX REPLACE "<" "" _cppname ${_cppname})
    string(REGEX REPLACE ">" "" _cppname ${_cppname})
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

  # Process generic options
  set(_opt_cmd )
  set(_opt_name "")
  foreach(_opt_key ${QAT_OPTIONS})
    string(REGEX REPLACE ":" " " _opt ${_opt_key})
    list(APPEND _opt_cmd "--${_opt}")
    string(REGEX REPLACE ":" "-" _opt ${_opt_key})
    set(_opt_name "${_opt_name}_${_opt}")
  endforeach()
  message(DEBUG "_opt_cmd: ${_opt_cmd}")
  message(DEBUG "_opt_name: ${_opt_name}")

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
    set(_testname "${QAT_COMMAND}_${_keyword}${_ulp_name}${_opt_name}" )
    message(DEBUG "_testname: ${_testname}")
    __add_test(${_testname}
      COMM ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests ${_opt_cmd} ${_ulp_cmd}
      DIS ${QAT_DISABLED}
      STP ${QAT_STEPS}
      PRF ${QAT_PERFONLY}
    )
  elseif(${_ids_len} LESS 1)
    foreach(_type ${QAT_TYPES})
      set(_testname "${QAT_COMMAND}_${_keyword}_type${_type}${_ulp_name}${_opt_name}" )
      message(DEBUG "_testname: ${_testname}")
      __add_test(${_testname}
        COMM ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests ${_opt_cmd} --type=${_type} ${_ulp_cmd}
        DIS ${QAT_DISABLED}
        STP ${QAT_STEPS}
        PRF ${QAT_PERFONLY}
      )
    endforeach()
  elseif(${_types_len} LESS 1)
    list(LENGTH QAT_ULPS _ulps_len)
    list(LENGTH QAT_SPLITS _splits_len)
    # init loop index
    set(_counter "0")
    foreach(_id ${QAT_IDS})
      message(DEBUG "_id: ${_id}")
      # check for custom ulp
      if(_ulps_len GREATER 0)
        if(NOT _ulps_len EQUAL _ids_len)
          message(SEND_ERROR "ULPS and IDS must have the same length")
        endif()
        list(GET QAT_ULPS ${_counter} _ulp)
        list(GET QAT_STEPS ${_counter} _steps)
        set(_ulp_cmd "--ulp ${_ulp}")
        set(_ulp_name "_ulp${_ulp}")
      endif()
      # check for custom split
      if(_splits_len GREATER 0)
        if(NOT _splits_len EQUAL _ids_len)
          message(SEND_ERROR "SPLITS and IDS must have the same length")
        endif()
        list(GET QAT_SPLITS ${_counter} _split)
        set(_split_name "_split${_split}")
        string(REGEX REPLACE ":" "_" _split_name ${_split_name})
        # parse split
        string(REGEX REPLACE ":" ";" _split ${_split})
        list(GET _split 0 _np)
        list(GET _split 1 _rank)
        set(_split_cmd "--np ${_np}" "--rank ${_rank}")
      endif()


      set(_testname "${QAT_COMMAND}_${_keyword}_id${_id}${_ulp_name}${_split_name}${_opt_name}" )
      message(DEBUG "_testname: ${_testname}")
      __add_test(${_testname}
        COMM ${QAT_COMMAND} [${${QAT_KEYWORD}}] -w NoTests ${_opt_cmd} --id=${_id} ${_ulp_cmd} ${_split_cmd}
        DIS ${QAT_DISABLED}
        STP ${_steps}
        PRF ${QAT_PERFONLY}
      )

      # loop index
      math(EXPR _counter "${_counter}+1")
    endforeach()
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
