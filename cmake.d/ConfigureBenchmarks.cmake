#
# Utility to add benchmark tests
#
# target
#     target implementation
# MODEL
#     name of the model
# ARCHIVEDIR
#     directory for storing the archive
# WORKDIR
#     working directory
# STARTFILES
#     list of files required to start
# TOOLS
#     list of tools
# VARIANTS
#     list of paths to edit in parameters.cfg
#     format: xmlpath:value
#
function(quicc_add_benchmark target)
  # parse inputs
  set(oneValueArgs MODEL ARCHIVEDIR WORKDIR TIMEOUT GITTAG MPIRANKS)
  set(multiValueArgs STARTFILES TOOLS VARIANTS)
  cmake_parse_arguments(QAB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_benchmark")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAB_MODEL: ${QAB_MODEL}")
  message(DEBUG "QAB_WORKDIR: ${QAB_WORKDIR}")
  message(DEBUG "QAB_ARCHIVEDIR: ${QAB_ARCHIVEDIR}")
  if(NOT QAB_MPIRANKS)
    set(QAB_MPIRANKS 4)
  endif()
  message(DEBUG "QAB_MPIRANKS: ${QAB_MPIRANKS}")

  if(NOT QAB_TIMEOUT)
    set(QAB_TIMEOUT 300)
  endif()
  message(DEBUG "QAB_TIMEOUT: ${QAB_TIMEOUT}")

  if(NOT QAB_STARTFILES)
    set(QAB_STARTFILES "parameters.cfg" "state_initial.hdf5")
  endif()
  message(DEBUG "QAB_STARTFILES: ${QAB_STARTFILES}")

  if(NOT QAB_TOOLS)
    set(QAB_TOOLS "validation_tools.py" "colorcodes.py")
  endif()
  message(DEBUG "QAB_TOOLS: ${QAB_TOOLS}")

  # default configs
  if(QUICC_MPI)
    set(_mpi_ranks ${QAB_MPIRANKS})
    set(_comm_algo "tubular")
  else()
    set(_mpi_ranks 1)
    set(_comm_algo "serial")
  endif()

  # Active variant filter
  set(_filterid "algorithm")
  # Check if there is an active variant or if we need to set the default
  set(_no_active_variant "True")
  foreach(_variant IN ITEMS ${QAB_VARIANTS})
    string(REGEX REPLACE ":" ";" _item "${_variant}")
    list(POP_BACK _item _value)
    string(REGEX REPLACE "/" ";" _item "${_item}")
    list(POP_BACK _item _name)
    list(FIND _filterid ${_name} _pos)
    if(_pos GREATER -1)
      set(_no_active_variant "False")
    endif()
  endforeach()

  message(DEBUG "QAB_VARIANTS: ${QAB_VARIANTS}")
  if(_no_active_variant)
    list(PREPEND QAB_VARIANTS "framework/parallel/algorithm:${_comm_algo}")
  endif()

  list(PREPEND QAB_VARIANTS "framework/parallel/grouper:transform")
  list(PREPEND QAB_VARIANTS "framework/parallel/cpus:${_mpi_ranks}")
  message(DEBUG "QAB_VARIANTS: ${QAB_VARIANTS}")

  foreach(_variant IN ITEMS ${QAB_VARIANTS})
    string(REGEX REPLACE ":" ";" _item "${_variant}")
    list(POP_BACK _item _value)
    string(REGEX REPLACE "/" ";" _item "${_item}")
    list(POP_BACK _item _name)
    list(FIND _filterid ${_name} _pos)
    if(_pos GREATER -1)
      string(APPEND _runid "_${_value}")
    endif()
  endforeach()
  message(DEBUG "_runid: ${_runid}")

  set(_exe "${QAB_MODEL}${target}Model")
  if(TARGET ${_exe})
    set(_bench "Benchmark${_exe}${_runid}")

    set(_refdir "${QAB_WORKDIR}/_refdata/${target}")
    set(_rundir "${QAB_WORKDIR}/_data/${target}${_runid}")
    message(VERBOSE "_rundir: ${_rundir}")
    set(_binsdir "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/Executables")
    message(VERBOSE "_binsdir: ${_binsdir}")
    set(_toolsdir "${PROJECT_SOURCE_DIR}/${QUICC_MODEL_PATH}/TestSuite")

    set(_args )
    foreach(_file IN LISTS QAB_STARTFILES)
      list(APPEND _args "COMMAND" ${CMAKE_COMMAND} -E copy
        "${_refdir}/${_file}"
        "${_rundir}/${_file}"
        )
    endforeach()

    foreach(_file IN LISTS QAB_TOOLS)
      list(APPEND _args "COMMAND" ${CMAKE_COMMAND} -E create_symlink
        "${_toolsdir}/${_file}"
        "${_rundir}/${_file}"
        )
    endforeach()
    message(DEBUG "custom target args: ${_args}")

    add_custom_target(${_bench} ALL
      COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/validate_benchmark_${target}.py"
        "${_rundir}/validate_benchmark.py"
      ${_args}
      )
    add_dependencies(${_bench} ${_exe})

    # Modify parameters.cfg for variants
    foreach(_variant IN ITEMS ${QAB_VARIANTS})
      string(REGEX REPLACE ":" ";" _vlist "${_variant}")
      list(GET _vlist 0 _path)
      list(GET _vlist 1 _value)
      add_custom_command(TARGET ${_bench} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove *.dat state0*.hdf5 *.gxl *.vtp
        COMMAND ${CMAKE_COMMAND} -E rename parameters.cfg parameters_orig.cfg
        COMMAND ${Python_EXECUTABLE} ${_toolsdir}/modify_xml.py -i parameters_orig.cfg -p ${_path} -v ${_value} -o parameters.cfg
        WORKING_DIRECTORY ${_rundir}
        )
    endforeach()

    # Fetch reference data
    include(FetchBenchmarkReference)
    quicc_fetch_benchmark_reference(
      ${_bench}
      MODEL ${QAB_MODEL}
      FILENAME "${target}.tar.gz"
      ARCHIVEDIR ${QAB_ARCHIVEDIR}
      DATADIR ${QAB_WORKDIR}
      GITTAG ${QAB_GITTAG}
    )

    set(_run "Run${_bench}")
    if(QUICC_MPI AND NOT QUICC_MPI_CI)
      # check which command is available
      foreach(_mpiexe IN ITEMS srun mpirun)
        message(VERBOSE "_mpiexe: ${_mpiexe}")
        find_program(mpiexe ${_mpiexe})
        if(mpiexe STREQUAL "mpiexe-NOTFOUND")
          message(VERBOSE "not found")
        else()
          message(VERBOSE "found")
          break()
        endif()
      endforeach()
      # check that we actually found something
      if(mpiexe STREQUAL "mpiexe-NOTFOUND")
        message(SEND_ERROR "could not find mpi executable.")
      endif()
      set(_test_param -n ${_mpi_ranks} "${_binsdir}/${_exe}")
      set(_command ${mpiexe} ${_test_param})
    else()
      set(_command ${_exe})
    endif()

    add_test(
      NAME ${_run}
      COMMAND ${_command}
      WORKING_DIRECTORY "${_rundir}"
      )
    set_tests_properties(${_run} PROPERTIES
      TIMEOUT ${QAB_TIMEOUT}
      )


    set(_validate "Validate${_bench}")
    add_test(
      NAME ${_validate}
      COMMAND "${Python_EXECUTABLE}" validate_benchmark.py
        -d "${_rundir}"
        -r "${_refdir}"
      WORKING_DIRECTORY "${_rundir}"
      )
    set_tests_properties(${_validate} PROPERTIES
      PASS_REGULAR_EXPRESSION "All benchmark validation tests passed!"
      TIMEOUT 60
      )

    set_tests_properties(${_validate} PROPERTIES DEPENDS "${_run}")
  else()
    message(WARNING "Tried to add benchmark but ${_exe} target is missing!")
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
