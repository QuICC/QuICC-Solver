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
#
function(quicc_add_benchmark target)
  # parse inputs
  set(oneValueArgs MODEL ARCHIVEDIR WORKDIR)
  set(multiValueArgs STARTFILES TOOLS)
  cmake_parse_arguments(QAB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_add_benchmark")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "QAB_MODEL: ${QAB_MODEL}")
  message(DEBUG "QAB_WORKDIR: ${QAB_WORKDIR}")
  message(DEBUG "QAB_ARCHIVEDIR: ${QAB_ARCHIVEDIR}")

  if(NOT QAB_STARTFILES)
    set(QAB_STARTFILES "parameters.cfg" "state_initial.hdf5")
  endif()
  message(DEBUG "QAB_STARTFILES: ${QAB_STARTFILES}")

  if(NOT QAB_TOOLS)
    set(QAB_TOOLS "validation_tools.py" "colorcodes.py")
  endif()
  message(DEBUG "QAB_TOOLS: ${QAB_TOOLS}")

  set(_exe "${QAB_MODEL}${target}Model")
  set(_bench "Benchmark${_exe}")

  set(_refdir "${QAB_WORKDIR}/_refdata/${target}")
  set(_rundir "${QAB_WORKDIR}/_data/${target}")
  message(VERBOSE "_rundir: ${_rundir}")
  set(_binsdir "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_MODEL_DIR}/Executables")
  message(VERBOSE "_binsdir: ${_binsdir}")
  set(_toolsdir "${PROJECT_SOURCE_DIR}/${QUICC_MODEL_PATH}/TestSuite")

  set(_args )
  foreach(_file IN LISTS QAB_STARTFILES)
    list(APPEND _args "COMMAND" ${CMAKE_COMMAND} -E create_symlink
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

  if(QUICC_MPI)
    set(_mpi_ranks 4)
    add_custom_command(TARGET ${_bench} POST_BUILD
      COMMAND sed -i "'s/<cpus>1<\\/cpus>/<cpus>${_mpi_ranks}<\\/cpus>/g'" parameters.cfg
      WORKING_DIRECTORY ${_rundir}
      )
  endif()

  # Fetch reference data
  include(FetchBenchmarkReference)
  quicc_fetch_benchmark_reference(
    ${_bench}
    MODEL ${QAB_MODEL}
    FILENAME "${target}.tar.gz"
    ARCHIVEDIR ${QAB_ARCHIVEDIR}
    DATADIR ${QAB_WORKDIR}
  )

  set(_run "RunBenchmark${_exe}")
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
    TIMEOUT 300
    )


  set(_validate "ValidateBenchmark${_exe}")
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

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
