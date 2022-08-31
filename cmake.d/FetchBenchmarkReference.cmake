include(FetchReference)

#
# Utility to fetch reference data for benchmarks
#
# target
#     name of the test to attach custom command to
# MODEL
#     model to which belong the benchmark
# FILENAME
#     filename of the archive
# ARCHIVEDIR
#     location where archive is stored
# DATADIR
#     location where data is extracted to ${DATADIR}/_refdata
function(quicc_fetch_benchmark_reference target)
  # parse inputs
  set(oneValueArgs MODEL FILENAME ARCHIVEDIR DATADIR GITTAG)
  set(multiValueArgs )
  cmake_parse_arguments(QFB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_fetch_benchmark_reference")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QFB_MODEL: ${QFB_MODEL}")
  message(DEBUG "QFB_FILENAME: ${QFB_FILENAME}")
  message(DEBUG "QFB_ARCHIVEDIR: ${QFB_ARCHIVEDIR}")
  message(DEBUG "QFB_DATADIR: ${QFB_DATADIR}")

  if(NOT QFB_GITTAG)
    set(QFB_GITTAG "v0.0.19")
  endif()
  message(DEBUG "QFB_GITTAG: ${QFB_GITTAG}")

  quicc_fetch_reference("${target}"
    NAME "${QFB_MODEL}"
    FILENAME "${QFB_FILENAME}"
    ARCHIVEDIR "${QFB_ARCHIVEDIR}"
    DATADIR "${QFB_DATADIR}"
    GITTAG "${QFB_GITTAG}"
    GITURL "https://gitlab.ethz.ch/quicc/test-benchmarks"
    )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
