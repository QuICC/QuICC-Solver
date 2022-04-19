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
  set(oneValueArgs MODEL FILENAME ARCHIVEDIR DATADIR)
  set(multiValueArgs )
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_fetch_benchmark_reference")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAT_MODEL: ${QAT_MODEL}")
  message(DEBUG "QAT_FILENAME: ${QAT_FILENAME}")
  message(DEBUG "QAT_ARCHIVEDIR: ${QAT_ARCHIVEDIR}")
  message(DEBUG "QAT_DATADIR: ${QAT_DATADIR}")

  quicc_fetch_reference("${target}"
    NAME "${QAT_MODEL}"
    FILENAME "${QAT_FILENAME}"
    ARCHIVEDIR "${QAT_ARCHIVEDIR}"
    DATADIR "${QAT_DATADIR}"
    GITTAG "v0.0.11"
    GITURL "https://gitlab.ethz.ch/quicc/test-benchmarks"
    )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
