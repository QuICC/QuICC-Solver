include(FetchReference)

#
# Utility to fetch reference data for tests
#
# target
#     name of the test to attach custom command to
# COMPONENT
#     component to which belong the tests
# FILENAME
#     filename of the archive
# ARCHIVEDIR
#     location where archive is stored
# DATADIR
#     location where data is extracted to ${DATADIR}/_refdata
function(quicc_fetch_test_reference target)
  # parse inputs
  set(oneValueArgs COMPONENT FILENAME ARCHIVEDIR DATADIR)
  set(multiValueArgs )
  cmake_parse_arguments(QAT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  message(DEBUG "quicc_fetch_test_reference")
  list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")
  message(DEBUG "target: ${target}")
  message(DEBUG "QAT_COMPONENT: ${QAT_COMPONENT}")
  message(DEBUG "QAT_FILENAME: ${QAT_FILENAME}")
  message(DEBUG "QAT_ARCHIVEDIR: ${QAT_ARCHIVEDIR}")
  message(DEBUG "QAT_DATADIR: ${QAT_DATADIR}")

  quicc_fetch_reference("${target}"
    NAME "${QAT_COMPONENT}"
    FILENAME "${QAT_FILENAME}"
    ARCHIVEDIR "${QAT_ARCHIVEDIR}"
    DATADIR "${QAT_DATADIR}"
    GITTAG "v1.0.28"
    GITURL "https://gitlab.ethz.ch/quicc/test-testdata"
    )

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()
