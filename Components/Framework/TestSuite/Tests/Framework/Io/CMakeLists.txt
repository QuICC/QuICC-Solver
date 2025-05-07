message(VERBOSE "Enabling Framework IO tests:")

set(TestExe FrameworkIoTests)

target_link_libraries(${QUICC_CURRENT_COMPONENT_TEST_LIB}
  QuICC::TestSuite
)

# Add target for all tests
add_executable(${TestExe} IoTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

if(QUICC_USE_MPI)
  # Add target for all MPI tests
  add_executable(${TestExe}Mpi IoTests.cpp)

  target_link_libraries(${TestExe}Mpi
    ${QUICC_CURRENT_COMPONENT_TEST_LIB}
  )
endif()

set(QUICC_WORK_DIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite")

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${QUICC_WORK_DIR}"
)
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Base dir for reference data archive
set(QUICC_REF_ARCHIVE_DIR "${CMAKE_BINARY_DIR}/_refdata")

# Fetch reference data
include(FetchTestReference)
quicc_fetch_test_reference(
  ${TestExe}
  COMPONENT Framework
  FILENAME "Io.tar.gz"
  ARCHIVEDIR ${QUICC_REF_ARCHIVE_DIR}
  DATADIR ${QUICC_WORK_DIR}
)

include(ConfigureTests)

set(_tests
  WLFl:Temperature:SphereNusseltWriter
  WLFl:Temperature:SphereScalarEnergyWriter
  WLFl:Temperature:SphereScalarLSpectrumWriter
  WLFl:Temperature:SphereScalarMSpectrumWriter
  WLFl:Temperature:SphereScalarMeanWriter
  WLFl:Temperature:SphereScalarModeSpectrumWriter
  WLFl:Temperature:SphereScalarNSpectrumWriter
  WLFl:Temperature:SphereScalarRSpectrumWriter
  WLFl:Velocity:SphereAngularMomentumWriter
  WLFl:Velocity:SphereTorPolEnergyWriter
  WLFl:Velocity:SphereTorPolLSpectrumWriter
  WLFl:Velocity:SphereTorPolMSpectrumWriter
  WLFl:Velocity:SphereTorPolNSpectrumWriter
  WLFl:Velocity:SphereTorPolRSpectrumWriter
  WLFl:Velocity:SphereTorPolModeSpectrumWriter
  WLFl:Velocity:SphereTorPolEnstrophyWriter
  WLFl:Velocity:SphereTorPolEnstrophyLSpectrumWriter
  WLFl:Velocity:SphereTorPolEnstrophyMSpectrumWriter
  WLFl:Magnetic:SphereDipolarityWriter
  WLFm:Temperature:SphereNusseltWriter
  WLFm:Temperature:SphereScalarEnergyWriter
  WLFm:Temperature:SphereScalarLSpectrumWriter
  WLFm:Temperature:SphereScalarMSpectrumWriter
  WLFm:Temperature:SphereScalarMeanWriter
  WLFm:Temperature:SphereScalarModeSpectrumWriter
  WLFm:Temperature:SphereScalarNSpectrumWriter
  WLFm:Temperature:SphereScalarRSpectrumWriter
  WLFm:Velocity:SphereAngularMomentumWriter
  WLFm:Velocity:SphereTorPolEnergyWriter
  WLFm:Velocity:SphereTorPolLSpectrumWriter
  WLFm:Velocity:SphereTorPolMSpectrumWriter
  WLFm:Velocity:SphereTorPolNSpectrumWriter
  WLFm:Velocity:SphereTorPolRSpectrumWriter
  WLFm:Velocity:SphereTorPolModeSpectrumWriter
  WLFm:Velocity:SphereTorPolEnstrophyWriter
  WLFm:Velocity:SphereTorPolEnstrophyLSpectrumWriter
  WLFm:Velocity:SphereTorPolEnstrophyMSpectrumWriter
  WLFm:Magnetic:SphereDipolarityWriter
  SLFl:Temperature:ShellNusseltWriter
  SLFl:Temperature:ShellScalarEnergyWriter
  SLFl:Temperature:ShellScalarLSpectrumWriter
  SLFl:Temperature:ShellScalarMSpectrumWriter
  SLFl:Temperature:ShellScalarModeSpectrumWriter
  SLFl:Temperature:ShellScalarRSpectrumWriter
  SLFl:Velocity:ShellTorPolEnergyWriter
  SLFl:Velocity:ShellTorPolLSpectrumWriter
  SLFl:Velocity:ShellTorPolMSpectrumWriter
  SLFl:Velocity:ShellTorPolModeSpectrumWriter
  SLFl:Velocity:ShellTorPolRSpectrumWriter
  SLFl:Velocity:ShellTorPolEnstrophyWriter
  SLFl:Velocity:ShellTorPolEnstrophyLSpectrumWriter
  SLFl:Velocity:ShellTorPolEnstrophyMSpectrumWriter
  SLFl:Magnetic:ShellDipolarityWriter
  SLFm:Temperature:ShellNusseltWriter
  SLFm:Temperature:ShellScalarEnergyWriter
  SLFm:Temperature:ShellScalarLSpectrumWriter
  SLFm:Temperature:ShellScalarMSpectrumWriter
  SLFm:Temperature:ShellScalarModeSpectrumWriter
  SLFm:Temperature:ShellScalarRSpectrumWriter
  SLFm:Velocity:ShellTorPolEnergyWriter
  SLFm:Velocity:ShellTorPolLSpectrumWriter
  SLFm:Velocity:ShellTorPolMSpectrumWriter
  SLFm:Velocity:ShellTorPolModeSpectrumWriter
  SLFm:Velocity:ShellTorPolRSpectrumWriter
  SLFm:Velocity:ShellTorPolEnstrophyWriter
  SLFm:Velocity:ShellTorPolEnstrophyLSpectrumWriter
  SLFm:Velocity:ShellTorPolEnstrophyMSpectrumWriter
  SLFm:Magnetic:ShellDipolarityWriter
)

# Create list of unique schemes
set(_schemes )
foreach(_test IN LISTS _tests)
  string(REPLACE ":" ";" _tags ${_test})
  list(GET _tags 0 _scheme)
  list(APPEND _schemes "${_scheme}")
endforeach()
list(REMOVE_DUPLICATES _schemes)

# Setup test data through custom command
set(_mod_cfg )
set(_copy_data )
foreach(_scheme IN LISTS _schemes)
  # Prepare serial tests
  list(APPEND _copy_data "COMMAND" "${CMAKE_COMMAND}" "-E" "remove" "_data/Framework/Io/${_scheme}/*")
  list(APPEND _copy_data "COMMAND" "${CMAKE_COMMAND}" "-E" "copy" "_refdata/Framework/Io/${_scheme}/parameters.cfg" "_data/Framework/Io/${_scheme}/")

  # Prepare MPI tests
  if(QUICC_USE_MPI)
    list(APPEND _copy_data "COMMAND" "${CMAKE_COMMAND}" "-E" "remove" "_data/Framework/Io/${_scheme}_mpi/*")
    list(APPEND _copy_data "COMMAND" "${CMAKE_COMMAND}" "-E" "copy_directory" "_data/Framework/Io/${_scheme}" "_data/Framework/Io/${_scheme}_mpi")
    list(APPEND _mod_cfg "COMMAND" "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/Models/TestSuite/modify_xml.py" "-i" "_data/Framework/Io/${_scheme}_mpi/parameters.cfg" "-p" "framework/parallel/algorithm" "-v" "tubular"  "-o" "_data/Framework/Io/${_scheme}_mpi/parameters_tmp.cfg")
    list(APPEND _mod_cfg "COMMAND" "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/Models/TestSuite/modify_xml.py" "-i" "_data/Framework/Io/${_scheme}_mpi/parameters_tmp.cfg" "-p" "framework/parallel/cpus" "-v" "4"  "-o" "_data/Framework/Io/${_scheme}_mpi/parameters.cfg")
    list(APPEND _mod_cfg "COMMAND" "${CMAKE_COMMAND}" "-E" "remove" "_data/Framework/Io/${_scheme}_mpi/parameters_tmp.cfg")
  endif()
endforeach()
add_custom_command(TARGET ${TestExe} POST_BUILD
  ${_copy_data}
  ${_mod_cfg}
  WORKING_DIRECTORY ${QUICC_WORK_DIR}
  )

# Create serial tests
foreach(_test IN LISTS _tests)
  string(REPLACE ":" ";" _tags ${_test})
  list(GET _tags 0 _scheme)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Scheme:Field:Writer
    ID 0
    ULP 500
    WORKDIR "${QUICC_WORK_DIR}/_data/Framework/Io/${_scheme}"
  )
endforeach()

# Create MPI tests
if(QUICC_USE_MPI)
  foreach(_test IN LISTS _tests)
    string(REPLACE ":" ";" _tags ${_test})
    list(GET _tags 0 _scheme)
    quicc_add_test(${_test}
      COMMAND ${TestExe}Mpi
      KEYWORD Scheme:Field:Writer
      ID 0
      ULP 500
      RANKS 4
      WORKDIR "${QUICC_WORK_DIR}/_data/Framework/Io/${_scheme}_mpi"
    )
  endforeach()
endif()

set(_tests
  WLFl:Velocity:SphereMaxAbsoluteFieldValueWriter
  WLFm:Velocity:SphereMaxAbsoluteFieldValueWriter
)

foreach(_test IN LISTS _tests)
  string(REPLACE ":" ";" _tags ${_test})
  list(GET _tags 0 _scheme)
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Scheme:Field:Writer
    ID 0
    ULP 500
    WORKDIR "${QUICC_WORK_DIR}/_data/Framework/Io/${_scheme}"
    DISABLED ON
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
