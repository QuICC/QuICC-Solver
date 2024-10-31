message(VERBOSE "Enabling TransformCoordinator tests:")

set(TestExe FrameworkTransformCoordinatorTests)

# Add target for all tests
add_executable(${TestExe} TransformCoordinatorTests.cpp)

target_link_libraries(${TestExe}
  ${QUICC_CURRENT_COMPONENT_TEST_LIB}
)

set(QUICC_WORK_DIR "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/TestSuite")

set_target_properties(
  ${TestExe}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${QUICC_WORK_DIR}"
)
list(APPEND CMAKE_MESSAGE_INDENT "${QUICC_CMAKE_INDENT}")

# Set different ULP for normal and MP implementation
include(SetUlp)

include(ConfigureTests)

set(_tests
  WLFl
  WLFm
  )
quicc_set_ulp(_ulp
  ULP 775 868 215000 107000 107000
  MPULP 11 11 13700 8544 8544
  )


foreach(_test ${_tests})
  foreach(_impl1d uniform triangular)
    foreach(_alg tubular single1d single2d)
      set(_dis False)
      if(_impl1d STREQUAL "triangular")
        # if(TARGET QuICC::Cuda)
          set(_dis True)
        # endif()
        # if(QUICC_USE_KOKKOS)
        #   set(_dis True)
        # endif()
      endif()
      quicc_add_test(${_test}
        COMMAND ${TestExe}
        KEYWORD Scheme
        ULPS ${_ulp}
        IDS 0 1 2 3 4
        OPTIONS algorithm:${_alg} truncation:${_impl1d}
        DISABLED ${_dis}
      )
    endforeach()
endforeach()
endforeach()

set(_tests
  SLFl
  SLFm
  )
quicc_set_ulp(_ulp
  ULP 155 110 53500 25090 25090
  MPULP 11 11 1850 1919 1919
  )
foreach(_test ${_tests})
  foreach(_alg tubular single1d single2d)
    quicc_add_test(${_test}
      COMMAND ${TestExe}
      KEYWORD Scheme
      ULPS ${_ulp}
      IDS 0 1 2 3 4
      OPTIONS algorithm:${_alg} truncation:uniform
    )
  endforeach()
endforeach()

set(_tests
  TFF
  )
quicc_set_ulp(_ulp
  ULP 11 11 1 1 1
  MPULP 11 11 11 11 11
  )
foreach(_test ${_tests})
  quicc_add_test(${_test}
    COMMAND ${TestExe}
    KEYWORD Scheme
    ULPS ${_ulp}
    IDS 0 1 2 3 4
    DISABLED ON
  )
endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
