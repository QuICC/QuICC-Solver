option(QUICC_TESTSUITE_TRANSFORM "Enable transform component testsuite?" OFF)
if(QUICC_TESTSUITE_TRANSFORM)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_TRANSFORM)
