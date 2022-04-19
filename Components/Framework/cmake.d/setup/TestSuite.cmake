option(QUICC_TESTSUITE_FRAMEWORK "Enable framework component testsuite?" OFF)
if(QUICC_TESTSUITE_FRAMEWORK)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_FRAMEWORK)
