option(QUICC_TESTSUITE_POLYNOMIAL "Enable polynomial component testsuite?" OFF)
if(QUICC_TESTSUITE_POLYNOMIAL)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_POLYNOMIAL)
