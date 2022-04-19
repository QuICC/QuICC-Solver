option(QUICC_TESTSUITE_COMMON "Enable Common component testsuite?" OFF)
if(QUICC_TESTSUITE_COMMON)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_COMMON)
