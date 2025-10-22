option(QUICC_TESTSUITE_DENSESM "Enable DenseSM component testsuite?" OFF)
if(QUICC_TESTSUITE_DENSESM)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_DENSESM)
