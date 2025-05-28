option(QUICC_TESTSUITE_FINITEDIFF "Enable Finite Differences component testsuite?" OFF)
if(QUICC_TESTSUITE_FINITEDIFF)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_FINITEDIFF)
