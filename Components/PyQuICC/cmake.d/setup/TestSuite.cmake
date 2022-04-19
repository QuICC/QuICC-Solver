option(QUICC_TESTSUITE_PYQUICC "Enable PyQuICC component testsuite?" OFF)
if(QUICC_TESTSUITE_PYQUICC)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_PYQUICC)
