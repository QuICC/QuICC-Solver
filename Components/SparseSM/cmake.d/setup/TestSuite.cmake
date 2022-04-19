option(QUICC_TESTSUITE_SPARSESM "Enable SparseSM component testsuite?" OFF)
if(QUICC_TESTSUITE_SPARSESM)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_SPARSESM)
