option(QUICC_TESTSUITE_IO "Enable I/O component testsuite?" OFF)
if(QUICC_TESTSUITE_IO)
  add_subdirectory(TestSuite)
endif(QUICC_TESTSUITE_IO)
