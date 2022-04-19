if(QUICC_NEEDS_CUDA)
  message(STATUS "***********************************************")
  message(STATUS "****************** CUDA setup *****************")
  message(STATUS "***********************************************")

  message(STATUS "--- Enabling CUDA")
  #enable_language(CUDA)

  set(QUICC_LIBRARIES_CUFFT "cufft" "cublas" "cusparse" PARENT_SCOPE)
endif()
