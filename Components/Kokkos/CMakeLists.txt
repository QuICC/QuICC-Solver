include(ConfigureLibrary)
quicc_add_library(KokkosWrapper
  INCLUDE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  SRC
    Kokkos.cpp
  KOKKOSLIB
)
