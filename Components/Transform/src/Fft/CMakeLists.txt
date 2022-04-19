quicc_target_sources(${QUICC_CURRENT_COMPONENT_LIB} ${QUICC_CMAKE_SRC_VISIBILITY}
   IFftOperator.cpp
   Setup.cpp
   Tools.cpp
)

add_subdirectory(Backend)
add_subdirectory(Fourier)
add_subdirectory(Chebyshev)
add_subdirectory(Worland)
