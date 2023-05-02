# Create component library
target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} PRIVATE
  Projector/D.cpp
  Integrator/D.cpp
  Cpu/Diff.cpp
)

if(TARGET CUDA::cufft AND QUICC_USE_CUFFT)
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} PRIVATE
    Cuda/Diff.cu
  )
endif()

# Add public header
set(PHS
  Mixed/Projector/D.hpp
  Mixed/Integrator/D.hpp
  Mixed/Diff.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)