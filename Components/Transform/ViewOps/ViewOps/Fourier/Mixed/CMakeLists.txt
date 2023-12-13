# Create component library
target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PRIVATE
    Cpu/Diff.cpp
)

if((TARGET QuICC::Fft::CuFft AND QUICC_USE_CUFFT) OR (TARGET QuICC::Fft::VkFft AND QUICC_USE_VKFFT))
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
    PRIVATE
      Cuda/Diff.cu
  )
endif()

# Add public header
set(PHS
  Mixed/Projector/D.hpp
  Mixed/Integrator/D.hpp
  Mixed/Diff.hpp
  Mixed/Types.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)