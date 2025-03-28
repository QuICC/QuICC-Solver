# Create component library
target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PRIVATE
    Cpu/Spec.cpp
    Cpu/Grid.cpp
)

if((TARGET QuICC::Fft::CuFft AND QUICC_USE_CUFFT) OR (TARGET QuICC::Fft::VkFft AND QUICC_USE_VKFFT))
  target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
    PRIVATE
      Cuda/Spec.cu
      Cuda/Grid.cu
  )
endif()

# Add public header
set(PHS
  Projector/SFG.hpp
  Integrator/GFS.hpp
  Spec.hpp
  Grid.hpp
  Tags.hpp
  Types.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
