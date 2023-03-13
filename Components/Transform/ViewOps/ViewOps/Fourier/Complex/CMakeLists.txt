target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} PRIVATE
  Projector/D.cpp
  Integrator/D.cpp
  Cpu/Diff.cpp
  Cpu/Diff2D.cpp
  Cpu/Mean.cpp
  Cuda/Diff.cu
  Cuda/Diff2D.cu
  Cuda/Mean.cu
)

# Add public header
set(PHS
  Complex/Projector/D.hpp
  Complex/Integrator/D.hpp
  Complex/Diff.hpp
  Complex/Diff2D.hpp
  Complex/Mean.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
