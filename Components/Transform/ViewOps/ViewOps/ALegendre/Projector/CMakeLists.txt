# Create component library
target_sources(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB} PRIVATE
  D.cpp
)

# Add public header
set(PHS
  Projector/D.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}_${QUICC_CURRENT_SUBCOMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
