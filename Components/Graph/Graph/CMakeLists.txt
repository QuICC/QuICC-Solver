# # Add public header
# set(PHS
#   Graph/MlirShims.hpp
# )

# set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
#   PROPERTIES
#     PUBLIC_HEADER "${PHS}"
# )

target_sources(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
    MlirShims.cpp
)