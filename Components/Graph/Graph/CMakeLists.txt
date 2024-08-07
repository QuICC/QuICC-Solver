# Add public header
set(PHS
  ./Graph/OpsMap.hpp
  ./Graph/Tags.hpp
  ./Graph/Types.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

target_sources(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
    OpsMap.cpp
)

add_subdirectory(Maps)
add_subdirectory(Shims)
