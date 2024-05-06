# Add public header
set(PHS
  ../Graph/FourierMap.hpp
  ../Graph/OpsMap.hpp
  ../Graph/Tags.hpp
  ../Graph/Types.hpp
)

set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

target_sources(${QUICC_CURRENT_COMPONENT_LIB}
  PRIVATE
    OpsMap.cpp
    Map/FourierMap.cpp
    Map/ALegendreMap.cpp
    Map/WorlandMap.cpp
    Map/PointwiseMap.cpp
    Map/TransposeMap.cpp
)

add_subdirectory(Shims)
