set(PHS
  Environment/QuICCEnv.hpp
  Environment/Typedefs.hpp
  Environment/MpiTypes.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)

target_sources(${QUICC_CURRENT_COMPONENT_LIB}
    PRIVATE
        QuICCEnv.cpp
        IEnvironment.cpp
)

if(QUICC_USE_MPI)
    target_sources(${QUICC_CURRENT_COMPONENT_LIB}
        PRIVATE
            Mpi.cpp
    )
else()
    target_sources(${QUICC_CURRENT_COMPONENT_LIB}
        PRIVATE
            Serial.cpp
    )
endif()