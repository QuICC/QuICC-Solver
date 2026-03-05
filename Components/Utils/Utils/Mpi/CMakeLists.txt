if(QUICC_USE_MPI)
    target_sources(${QUICC_CURRENT_COMPONENT_LIB}
        PRIVATE
            Utils.cpp
    )
endif()
