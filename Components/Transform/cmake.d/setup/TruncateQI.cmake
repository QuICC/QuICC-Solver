# Add option to truncate quasi-inverse Worland operator
set( _truncateQI "QUICC_TRANSFORM_WORLAND_TRUNCATE_QI")
option(${_truncateQI} "Truncate Worland quasi-inverse operators" OFF)
if(${_truncateQI})
  target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC ${_truncateQI})
endif()

# Add option to truncate quasi-inverse Chebyshev operator
set( _truncateQI "QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI")
option(${_truncateQI} "Truncate Worland quasi-inverse operators" OFF)
if(${_truncateQI})
  target_compile_definitions(${QUICC_CURRENT_COMPONENT_LIB} PUBLIC ${_truncateQI})
endif()
