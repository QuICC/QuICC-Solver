set(PHS
  DenseSM/Worland/GetSuperN.hpp
  DenseSM/Worland/Operator.hpp
  DenseSM/Worland/OperatorWithMean.hpp
  DenseSM/Worland/OperatorIN.hpp
  DenseSM/Worland/OperatorINWithMean.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
