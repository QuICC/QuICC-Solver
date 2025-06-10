set(PHS
  DenseOp/Worland/GetSuperN.hpp
  DenseOp/Worland/Operator.hpp
  DenseOp/Worland/OperatorWithMean.hpp
  DenseOp/Worland/OperatorIN.hpp
  DenseOp/Worland/OperatorINWithMean.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
