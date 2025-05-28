set(PHS
  SparseOp/FiniteDiff/Operator.hpp
  SparseOp/FiniteDiff/OperatorWithMean.hpp
)

# Add public header
set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
  PROPERTIES
    PUBLIC_HEADER "${PHS}"
)
