set(tags
  FieldToRhs
  SolverHasBc
  SolverNoTau
  Stencil
)

include(RegisterTags)
quicc_register_tags(
  NAMESPACE "ModelOperatorBoundary"
  BASECLASS "IModelOperatorBoundary"
  TAGS ${tags}
  )
